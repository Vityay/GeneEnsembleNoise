# Gamlss estimation of gene noise for each gene for COVID19

library(edgeR)
library(gamlss)
library(gtools)


setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("COVID19_RNAseq.Rdat")


################################################################################
# Functions for GAMLSS analysis
# gNBI.RS - the Rigby and Stasinopoulos algorithm for GAMLSS
gNBI.RS = function(dat, f1, f2){
  res = tryCatch(
    gamlss::gamlss(f1, sigma.fo = f2, data = dat, family = NBI(), method = RS(), n.cyc = 500),
  warning = function(w) NULL, error = function(e) NULL)
  res
}

# gNBI.CG -  the Cole and Green algorithm for GAMLSS
gNBI.CG = function(dat, f1, f2){
  res = tryCatch(
    gamlss::gamlss(f1, sigma.fo = f2, data = dat, family = NBI(), method = CG(), n.cyc = 500),
  warning = function(w) NULL, error = function(e) NULL)
  res
}

# Gamlss estimation of gene noise for each gene for COVID19
gene_noise_COVID <- mclapply(seq(nrow(DGE_RNAcounts$counts)), function(i) {

    cat(i, "\r")

    # COVID - data.frame
    # y - RNA counts
    # x - fixed effect factor (COVID19 mild/severe)
    # ofs - offset log(library size x normalization factor), see edgeR
    # Sex - gender
    # age - age
    # rnd - Sex:age interaction was used as random effect

    COVID = data.frame(y = DGE_RNAcounts$counts[i,],
      x = factor(DGE_RNAcounts$samples$stage),
      ofs = DGE_RNAcounts$samples$offset,
      Sex = DGE_RNAcounts$samples$Sex,
      age = DGE_RNAcounts$samples$age
    )
    COVID$age[is.na(COVID$age)] = median(COVID$age[!is.na(COVID$age)])
    COVID$age = quantcut(COVID$age, q=10)
    COVID$rnd = paste0(COVID$Sex, ":", COVID$age)

    # f1 - gamlss formula for mu
    f1 = formula(y ~ x + random(factor(rnd)) + offset(ofs) )
    # f2 - gamlss formula for sigma
    f2 = formula(y ~ x + random(factor(rnd)))

    # gamlss model
    # m0 - model for mu and sigma with x fixed factor effects: COVID stage and random factor effects rnd: Sex:age
    sink("/dev/null")

    m0 = gNBI.RS(COVID, f1, f2)

    if(is.null(m0)) { m0 = gNBI.CG(COVID, f1, f2)  }

    sink()

    # results:
    # mu_mild - mean gene expression for mild COVID19
    # mu_severe - mean gene expression for severe COVID19
    # bcv_mild - biological coefficient of variation of gene expression for mild COVID19
    # bcv_severe - biological coefficient of variation of gene expression for severe COVID19
    # p_mu - significance of COVID19 effect (mild vs severe) on mean gene expression
    # p_bcv - significance of COVID19 effect (mild vs severe) on gene expression noise (bcv - biological coefficient of variation)

    res = data.frame(mu_mild = NA, mu_severe = NA, bcv_mild = NA, bcv_severe = NA, p_mu = NA, p_bcv=NA)

    if(!is.null(m0)) {
      sink("/dev/null")
        m0s = summary(m0, "qr", save = T)
        mu = exp(cumsum(m0s$mu.coef.table[,1]))*1e06
        bcv = sqrt(exp(cumsum(m0s$sigma.coef.table[,1])))
        p_mu = m0s$mu.coef.table[2,4]
        p_bcv = m0s$sigma.coef.table[2,4]
        res = data.frame(mu_mild = mu[[1]], mu_severe = mu[[2]],
            bcv_mild = bcv[[1]], bcv_severe = bcv[[2]],
            p_mu = p_mu, p_bcv=p_bcv)
      sink()
    }

    res

}, mc.cores = 4 )

gene_noise_COVID <- do.call(rbind, gene_noise_COVID)
rm(f1, f2, gNBI.CG, gNBI.RS)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("COVID19_gene_noise_Figure6A.Rdat")
