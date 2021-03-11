# Ensemble gene noise for the COVID19

library(edgeR)
library(gamlss)
library(gprofiler2)
library(parallel)

################################################################################
# gamlss function to estimate ensemble gene noise from RNA-seq counts
NBI_ens_noise <- function(counts, ensemble, geneL, min.genes.n = 10, cores = 4) {

  res <- mclapply(seq_along(ensemble), function(i) {
    cat(i, "\r")
    # From RNA-seq select genes constituting gene ensemble
    idx <- rownames(counts) %in% ensemble[[i]]
    dat <- counts[idx,]
    dat <- merge(geneL, dat, by = 'row.names')
    # ofs - ofset log(transcript length)
    ofs <- log(dat$length)
    dat <- dat[,-c(1:3)]

    if(length(ofs) >=  min.genes.n) {

      sink(file="/dev/null")
      # Gamlss estimation of ensemble gene noise from RNA counts using Negative Binomial model (NBI) and ofset - log(transcript length)
      # Try different Gamlss algorithms the Cole and Green (CG) and the Rigby and Stasinopoulos (RS) if CG fails
      res <- tryCatch(
        apply(dat, 2, function(x) exp(gamlss(x ~ offset(ofs),
          data = na.omit(data.frame(x=x, ofs=ofs)), family = NBI, method = CG(), n.cyc = 100)$sigma.coefficients) ),
      warning = function(w) NULL, error = function(e) NULL)

      if(is.null(res)) {
        res <- tryCatch(
          apply(dat, 2, function(x) exp(gamlss(x ~ offset(ofs),
            data = na.omit(data.frame(x=x, ofs=ofs)), family = NBI, method = RS(), n.cyc = 100)$sigma.coefficients) ),
        warning = function(w) NULL, error = function(e) NULL)
      }

      sink()

    } else { res <- NULL }
    res
  }, mc.cores = cores)

  names(res) <- names(ensemble)
  res <- do.call(rbind, res)

}

###############################################################################
# Read CORUM and KEGG gene ensembles annotations
setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
load("org.Hs.core_paths.Rdat")

names(corum_core) <- paste0("complex:", names(corum_core))
names(kegg_path) <- paste0("kegg:", names(kegg_path))

ENS <- c(corum_core, kegg_path)
ENS <- ENS[sapply(ENS, length) >= 5]

rm(corum_core, kegg_path, org.Hs)

ENSG_to_gene <- gconvert(rownames(trxL), organism = "hsapiens", target = c("ENTREZGENE"))[,c('input', 'target')]
ENSG_to_gene <- ENSG_to_gene[ENSG_to_gene$target != "nan",]
ENSG_to_gene <- ENSG_to_gene[!duplicated(ENSG_to_gene$target), ]


trxL <- merge(trxL, ENSG_to_gene, by.x = 'row.names', by.y = 'input')
rownames(trxL) <- trxL$target
trxL <- trxL[,-c(1,4)]

################################################################################
# Read COVID19 RNA-seq counts
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("COVID19_RNAseq.Rdat")

rm(samples, samples_pData)

################################################################################
# Calculate ensemble gene noise
COVID_ENS_noise <- NBI_ens_noise(DGE_RNAcounts$counts, ensemble = ENS, geneL = trxL, min.genes.n = 5, cores = 4)

################################################################################
# Calculate significant changes in ensemble gene noise for the COVID19 mild/severe patients
COVID_ENS_noise_sig <- mclapply(seq(nrow(COVID_ENS_noise)), function(i) {

  cat(i, "\r")

  y = c(COVID_ENS_noise[i,])
  x = DGE_RNAcounts$samples$stage
  Sex = DGE_RNAcounts$samples$Sex
  age = DGE_RNAcounts$samples$age

  dat <- data.frame(y = y, x = x, Sex = Sex, age = age)
  dat$age[is.na(dat$age)] <- median(dat$age[!is.na(dat$age)])
  # boxplot(y ~ x, dat)

  sink("/dev/null")
  # gamlss model accounting for Sex and age as random effects and COVID19 mild vs sever as fixed effect factor
  # Try different Gamlss algorithms the Cole and Green (CG) and the Rigby and Stasinopoulos (RS) if CG fails
  m0 <-  tryCatch(
    gamlss::gamlss(y ~ x +  re(random = ~ 1|Sex) + re(random = ~ 1|pb(age)), data = dat, n.cyc = 1000),
    error = function(e) NULL, warning = function(w) NULL)

  if(is.null(m0)) {
    m0 <-  tryCatch(
      gamlss::gamlss(y ~ x +  re(random = ~ 1|Sex) + re(random = ~ 1|pb(age)), data = dat, method = CG(), n.cyc = 1000),
      error = function(e) NULL, warning = function(w) NULL)
  }
  # If both RS and CG fail drop random effects
  if(is.null(m0)) {
    m0 <-  tryCatch(
      gamlss::gamlss(y ~ x, data = dat, n.cyc = 1000),
      error = function(e) NULL, warning = function(w) NULL)
  }

  if(!is.null(m0)) {
    m0s <-  summary(m0, method = "vcov", save = TRUE)
    res <-  t(as.data.frame(m0s$coef.table[2,3:4]))
    colnames(res) <- c("t", "p")
    rownames(res) <- NULL
  } else {
    res <- data.frame(t = NA, p = NA)
  }

  sink()

  res

}, mc.cores = 6)

COVID_ENS_noise_sig <- do.call(rbind, COVID_ENS_noise_sig)
rownames(COVID_ENS_noise_sig) <- rownames(COVID_ENS_noise)
COVID_ENS_noise_sig <- as.data.frame(COVID_ENS_noise_sig)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("COVID19_Ensemble_noise.RDat")

################################################################################
