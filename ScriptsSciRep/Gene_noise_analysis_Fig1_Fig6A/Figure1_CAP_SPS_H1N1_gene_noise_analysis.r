library(gamlss)
library(gtools)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("GSE65682_series_all_mortality.Rdat")

################################################################################
# Gamlss estimation of gene noise for each gene for CAP

# f - factor:
# ctl - control;
# cap:0 - CAP survived; cap:1 - CAP decieased;
# other:0 - Sepsis survived; other:1 - Sepsis decieased;
f <- gsub(":NA", "", paste(pheno$d, pheno$mort, sep=":"))
f <- relevel(factor(f), ref="ctl")

gene_noise_CAP <- mclapply(seq(nrow(affyData_out)), function(i) {
  cat(i, "\n")

  # df - data.frame
  # y - log2 hybridizations
  # x - factor (ctl - control; cap:0 - CAP survived; cap:1 - CAP decieased)
  # x2 - factor (ctl - control; cap - CAP)
  df <- data.frame(y = as.numeric(affyData_out[i,]), x = f, age = pheno$age)
  df <- subset(df, x != "other:0" & x != "other:1")
  df$x <- droplevels(df$x)
  df$x2 <- factor(substring(as.character(df$x), 1,3), levels=c("ctl", "cap"))
  df$age <- quantcut(df$age, q=10)

  sink("/dev/null")
   # df$yr <- resid(gamlss(y ~ age, sigma.fo = ~ age, data = df, n.cyc=1000))
   # df$yr <- df$y

   # gamlss models
   # m0 - model for mu and sigma with x factor effects: ctl - control; cap:0 - CAP survived; cap:1 - CAP decieased
   # m1 - model for mu with x factor effects: ctl - control; cap:0 - CAP survived; cap:1 - CAP decieased
   # m2 - model for sigma with x factor effects: ctl - control; cap:0 - CAP survived; cap:1 - CAP decieased

   m0 <- tryCatch( gamlss(y ~ 0 + x , sigma.fo = ~ 0 + x , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m1 <- tryCatch( gamlss(y ~ 0 + x , sigma.fo = ~ 1, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m2 <- tryCatch( gamlss(y ~ 1, sigma.fo = ~ 0 + x , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)

   # m0a - model for mu and sigma with x2 factor effects: ctl - control; cap - CAP
   # m1a - model for mu with x2 factor effects: ctl - control; cap - CAP
   # m2a - model for sigma with x2 factor effects: ctl - control; cap - CAP

   m0a <- tryCatch( gamlss(y ~ 0 + x2 , sigma.fo = ~ 0 + x2 , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m1a <- tryCatch( gamlss(y ~ 0 + x2 , sigma.fo = ~ 1, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m2a <- tryCatch( gamlss(y ~ 1, sigma.fo = ~ 0 + x2 , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
  sink()

  if(any(sapply(list(m0,m1,m2), is.null))) {
    res1 <- data.frame(mu_ctl = NA, mu_cap0 = NA, mu_cap1 = NA, v_ctl = NA, v_cap0 = NA, v_cap1 = NA, V_ave = NA, p_mu = NA, p_sigma=NA)
  } else {
      # LR test for the significance of CAP effect (ctl, cap:0, cap:1) on mean gene expression
      p_mu <- pchisq(2*as.numeric(logLik(m0)-logLik(m2)), 2, lower.tail=F)
      # LR test for the significance of CAP effect on gene expression noise (sigma - standard deviation)
      p_sigma <- pchisq(2*as.numeric(logLik(m0)-logLik(m1)), 2, lower.tail=F)
      # Estimation of mean gene expression
      mu <- mean(df$y) + m0$mu.coefficients[1:3]
      # Estimation of gene expression noise (sigma^2 - variance)
      sigma <- exp(m0$sigma.coefficients[1:3])^2
      # results:
      # mu_ctl - mean gene expression for control
      # mu_cap0 - mean gene expression for cap0 (CAP survived)
      # mu_cap1 - mean gene expression for cap1 (CAP decieased)
      # v_ctl - gene expression variance for control
      # v_cap0 - gene expression variance for cap0 (CAP survived)
      # v_cap1 - gene expression variance for cap1 (CAP decieased)
      # V_ave - pooled gene expression variance
      # p_mu - significance of CAP effect on mean gene expression
      # p_sigma - significance of CAP effect on gene expression noise (sigma - standard deviation)
      res1 <- data.frame(mu_ctl = mu[[1]], mu_cap0 = mu[[2]], mu_cap1 = mu[[3]],
        v_ctl = sigma[[1]], v_cap0 = sigma[[2]], v_cap1 = sigma[[3]], V_ave = exp(m1$sigma.coefficients[[1]])^2,
        p_mu = p_mu, p_sigma=p_sigma)
  }

  if(any(sapply(list(m0a,m1a,m2a), is.null))) {
    res2 <- data.frame(mu2_ctl = NA, mu2_cap = NA, v2_ctl = NA, v2_cap = NA, p2_mu = NA, p2_sigma=NA)
  } else {
      # LR test for the significance of CAP effect (ctl, cap) on mean gene expression and sigma
      p2_mu  <- pchisq(2*as.numeric(logLik(m0a)-logLik(m2a)), 1, lower.tail=F)
      p2_sigma <- pchisq(2*as.numeric(logLik(m0a)-logLik(m1a)), 1, lower.tail=F)
      mu2 <- mean(df$y) + m0a$mu.coefficients[1:3]
      sigma2 <- exp(m0a$sigma.coefficients[1:3])^2
      # results:
      # mu2_ctl - mean gene expression for control
      # mu2_cap - mean gene expression for CAP (survived/decieased)
      # v2_ctl - gene expression variance for control
      # v2_cap - gene expression variance for CAP (survived/decieased)
      # p2_mu - significance of ctl vs CAP effect on mean gene expression
      # p2_sigma - significance of ctl vs CAP effect on gene expression noise (sigma - standard deviation)
      res2 <- data.frame(mu2_ctl = mu2[[1]], mu2_cap = mu2[[2]],
        v2_ctl = sigma2[[1]], v2_cap = sigma2[[2]],
        p2_mu = p2_mu, p2_sigma=p2_sigma)
  }

  res <- cbind(res1, res2)
}, mc.cores=6)
gene_noise_CAP <- do.call(rbind, gene_noise_CAP)
gene_noise_CAP$padj_sigma <- p.adjust(gene_noise_CAP$p_sigma, "bonferroni")
gene_noise_CAP$padj_mu <- p.adjust(gene_noise_CAP$p_mu, "bonferroni")
rownames(gene_noise_CAP) <- rownames(affyData_out)

################################################################################
# Gamlss estimation of gene noise for each gene for Sepsis (SPS) (same as above)

gene_noise_SPS <- mclapply(seq(nrow(affyData_out)), function(i) {
  cat(i, "\n")
  df <- data.frame(y = as.numeric(affyData_out[i,]), x = f, age = pheno$age)
  df <- subset(df, x != "cap:0" & x != "cap:1")
  df$x <- droplevels(df$x)
  df$x2 <- factor(substring(as.character(df$x), 1,3), levels=c("ctl", "oth"))
  df$age <- quantcut(df$age, q=10)

  sink("/dev/null")
   # df$yr <- resid(gamlss(y ~ age, sigma.fo = ~ age, data = df, n.cyc=1000))
   # df$yr <- df$y

   m0 <- tryCatch( gamlss(y ~ 0 + x , sigma.fo = ~ 0 + x , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m1 <- tryCatch( gamlss(y ~ 0 + x , sigma.fo = ~ 1, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m2 <- tryCatch( gamlss(y ~ 1, sigma.fo = ~ 0 + x , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)

   m0a <- tryCatch( gamlss(y ~ 0 + x2 , sigma.fo = ~ 0 + x2 , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m1a <- tryCatch( gamlss(y ~ 0 + x2 , sigma.fo = ~ 1, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m2a <- tryCatch( gamlss(y ~ 1, sigma.fo = ~ 0 + x2 , data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
  sink()

  if(any(sapply(list(m0,m1,m2), is.null))) {
    res1 <- data.frame(mu_ctl = NA, mu_sps0 = NA, mu_sps1 = NA, v_ctl = NA, v_sps0 = NA, v_sps1 = NA, V_ave = NA, p_mu = NA, p_sigma=NA)
  } else {
      p_mu <- pchisq(2*as.numeric(logLik(m0)-logLik(m2)), 2, lower.tail=F)
      p_sigma <- pchisq(2*as.numeric(logLik(m0)-logLik(m1)), 2, lower.tail=F)
      mu <- mean(df$y) + m0$mu.coefficients[1:3]
      sigma <- exp(m0$sigma.coefficients[1:3])^2
      res1 <- data.frame(mu_ctl = mu[[1]], mu_sps0 = mu[[2]], mu_sps1 = mu[[3]],
        v_ctl = sigma[[1]], v_sps0 = sigma[[2]], v_sps1 = sigma[[3]], V_ave = exp(m1$sigma.coefficients[[1]])^2,
        p_mu = p_mu, p_sigma=p_sigma)
  }

  if(any(sapply(list(m0a,m1a,m2a), is.null))) {
    res2 <- data.frame(mu2_ctl = NA, mu2_sps = NA, v2_ctl = NA, v2_sps = NA, p2_mu = NA, p2_sigma=NA)
  } else {
      p2_mu  <- pchisq(2*as.numeric(logLik(m0a)-logLik(m2a)), 1, lower.tail=F)
      p2_sigma <- pchisq(2*as.numeric(logLik(m0a)-logLik(m1a)), 1, lower.tail=F)
      mu2 <- mean(df$y) + m0a$mu.coefficients[1:3]
      sigma2 <- exp(m0a$sigma.coefficients[1:3])^2
      res2 <- data.frame(mu2_ctl = mu2[[1]], mu2_sps = mu2[[2]],
        v2_ctl = sigma2[[1]], v2_sps = sigma2[[2]],
        p2_mu = p2_mu, p2_sigma=p2_sigma)
  }

  res <- cbind(res1, res2)
}, mc.cores=6)
gene_noise_SPS <- do.call(rbind, gene_noise_SPS)
gene_noise_SPS$padj_sigma <- p.adjust(gene_noise_SPS$p_sigma, "bonferroni")
gene_noise_SPS$padj_mu <- p.adjust(gene_noise_SPS$p_mu, "bonferroni")
rownames(gene_noise_SPS) <- rownames(affyData_out)

################################################################################
# Gamlss estimation of gene noise for each gene for H1N1 (same as above)

load("eset_H1N1.Rdat")
H1N1_stage <- samples_H1N1$'disease phase:ch1'

gene_noise_H1N <- mclapply(seq(nrow(eset_H1N1)), function(i) {
  cat(i, "\n")
  df <- data.frame(y = as.numeric(eset_H1N1[i,]), x = factor(H1N1_stage))
  df$x2 <- factor(gsub("late ", "", gsub("early ", "", as.character(df$x))), levels=c("control", "period"))

  sink("/dev/null")

   m0 <- tryCatch( gamlss(y ~ 0 + x, sigma.fo = ~ 0 + x, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m1 <- tryCatch( gamlss(y ~ 0 + x, sigma.fo = ~ 1, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m2 <- tryCatch( gamlss(y ~ 1, sigma.fo = ~ 0 + x, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)

   m0a <- tryCatch( gamlss(y ~ 0 + x2, sigma.fo = ~ 0 + x2, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m1a <- tryCatch( gamlss(y ~ 0 + x2, sigma.fo = ~ 1, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
   m2a <- tryCatch( gamlss(y ~ 1, sigma.fo = ~ 0 + x2, data = df, n.cyc=1000), warning = function(w) NULL, error = function(e) NULL)
  sink()

  if(any(sapply(list(m0,m1,m2), is.null))) {
    res1 <- data.frame(mu_ctl = NA, mu_h1n0 = NA, mu_h1n1 = NA, v_ctl = NA, v_h1n0 = NA, v_h1n1 = NA, V_ave = NA, p_mu = NA, p_sigma=NA)
  } else {
      p_mu <- pchisq(2*as.numeric(logLik(m0)-logLik(m2)), m0$df.fit-m2$df.fit, lower.tail=F)
      p_sigma <- pchisq(2*as.numeric(logLik(m0)-logLik(m1)), m0$df.fit-m1$df.fit, lower.tail=F)
      mu <- mean(df$y) + m0$mu.coefficients
      sigma <- exp(m0$sigma.coefficients)^2
      # results:
      # mu_ctl - mean gene expression for control
      # mu_h1n0 - mean gene expression for H1N1 early
      # mu_h1n1 - mean gene expression for H1N1 late
      # v_ctl - gene expression variance for control
      # v_h1n0 - gene expression variance for H1N1 early
      # v_h1n1 - gene expression variance for H1N1 late
      # V_ave - pooled gene expression variance
      # p_mu - significance of H1N1 effect on mean gene expression
      # p_sigma - significance of H1N1 effect on gene expression noise (sigma - standard deviation)
      res1 <- data.frame(mu_ctl = mu[[1]], mu_h1n0 = mu[[2]], mu_h1n1 = mu[[3]],
        v_ctl = sigma[[1]], v_h1n0 = sigma[[2]], v_h1n1 = sigma[[3]], V_ave = exp(m1$sigma.coefficients)^2,
        p_mu = p_mu, p_sigma=p_sigma)
  }

  if(any(sapply(list(m0a,m1a,m2a), is.null))) {
    res2 <- data.frame(mu2_ctl = NA, mu2_h1n = NA, v2_ctl = NA, v2_h1n = NA, p2_mu = NA, p2_sigma=NA)
  } else {
      p2_mu  <- pchisq(2*as.numeric(logLik(m0a)-logLik(m2a)), m0a$df.fit-m2a$df.fit, lower.tail=F)
      p2_sigma <- pchisq(2*as.numeric(logLik(m0a)-logLik(m1a)), m0a$df.fit-m1a$df.fit, lower.tail=F)
      mu2 <- mean(df$y) + m0a$mu.coefficients
      sigma2 <- exp(m0a$sigma.coefficients)^2
      # results:
      # mu2_ctl - mean gene expression for control
      # mu2_h1n - mean gene expression for H1N1 (early/late)
      # v2_ctl - gene expression variance for control
      # v2_h1n - gene expression variance for H1N1 (early/late)
      # p2_mu - significance of ctl vs H1N1 effect on mean gene expression
      # p2_sigma - significance of ctl vs H1N1 effect on gene expression noise (sigma - standard deviation)
      res2 <- data.frame(mu2_ctl = mu2[[1]], mu2_h1n = mu2[[2]],
        v2_ctl = sigma2[[1]], v2_h1n = sigma2[[2]],
        p2_mu = p2_mu, p2_sigma=p2_sigma)
  }

  res <- cbind(res1, res2)
}, mc.cores=6)
gene_noise_H1N <- do.call(rbind, gene_noise_H1N)
gene_noise_H1N$padj_sigma <- p.adjust(gene_noise_H1N$p_sigma, "BH")
gene_noise_H1N$padj_mu <- p.adjust(gene_noise_H1N$p_mu, "BH")
rownames(gene_noise_H1N) <- rownames(eset_H1N1)

rm(f, affyData, affyData_out,eset_H1N1)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("CAP_Sepsis_H1N1_gene_noise_Figure1.Rdat")
