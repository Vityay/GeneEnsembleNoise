# Figure S1

library(deming)
library(quantreg)
library(mvoutlier)
library(gamlss)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("GSE65682_series_all_mortality.Rdat")

# pool CAP and Sepsis patients
pheno$d2 <- as.character(pheno$d)
pheno$d2[pheno$d2 == "other"] <- "cap"

# Let's further trim lowly expressed / noisy genes
out <- sign2(affyData, explvar = 0.888265, qcrit = 1-1e-8)
affyData_out2 <- affyData[out$wfinal01 == 0, ]
# Compare
# plot(density(unlist(affyData[out$wfinal01 == 1, ])), xlim = c(1,14))
# lines(density(unlist(affyData[out$wfinal01 == 0, ])), col=2)

delta <- apply(affyData_out2, 1, function(x) {
  # calculate value of t test for each gene (CAP/Sepsis vs healthy)
  df <- data.frame( y = as.numeric(x), x =  pheno$d2)
  res <- t.test(y ~ x, df)
  # calculate difference for each gene (CAP/Sepsis vs healthy)
  d <- as.numeric(-diff(res$estimate))
  # estimate difference in variances for each gene (CAP/Sepsis vs healthy)
  m0 <- gamlss(y ~ 0 + x, sigma.fo = ~ 0 + x, data=df, family=NO)
  dv <- -diff(exp(m0$sigma.coefficients)^2)
  res <- c(d = d, p = res$p.value, t = res$statistic, dv = dv, v = var(df$y), mu = mean(df$y))
  res
})
delta <- data.frame(t(delta))
delta$padj <- p.adjust(delta$p, "bonferroni")

################################################################################

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("FigureS1.pdf", width=11, height=8.5)

par(mfrow=c(2,3))
dfd <-  data.frame(v = delta$v, d = abs(delta$d))
R_dfd <- cor.test(dfd$d, dfd$v, method="spearman")
plot(d ~ v, data=dfd, las=1, xlab="variance", ylab="difference", xlim=c(0.05,10), ylim=c(0.01,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfd)[ sign2(log10(dfd), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dft <- data.frame(v = delta$v, t = abs(delta$t))
R_dft <- cor.test(dft$t, dft$v, method="spearman")
plot(t ~ v, data=dft, las=1, xlab="variance", ylab="t value", xlim=c(0.05,10), ylim=c(0.05,100), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(t ~ v,, data =  log10(dft)[ sign2(log10(dft), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dft$estimate, 3)), paste0("p = ", round(R_dft$p.value, 5)) ) )

dfp <- na.omit(data.frame(v = delta$v, p = -log10(delta$p)))
dfp <- dfp[dfp$p >0,]
R_dfp <- cor.test(dfp$p, dfp$v, method="spearman")
plot(p ~ v, data=dfp, las=1, xlab="variance", ylab="-log10(p)", xlim=c(0.05,10), ylim=c(0.005, 500), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(p ~ v,, data =  log10(dfp)[ sign2(log10(dfp), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfp$estimate, 3)), paste0("p = ", round(R_dfp$p.value, 5)) ) )


################################################################################

dfd <-  data.frame(mu = delta$mu, d = abs(delta$d))
R_dfd <- cor.test(dfd$d, dfd$mu, method="spearman")
plot(d ~ mu, data=dfd, las=1, xlab="mean", ylab="difference", xlim=c(2,15), ylim=c(0.01,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ mu, data = log10(dfd), 0.5), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dft <- data.frame(mu = delta$mu, t = abs(delta$t))
R_dft <- cor.test(dft$t, dft$mu, method="spearman")
plot(t ~ mu, data=dft, las=1, xlab="mean", ylab="t value", xlim=c(2,15), ylim=c(0.05,100), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(t ~ mu, data = log10(dft), 0.5), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dft$estimate, 3)), paste0("p = ", round(R_dft$p.value, 5)) ) )

dfp <- na.omit(data.frame(mu = delta$mu, p = -log10(delta$p)))
dfp <- dfp[dfp$p >0,]
R_dfp <- cor.test(dfp$p, dfp$mu, method="spearman")
plot(p ~ mu, data=dfp, las=1, xlab="mean", ylab="-log10(p)", xlim=c(2,15), ylim=c(0.005, 500), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(p ~ mu, data = log10(dfp), 0.5), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfp$estimate, 3)), paste0("p = ", round(R_dfp$p.value, 5)) ) )

dev.off()
