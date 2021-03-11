library(deming)
library(quantreg)
library(mvoutlier)
library(gamlss)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
load("CAP_Sepsis_Ensemble_noise.RDat")

pheno$d2 <- as.character(pheno$d)
pheno$d2[pheno$d2 == "other"] <- "cap"


GN_var <- apply(affyDat_noise_paths, 2, var)
GN_mu <- apply(affyDat_noise_paths, 2, mean)

GN_diff <- apply(affyDat_noise_paths, 2, function(x) {
  df <- data.frame(v = as.numeric(x), d = pheno$d2)
  res <- -diff(aggregate(v ~ d, df, mean)$v)
})

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("FigureS4.pdf", paper="a4r", width=7.3, height=4.25)

par(mfrow=c(1,2))
dfd <-  data.frame(v = GN_var, d = abs(GN_diff))
R_dfd <- cor.test(dfd$d, dfd$v, method="spearman")
plot(d ~ v, data=dfd, las=1, xlab="variance", ylab="difference",
     xlim=c(0.01,10), ylim=c(0.01,10), pch=20, col="#00000088", log="xy", cex=1)
abline(rq(d ~ v,, data =  log10(dfd)), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dfd <-  data.frame(mu = GN_mu, d = abs(GN_diff))
R_dfd <- cor.test(dfd$d, dfd$mu, method="spearman")
plot(d ~ mu, data=dfd, las=1, xlab="mean", ylab="difference",
     xlim=c(1,15), ylim=c(0.01,10), pch=20, col="#00000088", log="xy", cex=1)
abline(rq(d ~ mu, data = log10(dfd), 0.5), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dev.off()
