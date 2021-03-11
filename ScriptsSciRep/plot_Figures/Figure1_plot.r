# library(edgeR)
# library(gamlss.mx)
# library(parallel)
# library(tidyverse)
# library(SuppDists)
# library(quantreg)
# library(caret)
# library(MASS)
# library(viridis)
# library(gplots)
library(reshape2)


setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("CAP_Sepsis_H1N1_gene_noise_Figure1.Rdat")

################################################################################

# Mean gene expression noise (variance) +/- std.dev for CAP (see CAP_Sepsis_H1N1_gene_noise_Figure1.r)
df1_cap <- melt(sqrt(gene_noise_CAP[,c(4:6)]))
df1_cap_mu <- aggregate(value ~ variable, df1_cap, mean)$value
df1_cap_upp_sd <- aggregate(value ~ variable, data=df1_cap, function(x) mean(x) + sd(x))$value
df1_cap_lwr_sd <- aggregate(value ~ variable, data=df1_cap, function(x) mean(x) - sd(x))$value

# Mean gene expression noise (variance) +/- std.dev for genes for which variance changes significantly in CAP
# (see CAP_Sepsis_H1N1_gene_noise_Figure1.r)
df2_cap <- melt(sqrt(subset(gene_noise_CAP, padj_sigma <= 0.05)[,c(4:6)]))
df2_cap_mu <- aggregate(value ~ variable, df2_cap, mean)$value
df2_cap_upp_sd <- aggregate(value ~ variable, data=df2_cap, function(x) mean(x) + sd(x))$value
df2_cap_lwr_sd <- aggregate(value ~ variable, data=df2_cap, function(x) mean(x) - sd(x))$value

# Mean gene expression noise (variance) +/- std.dev for Sepsis (see CAP_Sepsis_H1N1_gene_noise_Figure1.r)
df1_sps <- melt(sqrt(gene_noise_SPS[,c(4:6)]))
df1_sps_mu <- aggregate(value ~ variable, df1_sps, mean)$value
df1_sps_upp_sd <- aggregate(value ~ variable, data=df1_sps, function(x) mean(x) + sd(x))$value
df1_sps_lwr_sd <- aggregate(value ~ variable, data=df1_sps, function(x) mean(x) - sd(x))$value

# Mean gene expression noise (variance) +/- std.dev for genes for which variance changes significantly in Sepsis
# (see CAP_Sepsis_H1N1_gene_noise_Figure1.r)
df2_sps <- melt(sqrt(subset(gene_noise_SPS, padj_sigma <= 0.05)[,c(4:6)]))
df2_sps_mu <- aggregate(value ~ variable, df2_sps, mean)$value
df2_sps_upp_sd <- aggregate(value ~ variable, data=df2_sps, function(x) mean(x) + sd(x))$value
df2_sps_lwr_sd <- aggregate(value ~ variable, data=df2_sps, function(x) mean(x) - sd(x))$value

# Mean gene expression noise (variance) +/- std.dev for H1N+ (see CAP_Sepsis_H1N1_gene_noise_Figure1.r)
df1_h1n <- melt(sqrt(gene_noise_H1N[,c(4:6)]))
df1_h1n_mu <- aggregate(value ~ variable, df1_h1n, mean)$value
df1_h1n_upp_sd <- aggregate(value ~ variable, data=df1_h1n, function(x) mean(x) + sd(x))$value
df1_h1n_lwr_sd <- aggregate(value ~ variable, data=df1_h1n, function(x) mean(x) - sd(x))$value

# Mean gene expression noise (variance) +/- std.dev for genes for which variance changes significantly in H1N1
# (see CAP_Sepsis_H1N1_gene_noise_Figure1.r)
df2_h1n <- melt(sqrt(subset(gene_noise_H1N, padj_sigma <= 0.05)[,c(4:6)]))
df2_h1n_mu <- aggregate(value ~ variable, df2_h1n, mean)$value
df2_h1n_upp_sd <- aggregate(value ~ variable, data=df2_h1n, function(x) mean(x) + sd(x))$value
df2_h1n_lwr_sd <- aggregate(value ~ variable, data=df2_h1n, function(x) mean(x) - sd(x))$value

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("Figure1A.pdf", paper="a4r", width=5.5, height=10)

par(mfrow=c(2,3))

stripchart(value ~ variable, data = df1_cap, method="jitter", vertical = TRUE,
  col = NULL, pch=NA, cex=0.3, ylab = "sigma", las=1, ylim=c(0,0.95))
lines(c(1:3), df1_cap_mu, col=1, lty=2, lwd=2  )
segments(c(1:3), df1_cap_upp_sd, c(1:3), df1_cap_lwr_sd, lwd=3, col = c("#21A883", "#4591CC", "#673AB7")  )
points(c(1:3), df1_cap_mu, pch=21, col = c("#21A883", "#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("bottomright", legend =
c(paste0("p(1-0) = ", round(t.test(sqrt(gene_noise_CAP$v_ctl),  sqrt(gene_noise_CAP$v_cap0))$p.value, 5)),
  paste0("p(2-1) = ", round(t.test(sqrt(gene_noise_CAP$v_cap1), sqrt(gene_noise_CAP$v_cap0))$p.value, 5)))
)

stripchart(value ~ variable, data = df1_sps, method="jitter", vertical = TRUE,
  col = NULL, pch=NA, cex=0.3, ylab = "sigma", las=1, ylim=c(0,0.95))
lines(c(1:3), df1_sps_mu, col=1, lty=2, lwd=2  )
segments(c(1:3), df1_sps_upp_sd, c(1:3), df1_sps_lwr_sd, lwd=3, col = c("#21A883", "#4591CC", "#673AB7")  )
points(c(1:3), df1_sps_mu, pch=21, col = c("#21A883", "#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("bottomright", legend =
c(paste0("p(1-0) = ", round(t.test(sqrt(gene_noise_SPS$v_ctl ), sqrt(gene_noise_SPS$v_sps0))$p.value, 5)),
  paste0("p(2-1) = ", round(t.test(sqrt(gene_noise_SPS$v_sps1), sqrt(gene_noise_SPS$v_sps0))$p.value, 5)) )
)

stripchart(value ~ variable, data = df1_h1n, method="jitter", vertical = TRUE,
col = NULL, pch=NA, cex=0.3, ylab = "sigma", las=1, ylim=c(0,0.65))
lines(c(1:3), df1_h1n_mu, col=1, lty=2, lwd=2  )
segments(c(1:3), df1_h1n_upp_sd, c(1:3), df1_h1n_lwr_sd, lwd=3, col = c("#21A883", "#4591CC", "#673AB7")  )
points(c(1:3), df1_h1n_mu, pch=21, col = c("#21A883", "#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("bottomright", legend =
c(paste0("p(1-0) = ", round(t.test(sqrt(gene_noise_H1N$v_ctl ), sqrt(gene_noise_H1N$v_h1n0))$p.value, 5)),
  paste0("p(2-1) = ", round(t.test(sqrt(gene_noise_H1N$v_h1n1), sqrt(gene_noise_H1N$v_h1n0))$p.value, 5)) )
)


stripchart(value ~ variable, data = df2_cap, method="jitter", vertical = TRUE,
col = NULL, pch=NA, cex=0.3, ylab = "sigma, padj <= 0.05", las=1, ylim=c(0,1.1))
lines(c(1:3), df2_cap_mu, col=1, lty=2, lwd=2  )
segments(c(1:3), df2_cap_upp_sd, c(1:3), df2_cap_lwr_sd, lwd=3, col = c("#21A883", "#4591CC", "#673AB7")  )
points(c(1:3), df2_cap_mu, pch=21, col = c("#21A883", "#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("bottomright", legend =
c(paste0("p(1-0) = ", round(t.test(sqrt(subset(gene_noise_CAP, padj_sigma <= 0.05)$v_ctl),  sqrt(subset(gene_noise_CAP, padj_sigma <= 0.05)$v_cap0))$p.value, 5)),
  paste0("p(2-1) = ", round(t.test(sqrt(subset(gene_noise_CAP, padj_sigma <= 0.05)$v_cap1), sqrt(subset(gene_noise_CAP, padj_sigma <= 0.05)$v_cap0))$p.value, 5)))
)

stripchart(value ~ variable, data = df2_sps, method="jitter", vertical = TRUE,
col = NULL, pch=NA, cex=0.3, ylab = "sigma, padj <= 0.05", las=1, ylim=c(0,1.1))
lines(c(1:3), df2_sps_mu, col=1, lty=2, lwd=2  )
segments(c(1:3), df2_sps_upp_sd, c(1:3), df2_sps_lwr_sd, lwd=3, col = c("#21A883", "#4591CC", "#673AB7")  )
points(c(1:3), df2_sps_mu, pch=21, col = c("#21A883", "#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("bottomright", legend =
c(paste0("p(1-0) = ", round(t.test(sqrt(subset(gene_noise_SPS, padj_sigma <= 0.05)$v_ctl),  sqrt(subset(gene_noise_SPS, padj_sigma <= 0.05)$v_sps0))$p.value, 5)),
  paste0("p(2-1) = ", round(t.test(sqrt(subset(gene_noise_SPS, padj_sigma <= 0.05)$v_sps1), sqrt(subset(gene_noise_SPS, padj_sigma <= 0.05)$v_sps0))$p.value, 5)))
)

stripchart(value ~ variable, data = df2_h1n, method="jitter", vertical = TRUE,
col = NULL, pch=NA, cex=0.3, ylab = "sigma, padj <= 0.05", las=1, ylim=c(0,0.9))
lines(c(1:3), df2_h1n_mu, col=1, lty=2, lwd=2  )
segments(c(1:3), df2_h1n_upp_sd, c(1:3), df1_h1n_lwr_sd, lwd=3, col = c("#21A883", "#4591CC", "#673AB7")  )
points(c(1:3), df2_h1n_mu, pch=21, col = c("#21A883", "#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("bottomright", legend =
c(paste0("p(1-0) = ", round(t.test(sqrt(subset(gene_noise_H1N, padj_sigma <= 0.05)$v_ctl),  sqrt(subset(gene_noise_H1N, padj_sigma <= 0.05)$v_h1n0))$p.value, 5)),
  paste0("p(2-1) = ", round(t.test(sqrt(subset(gene_noise_H1N, padj_sigma <= 0.05)$v_h1n1), sqrt(subset(gene_noise_H1N, padj_sigma <= 0.05)$v_h1n0))$p.value, 5)))
)

dev.off()


################################################################################


library(Rfit)
library(quantreg)
library(mvoutlier)

# Fluctuation-responce: variance in gene expression correlates with response in mean gene expression
# See main text for details

pdf("Figure1BCD.pdf", width=11, height=12.75)

par(mfrow=c(3,3))

# CAP
dfd <-  data.frame(v = gene_noise_CAP$V_ave, d = abs(gene_noise_CAP$mu2_cap - gene_noise_CAP$mu2_ctl))
R_dfd <- cor.test(dfd$d, dfd$v, method="spearman")
plot(d ~ v, data=dfd, las=1, xlab="variance", ylab="difference mu CAP - ctl", xlim=c(5e-03,5), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfd)[ sign2(log10(dfd), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dfv <-  data.frame(v = gene_noise_CAP$V_ave, d = abs(gene_noise_CAP$v2_cap - gene_noise_CAP$v2_ctl))
R_dfv <- cor.test(dfv$d, dfv$v, method="spearman")
plot(d ~ v, data=dfv, las=1, xlab="variance", ylab="difference v CAP - ctl", xlim=c(5e-03,5), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfv)[ sign2(log10(dfv), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfv$estimate, 3)), paste0("p = ", round(R_dfv$p.value, 5)) ) )

dfdv <-  data.frame(v = abs(gene_noise_CAP$v2_cap - gene_noise_CAP$v2_ctl), d = abs(gene_noise_CAP$mu2_cap - gene_noise_CAP$mu2_ctl))
R_dfdv <- cor.test(dfdv$d, dfdv$v, method="spearman")
plot(v ~ d, data=dfdv, las=1, xlab="difference mu CAP - ctl", ylab="difference v CAP - ctl",
     xlim=c(1e-03,10), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(v ~ d, data =  log10(dfdv)[ sign2(log10(dfdv), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfdv$estimate, 3)), paste0("p = ", round(R_dfdv$p.value, 5)) ) )

# Sepsis
dfd <-  data.frame(v = gene_noise_SPS$V_ave, d = abs(gene_noise_SPS$mu2_sps - gene_noise_SPS$mu2_ctl))
R_dfd <- cor.test(dfd$d, dfd$v, method="spearman")
plot(d ~ v, data=dfd, las=1, xlab="variance", ylab="difference mu SPS - ctl", xlim=c(5e-03,5), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfd)[ sign2(log10(dfd), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dfv <-  data.frame(v = gene_noise_SPS$V_ave, d = abs(gene_noise_SPS$v2_sps - gene_noise_SPS$v2_ctl))
R_dfv <- cor.test(dfv$d, dfv$v, method="spearman")
plot(d ~ v, data=dfv, las=1, xlab="variance", ylab="difference v SPS - ctl", xlim=c(5e-03,5), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfv)[ sign2(log10(dfv), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfv$estimate, 3)), paste0("p = ", round(R_dfv$p.value, 5)) ) )

dfdv <-  data.frame(v = abs(gene_noise_SPS$v2_sps - gene_noise_SPS$v2_ctl), d = abs(gene_noise_SPS$mu2_sps - gene_noise_SPS$mu2_ctl))
R_dfdv <- cor.test(dfdv$d, dfdv$v, method="spearman")
plot(v ~ d, data=dfdv, las=1, xlab="difference mu SPS - ctl", ylab="difference v SPS - ctl",
     xlim=c(1e-03,10), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(v ~ d, data =  log10(dfdv)[ sign2(log10(dfdv), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfdv$estimate, 3)), paste0("p = ", round(R_dfdv$p.value, 5)) ) )

# H1N1
dfd <-  data.frame(v = gene_noise_H1N$V_ave, d = abs(gene_noise_H1N$mu2_h1n - gene_noise_H1N$mu2_ctl))
R_dfd <- cor.test(dfd$d, dfd$v, method="spearman")
plot(d ~ v, data=dfd, las=1, xlab="variance", ylab="difference mu H1N1 - ctl", xlim=c(5e-03,5), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfd)[ sign2(log10(dfd), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfd$estimate, 3)), paste0("p = ", round(R_dfd$p.value, 5)) ) )

dfv <-  data.frame(v = gene_noise_H1N$V_ave, d = abs(gene_noise_H1N$v2_h1n - gene_noise_H1N$v2_ctl))
R_dfv <- cor.test(dfv$d, dfv$v, method="spearman")
plot(d ~ v, data=dfv, las=1, xlab="variance", ylab="difference v H1N1 - ctl", xlim=c(5e-03,5), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(d ~ v,, data =  log10(dfv)[ sign2(log10(dfv), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfv$estimate, 3)), paste0("p = ", round(R_dfv$p.value, 5)) ) )

dfdv <-  data.frame(v = abs(gene_noise_H1N$v2_h1n - gene_noise_H1N$v2_ctl), d = abs(gene_noise_H1N$mu2_h1n - gene_noise_H1N$mu2_ctl))
R_dfdv <- cor.test(dfdv$d, dfdv$v, method="spearman")
plot(v ~ d, data=dfdv, las=1, xlab="difference mu H1N1 - ctl", ylab="difference v H1N1 - ctl",
     xlim=c(1e-03,10), ylim=c(1e-03,10), pch=20, col="#00000022", log="xy", cex=0.7)
abline(rq(v ~ d, data =  log10(dfdv)[ sign2(log10(dfdv), explvar = 0.5, qcrit = 0.5)$wfinal01 == 0, ]), col=2, lwd =2)
legend("topleft", legend=c( paste0("r = ", round(R_dfdv$estimate, 3)), paste0("p = ", round(R_dfdv$p.value, 5)) ) )

dev.off()
