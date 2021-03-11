library(reshape2)


setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("COVID19_gene_noise_Figure6A.Rdat")

# Remove outlier genes from analysis with extremely high BCV
gene_noise_COVID$bcv_mild[gene_noise_COVID$bcv_mild >= 10] <- NA
gene_noise_COVID$bcv_severe[gene_noise_COVID$bcv_severe >= 10] <- NA

# data frame BCV mild/severe
dat <- gene_noise_COVID[,c(3:4)]
colnames(dat) <- c("mild", "severe")
dat <- na.omit(dat)

# Mean gene expression BCV +/- std.dev for COVID19
df1_COVID <- melt(dat)
df1_COVID_mu <- aggregate(value ~ variable, df1_COVID, mean)$value
df1_COVID_upp_sd <- aggregate(value ~ variable, data=df1_COVID, function(x) mean(x) + sd(x))$value
df1_COVID_lwr_sd <- aggregate(value ~ variable, data=df1_COVID, function(x) mean(x) - sd(x))$value

# Mean gene expression BCV +/- std.dev for COVID19 for which BCV changes significantly between mild/severe COVID19
dat2 <- gene_noise_COVID[p.adjust(gene_noise_COVID$p_bcv, "BH") < 0.05,c(3:4)]
colnames(dat2) <- c("mild", "severe")
dat2 <- na.omit(dat2)

df2_COVID <- melt(dat2)
df2_COVID_mu <- aggregate(value ~ variable, df2_COVID, mean)$value
df2_COVID_upp_sd <- aggregate(value ~ variable, data=df2_COVID, function(x) mean(x) + sd(x))$value
df2_COVID_lwr_sd <- aggregate(value ~ variable, data=df2_COVID, function(x) mean(x) - sd(x))$value

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("Figure6A.pdf", paper="a4r", width=5.5, height=8)

par(mfrow=c(1,2))

stripchart(value ~ variable, data = df1_COVID, method="jitter", vertical = TRUE,
  col = NULL, pch=NA, cex=0.3, ylab = "bcv", las=1, ylim=c(0.1,0.7))
lines(c(1:2), df1_COVID_mu, col=1, lty=2, lwd=2  )
segments(c(1:2), df1_COVID_upp_sd, c(1:2), df1_COVID_lwr_sd, lwd=3, col = c("#4591CC", "#673AB7")  )
points(c(1:2), df1_COVID_mu, pch=21, col = c("#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("topleft", legend = c(
  paste0("p(1-0) = ", round(t.test(sqrt(gene_noise_COVID$bcv_mild),  sqrt(gene_noise_COVID$bcv_severe))$p.value, 5))
))


stripchart(value ~ variable, data = df2_COVID, method="jitter", vertical = TRUE,
  col = NULL, pch=NA, cex=0.3, ylab = "bcv", las=1, ylim=c(0.1,0.7))
lines(c(1:2), df2_COVID_mu, col=1, lty=2, lwd=2  )
segments(c(1:2), df2_COVID_upp_sd, c(1:2), df2_COVID_lwr_sd, lwd=3, col = c("#4591CC", "#673AB7")  )
points(c(1:2), df2_COVID_mu, pch=21, col = c("#4591CC", "#673AB7"), bg = "white", cex=3, lwd=3  )
legend("topleft", legend = c(
  paste0("p(1-0) = ", round(t.test(sqrt(gene_noise_COVID$bcv_mild),  sqrt(gene_noise_COVID$bcv_severe))$p.value, 5))
))

dev.off()
