library(quantreg)
library(Rfit)

# plot function
plot_df <- function(df, rankscale = FALSE, ylim = NULL, stripcol = c("#21A883CC", "#4591CCCC", "#673AB7CC"), COL = c("#21A883", "#4591CC", "#673AB7") ) {
  path <- rownames(df$tau)
  tau <- df$tau
  df <- df$df; df$x <- df$x + 1

  if(rankscale == TRUE) {
    range_y <- range(df$y)
    rank_y <- rank(df$y)
    rank_y_sc <- (rank_y - min(rank_y))/(max(rank_y)-min(rank_y) )
    df$y <- (rank_y_sc)*(max(range_y)-min(range_y)) + min(range_y)
  }

  q0 = round(tau$tau, 3)
  ptau = round(tau$p, 5)

  X0 <- seq(1, 3, 0.1)
  mm <- model.matrix( ~ X0)
  m1 <- rfit(y ~ x, df)
  pred.m1 <- coef(m1)[[1]] + X0*coef(m1)[[2]]
  pred.se <- apply(mm, 1, function(Xrow) t(Xrow) %*% vcov(m1) %*% Xrow)
  upp <- pred.m1 + 1.96*sqrt(pred.se)
  lwr <- pred.m1 - 1.96*sqrt(pred.se)

  upp_sd <- aggregate(y ~ x, data=df, function(x) mean(x) + sd(x))$y
  lwr_sd <- aggregate(y ~ x, data=df, function(x) mean(x) - sd(x))$y

  if(is.null(ylim)) { ylim <- range(c(df$y, upp, lwr, upp_sd, lwr_sd))  }

  stripchart(df$y ~ I(df$x-1), vertical=T, method="jitter", pch=c(15,16, 17), col=stripcol, add=F, cex=2, ylab = colnames(H1N1_noise_paths)[H1N1_idx], ylim = ylim, las=1)
  polygon(c(X0, rev(X0)), c(lwr, rev(upp)), col="#00000033", border = NA);
  lines(X0, pred.m1, col=1, lty=2, lwd=3);
  segments(c(1:3), upp_sd, c(1:3), lwr_sd, col = COL, lwd=3  )
  points(c(1:3), aggregate(y ~ x, data=df, function(x) mean(x))$y, pch=21, col = COL, bg = "white", cex=3, lwd=3  )
  legend("topleft", legend = paste0("H1N1: tau = ", q0, ", p = ", ptau))

}

################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("CAP_Sepsis_Ensemble_noise_sig.RDat")
load("H1N1_Ensemble_noise_sig.RDat")

pheno_CAP$mort_num[pheno_CAP$mort == "NA"] <- -1
pheno_CAP$mort_num <- pheno_CAP$mort_num + 1
pheno_SPS$mort_num[pheno_SPS$mort == "NA"] <- -1
pheno_SPS$mort_num <- pheno_SPS$mort_num + 1


ALL_UP <- Reduce(intersect, list( H1N1 = rownames(padjR_H1N1_stage_sig_UP),
                 CAP = rownames(padjR_CAP_stage_sig_UP),
                 sepsis = rownames(padjR_SPS_stage_sig_UP)))
CAP_H1N1_UP <- Reduce(intersect, list( H1N1 = rownames(padjR_H1N1_stage_sig_UP), CAP = rownames(padjR_CAP_stage_sig_UP) ) )
CAP_H1N1_UP <-  CAP_H1N1_UP[!(CAP_H1N1_UP %in% ALL_UP)]

CAP_SPS_UP <- Reduce(intersect, list( SPS = rownames(padjR_SPS_stage_sig_UP), CAP = rownames(padjR_CAP_stage_sig_UP) ) )
CAP_SPS_UP <-  CAP_SPS_UP[!(CAP_SPS_UP %in% ALL_UP)]

################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("Figure2B.pdf", paper="a4r", width=8, height=5)

for(k in seq(ALL_UP)) {

idx <- ALL_UP[k]
H1N1_idx <- which(colnames(H1N1_noise_paths) == idx)
df0 <- data.frame(y=H1N1_noise_paths[,H1N1_idx] , x = as.numeric(factor(H1N1_stage))-1)

cap_idx <- which(colnames(affyDat_noise_paths_CAP) == idx)
df1 <- data.frame(y = affyDat_noise_paths_CAP[,cap_idx] , x = pheno_CAP$mort_num, age=pheno_CAP$age)
m00 <- lm(y ~ age, data=df1)
df1$y <-coefficients(m00)[[1]] + resid(m00)
df1$age <- NULL

sps_idx <- which(colnames(affyDat_noise_paths_SPS) == idx)
df2 <- data.frame(y = affyDat_noise_paths_SPS[,cap_idx] , x = pheno_SPS$mort_num, age=pheno_SPS$age)
m00 <- lm(y ~ age, data=df2)
df2$y <- coefficients(m00)[[1]] + resid(m00)
df2$age <- NULL

ylim = range(c(df1$y, df2$y))

par(mfrow = c(1,3))
plot_df(list(df = df0, tau = padjR_H1N1_stage[H1N1_idx,]), ylim = NULL, stripcol = c("#21A883CC", "#4591CCCC", "#673AB7CC"), COL=1)
plot_df(list(df = df1, tau = padjR_CAP_stage[cap_idx,]),   ylim = ylim, stripcol = c("#21A88388", "#4591CC88", "#673AB788"), COL=1)
plot_df(list(df = df2, tau = padjR_SPS_stage[sps_idx,]),   ylim = ylim, stripcol = c("#21A88388", "#4591CC88", "#673AB788"), COL=1)
}

dev.off()

################################################################################

pdf("Figure2C.pdf", paper="a4r", width=8, height=5)

for(k in seq(CAP_H1N1_UP)) {

idx <- CAP_H1N1_UP[k]
H1N1_idx <- which(colnames(H1N1_noise_paths) == idx)
df0 <- data.frame(y=H1N1_noise_paths[,H1N1_idx] , x = as.numeric(factor(H1N1_stage))-1)

cap_idx <- which(colnames(affyDat_noise_paths_CAP) == idx)
df1 <- data.frame(y = affyDat_noise_paths_CAP[,cap_idx] , x = pheno_CAP$mort_num, age=pheno_CAP$age)
m00 <- lm(y ~ age, data=df1)
df1$y <-coefficients(m00)[[1]] + resid(m00)
df1$age <- NULL

sps_idx <- which(colnames(affyDat_noise_paths_SPS) == idx)
df2 <- data.frame(y = affyDat_noise_paths_SPS[,cap_idx] , x = pheno_SPS$mort_num, age=pheno_SPS$age)
m00 <- lm(y ~ age, data=df2)
df2$y <- coefficients(m00)[[1]] + resid(m00)
df2$age <- NULL

ylim = range(c(df1$y, df2$y))

par(mfrow = c(1,3))
plot_df(list(df = df0, tau = padjR_H1N1_stage[H1N1_idx,]), ylim = NULL, stripcol = c("#21A883CC", "#4591CCCC", "#673AB7CC"), COL=1)
plot_df(list(df = df1, tau = padjR_CAP_stage[cap_idx,]),   ylim = ylim, stripcol = c("#21A88388", "#4591CC88", "#673AB788"), COL=1)
plot_df(list(df = df2, tau = padjR_SPS_stage[sps_idx,]),   ylim = ylim, stripcol = c("#21A88388", "#4591CC88", "#673AB788"), COL=1)
}

dev.off()


################################################################################
# Supplementary tables
# ALL_UP_table <- lapply(seq(ALL_UP), function(k) {
#   idx <- ALL_UP[k]
#   H1N1_idx <- which(colnames(H1N1_noise_paths) == idx)
#   cap_idx <- which(colnames(affyDat_noise_paths_CAP) == idx)
#   sps_idx <- which(colnames(affyDat_noise_paths_SPS) == idx)
#
#   path <- names(core_paths_subset)[grep(ALL_UP[k], names(core_paths_subset), fixed=T)]
#   genes <- sort(core_paths_subset[grep(ALL_UP[k], names(core_paths_subset), fixed=T)][[1]])
#   tau_h1n <- round(as.numeric(padjR_H1N1_stage[H1N1_idx,1]),3)
#   tau_cap <-   round(as.numeric(padjR_CAP_stage[cap_idx,1]),3)
#   tau_sps <-   round(as.numeric(padjR_SPS_stage[sps_idx,1]),3)
#
#   p_h1n <- round(as.numeric(padjR_H1N1_stage[H1N1_idx,2]),6)
#   p_cap <-   round(as.numeric(padjR_CAP_stage[cap_idx,2]),6)
#   p_sps <-   round(as.numeric(padjR_SPS_stage[sps_idx,2]),6)
#
#   res <- data.frame(path = path,
#     H1N1 = paste(tau_h1n, p_h1n, sep="; "),
#     CAP  = paste(tau_cap, p_cap, sep="; "),
#     SPS  = paste(tau_sps, p_sps, sep="; "),
#     genes = paste(genes, collapse = "; ") )
#   res
# })
# ALL_UP_table <- do.call(rbind, ALL_UP_table)
#
#
# CAP_H1N1_UP_table <- lapply(seq(CAP_H1N1_UP), function(k) {
#   idx <- CAP_H1N1_UP[k]
#   H1N1_idx <- which(colnames(H1N1_noise_paths) == idx)
#   cap_idx <- which(colnames(affyDat_noise_paths_CAP) == idx)
#   sps_idx <- which(colnames(affyDat_noise_paths_SPS) == idx)
#
#   path <- names(core_paths_subset)[grep(CAP_H1N1_UP[k], names(core_paths_subset), fixed=T)]
#   genes <- sort(core_paths_subset[grep(CAP_H1N1_UP[k], names(core_paths_subset), fixed=T)][[1]])
#   tau_h1n <- round(as.numeric(padjR_H1N1_stage[H1N1_idx,1]),3)
#   tau_cap <-   round(as.numeric(padjR_CAP_stage[cap_idx,1]),3)
#   tau_sps <-   round(as.numeric(padjR_SPS_stage[sps_idx,1]),3)
#
#   p_h1n <- round(as.numeric(padjR_H1N1_stage[H1N1_idx,2]),6)
#   p_cap <-   round(as.numeric(padjR_CAP_stage[cap_idx,2]),6)
#   p_sps <-   round(as.numeric(padjR_SPS_stage[sps_idx,2]),6)
#
#   res <- data.frame(path = path,
#     H1N1 = paste(tau_h1n, p_h1n, sep="; "),
#     CAP  = paste(tau_cap, p_cap, sep="; "),
#     SPS  = paste(tau_sps, p_sps, sep="; "),
#     genes = paste(genes, collapse = "; ") )
#   res
# })
# CAP_H1N1_UP_table <- do.call(rbind, CAP_H1N1_UP_table)
#
# write.table(ALL_UP_table, "ALL_UP_table.txt", row.names=F, col.names=T, quote=F, sep="\t")
# write.table(CAP_H1N1_UP_table, "CAP_H1N1_UP_table.txt", row.names=F, col.names=T, quote=F, sep="\t")
