library(edgeR)
library(gamlss)
library(parallel)

################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("COVID19_Ensemble_noise.RDat")

# Adjust Ensemble gene noise p values
COVID_ENS_noise_sig$padj <- p.adjust(COVID_ENS_noise_sig$p, "bonferroni")

################################################################################
# Select pathways to plot
pathways_to_plot <- c(
  grep("mitochondrial respiratory chain complex I assembly", rownames(COVID_ENS_noise_sig), fixed = T),
  # grep("Respiratory chain complex I (subcomplex I alpha)", rownames(COVID_ENS_noise_sig), fixed = T),
  grep("Respiratory chain complex I (nuclear encoded subunits)", rownames(COVID_ENS_noise_sig), fixed = T),
  grep("Osteoclast differentiation", rownames(COVID_ENS_noise_sig)),
  grep("Tight junction", rownames(COVID_ENS_noise_sig)),
  grep("Axon guidance", rownames(COVID_ENS_noise_sig)),
  grep("HIF-1 signaling pathway", rownames(COVID_ENS_noise_sig)),
  grep("Peroxisome", rownames(COVID_ENS_noise_sig)),
  grep("Necroptosis", rownames(COVID_ENS_noise_sig)),
  grep("NOD-like receptor signaling pathway", rownames(COVID_ENS_noise_sig)),
  grep("Fc epsilon RI signaling pathway", rownames(COVID_ENS_noise_sig)),
  grep("Autophagy - other", rownames(COVID_ENS_noise_sig)),
  grep("Biosynthesis of amino acids", rownames(COVID_ENS_noise_sig)),
  grep("TRAPP complex", rownames(COVID_ENS_noise_sig)),
  grep("Glucagon signaling pathway", rownames(COVID_ENS_noise_sig)),
  grep("Propanoate metabolism", rownames(COVID_ENS_noise_sig)),
  grep("Circadian rhythm", rownames(COVID_ENS_noise_sig)),
  grep("Dopaminergic synapse", rownames(COVID_ENS_noise_sig)),
  grep("Amyotrophic lateral sclerosis", rownames(COVID_ENS_noise_sig))
)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("Figure6B_S10.pdf", paper="a4r", width=8, height=5)

stripcol = c("#4591CCCC", "#673AB7CC"); COL = c("#4591CC", "#673AB7")

par(mfrow=c(1,4))

for(i in seq(pathways_to_plot)) {

  pthw <- rownames(COVID_ENS_noise)[pathways_to_plot[i]]
  st <- COVID_ENS_noise_sig[pathways_to_plot[i],]
  df <- data.frame(y = c(COVID_ENS_noise[pathways_to_plot[i],]), x = DGE_RNAcounts$samples$stage)

  mu <- aggregate(y ~ x, data=df, function(x) mean(x))$y
  upp_sd <- aggregate(y ~ x, data=df, function(x) mean(x) + sd(x))$y
  lwr_sd <- aggregate(y ~ x, data=df, function(x) mean(x) - sd(x))$y

  stripchart(df$y ~ df$x, vertical=T, method="jitter", pch=c(16, 17),
    col=stripcol, add=F, cex=2, ylab = pthw, las = 1)

  segments(1, mu[1], 2, mu[2], col = 1, lty = 2, lwd=2  )
  segments(c(1:2), upp_sd, c(1:2), lwr_sd, col = COL, lwd=3  )
  points(c(1:2), mu, pch=21, col = COL, bg = "white", cex=3, lwd=3  )
  legend("topleft", legend = paste0("t = ", round(st$t, 3), ", p = ", round(st$padj, 5)), cex=0.9)
}

dev.off()
