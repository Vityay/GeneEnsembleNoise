library(viridis)
library(heatmap3)
library(gplots)
library(superheat)
library(WGCNA)


setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("CAP_sepsis_WGCNA.Rdat")

hmap <- cor(datME)
idx <- stats::hclust(dist(hmap), method = "ward.D2")$order
idx <- rev(idx)
hmap <- hmap[idx, idx]
hmap <- hmap[1:nrow(hmap), nrow(hmap):1]
hmap_col <- col2hex(gsub("ME", "", rownames(hmap)))
rownames(hmap) <- NULL
colnames(hmap) <- NULL

datME_sig <- datME
datME_sig <- datME_sig[,idx]
datME_sig <- apply(datME_sig, 2, function(x) t.test(x ~ phenotype)$p.value)
datME_sig <- p.adjust(datME_sig, "BH")
datME_sig <- -log10(datME_sig)
# datME_sig[grep("skyblue", names(datME_sig))]

# datME_sig[datME_sig < -log10(0.05)] <- 0
# datME_sig <- apply(datME_sig, 2, function(x) dcor.test(x, phenotype, R = 1000)$p.value)

datME_sig_col <- col2hex(gsub("ME", "", names(datME_sig)))
names(datME_sig) <- NULL

p1 <- superheat(hmap,
  # set heatmap color map
    heat.pal = viridis(100),
    heat.pal.values = seq(0.025, 1, length.out = 100),
    # heat.col.scheme = "viridis",
    heat.lim = c(0.5, 1),
    heat.na.col = "black",
    grid.hline = FALSE,
    grid.vline = FALSE,
    left.label.text.size = -1,
    bottom.label.text.size = -1,
    left.label.col = hmap_col,
    bottom.label.col = rev(hmap_col),
  # right plot: HDI
    yr = datME_sig,
    yr.plot.type = "bar",
    yr.axis.name = "-log10(FDR)",
    yr.bar.col = 1,
    yr.obs.col = datME_sig_col,
    yr.lim = c(0, -log10(1e-03))
)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures/")
ggplot2::ggsave("FigureS5A.pdf", p1$plot, height = 12, width = 12)
