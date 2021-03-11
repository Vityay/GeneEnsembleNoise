# Figure S11B - WGCNA pathway enrichment

setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
load("org.Hs.core_paths.Rdat")

names(kegg_path) <- paste0("kegg:", names(kegg_path))
names(corum_core) <- paste0("complex:", names(corum_core))
corum_core$`complex:mitochondrial respiratory chain complex I assembly` <- NULL
core_paths <- c(corum_core, kegg_path)

################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
load("COVID_WGCNA.Rdat")

hmap <- cor(datME)
idx <- stats::hclust(dist(hmap), method = "ward.D2")$order
idx <- rev(idx)

datME_sig <- datME
datME_sig <- datME_sig[,idx]
datME_sig <- apply(datME_sig, 2, function(x) t.test(x ~ phenotype)$p.value)
datME_sig <- p.adjust(datME_sig, "fdr")
datME_sig <- -log10(datME_sig)

datME_sig <- datME_sig[datME_sig >= -log10(0.05)]

sigME <- names(datME_sig <- datME_sig[datME_sig >= -log10(0.05)])
sigME <- gsub("ME", "", sigME)

################################################################################
# For explanetion of gene set enrichment see
# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
N <- length(unlist(core_paths)[!duplicated(unlist(core_paths))])
m <- sapply(core_paths, length)
n <- N-m

sigME_ens <- lapply(sigME, function(module) {

  k <- sum(GeneNames[dynamicColors %in% module] %in% unlist(core_paths)[!duplicated(unlist(core_paths))])
  x <- sapply(core_paths, function(x) sum(x %in% GeneNames[dynamicColors %in% module]) )

  sig <- phyper(q = x -1, m=m, n=n, k=k, lower.tail=FALSE)
  sig <- p.adjust(sig, "fdr")
  res <- sig[sig < 0.001]

  if(length(res) > 0) {
    res <- data.frame(module = module, ensemble = names(res), fdr = -log10(res))
    res <- res[order(res$fdr, decreasing = T),]
    rownames(res) <- NULL
  } else {res <- NULL}

  res

})

sigME_ens <- do.call(rbind, sigME_ens)
sigME_ens <- sigME_ens[!duplicated(sigME_ens$ensemble),]

sigME_ens <- sigME_ens[nrow(sigME_ens):1,]
sigME_ens$fdr[sigME_ens$fdr > -log10(1e-10)] <- -log10(1e-10)
sigME_ens <- sigME_ens[nrow(sigME_ens):1,]

################################################################################

library(gplots)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")
pdf("FigureS11B.pdf", width=8, height=6.4)

par(mar=c(4,16,1,0.5))
barplot(sigME_ens[,'fdr'], horiz = T, names.arg = sigME_ens[,'ensemble'], col = col2hex(sigME_ens[,'module']),
      xlab = "-log10(fdr)", main = "WGCNA enrichment (fdr < 0.01)",las=1, cex.axis=0.7, cex.names = 0.7)

dev.off()
