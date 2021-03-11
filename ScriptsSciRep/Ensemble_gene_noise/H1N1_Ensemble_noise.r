# Ensemble gene noise for the H1N1

library(parallel)

setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
# Read CORUM and KEGG gene ensembles annotations
load("org.Hs.core_paths.Rdat")
names(kegg_path) <- paste0("kegg:", names(kegg_path))
# redundant to Respiratory chain complex I (subcomplex I alpha), mitochondrial complex
corum_core$`mitochondrial respiratory chain complex I assembly` <- NULL
names(corum_core) <- paste0("complex:", names(corum_core))
core_paths <- c(corum_core, kegg_path)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
# Read H1N1 Illumina data
load("eset_H1N1.Rdat")

# match genes between ensembles annotations and affymetrix data
core_paths_H1N1 <- lapply(seq_along(core_paths), function(i) core_paths[[i]][ core_paths[[i]] %in% rownames(eset_H1N1) ] )
names(core_paths_H1N1) <- names(core_paths)
# select ensembles represented by more than 10 genes
core_paths_H1N1 <- core_paths_H1N1[ sapply(core_paths_H1N1, length) >= 10 ]

# Calculate ensemble gene noise
H1N1_noise_paths <- lapply(seq_along(core_paths_H1N1), function(i) {
  cat(i, "\r")
  idx <- rownames(eset_H1N1) %in% core_paths_H1N1[[i]]
  dat <- eset_H1N1[idx,]
  res <- apply(dat, 2, var)
  res
})
names(H1N1_noise_paths) <- names(core_paths_H1N1)
H1N1_noise_paths <- do.call(cbind, H1N1_noise_paths)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("H1N1_Ensemble_noise.RDat")

################################################################################
# Estimate significance of association of ensemble gene noise with H1N1

library(SuppDists)
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("H1N1_Ensemble_noise.RDat")

H1N1_stage <- samples_H1N1$'disease phase:ch1'
pR_H1N1_stage <- apply(H1N1_noise_paths, 2, function(x) {

  q = cor(x, as.numeric(factor(H1N1_stage)), method="kendall")
  r = length(H1N1_stage)
  res <- c(tau = q, p = pKendall(q,r, lower=F))

})
pR_H1N1_stage <- data.frame(t(pR_H1N1_stage))
padjR_H1N1_stage <- pR_H1N1_stage
padjR_H1N1_stage$p <- p.adjust(pR_H1N1_stage$p, "BH")
# padjR_H1N1_stage[grep("mitoch", rownames(padjR_H1N1_stage)),]

padjR_H1N1_stage_sig <- padjR_H1N1_stage[padjR_H1N1_stage$p <= 0.05,]
padjR_H1N1_stage_sig_UP <- padjR_H1N1_stage_sig[padjR_H1N1_stage_sig$tau > 0,]
# There are no significantly down-regulated pathways
# padjR_H1N1_stage_sig_DOWN <- padjR_H1N1_stage_sig[padjR_H1N1_stage_sig$tau < 0,]

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("H1N1_Ensemble_noise_sig.RDat")
