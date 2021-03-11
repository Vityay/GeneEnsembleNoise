# Ensemble gene noise for the CAP/Sepsis

library(parallel)

setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
# Read CORUM and KEGG gene ensembles annotations
load("org.Hs.core_paths.Rdat")
load("pr_cl.Rdat")
names(kegg_path) <- paste0("kegg:", names(kegg_path))
# redundant to Respiratory chain complex I (subcomplex I alpha), mitochondrial complex
corum_core$`mitochondrial respiratory chain complex I assembly` <- NULL
names(corum_core) <- paste0("complex:", names(corum_core))
core_paths <- c(corum_core, kegg_path, promSim_cl)
# promSim_cl - gene ensembles based on promoter similarity


setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
# Read CAP/Sepsis affymetrix data
load("GSE65682_series_all_mortality.Rdat")

# match genes between ensembles annotations and affymetrix data
core_paths_subset <- lapply(seq_along(core_paths), function(i) core_paths[[i]][ core_paths[[i]] %in% rownames(affyData_out) ] )
names(core_paths_subset) <- names(core_paths)
# select ensembles represented by more than 10 genes
core_paths_subset <- core_paths_subset[ sapply(core_paths_subset, length) >= 10 ]

# Calculate ensemble gene noise
affyDat_noise_paths <- lapply(seq_along(core_paths_subset), function(i) {
  cat(i, "\r")
  idx <- rownames(affyData_out) %in% core_paths_subset[[i]]
  dat <- affyData_out[idx,]
  res <- apply(dat, 2, var)
  res
})
names(affyDat_noise_paths) <- names(core_paths_subset)
affyDat_noise_paths <- do.call(cbind, affyDat_noise_paths)

affyDat_noise_paths <- affyDat_noise_paths[pheno$d != "ctl", ]
pheno <- pheno[pheno$d != "ctl", ]

# save ensemble gene noise for CAP
affyDat_noise_paths_CAP <- affyDat_noise_paths[pheno$d != "other", ]
pheno_CAP <- pheno[pheno$d != "other", ]

# save ensemble gene noise for Sepsis
affyDat_noise_paths_SPS <- affyDat_noise_paths[pheno$d != "cap", ]
pheno_SPS <- pheno[pheno$d != "cap", ]

rm(affyData, affyData_out, core_paths, corum_core, kegg_path, org.Hs, promSim_cl, samples, samples_pData, samples_pData_subset, trxL)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("CAP_Sepsis_Ensemble_noise_model.RDat")

