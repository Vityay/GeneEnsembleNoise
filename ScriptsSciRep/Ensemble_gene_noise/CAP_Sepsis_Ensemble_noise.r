# Ensemble gene noise for the CAP/Sepsis

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

# save ensemble gene noise for CAP
affyDat_noise_paths_CAP <- affyDat_noise_paths[pheno$d != "other", ]
pheno_CAP <- pheno[pheno$d != "other", ]

# save ensemble gene noise for Sepsis
affyDat_noise_paths_SPS <- affyDat_noise_paths[pheno$d != "cap", ]
pheno_SPS <- pheno[pheno$d != "cap", ]

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("CAP_Sepsis_Ensemble_noise.RDat")

################################################################################
# Estimate significance of association of ensemble gene noise with CAP and Sepsis

library(SuppDists)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("CAP_Sepsis_Ensemble_noise.RDat")

# Estimate significance of association of ensemble gene noise with CAP
pheno_CAP$mort_num[pheno_CAP$mort == "NA"] <- -1

pR_CAP_stage <- apply(affyDat_noise_paths_CAP, 2, function(x) {
  df <- data.frame(x=x, age=pheno_CAP$age, d = pheno_CAP$d, mort = pheno_CAP$mort_num)

  # subtract the effect of age
  m00 <- lm(x ~ age, data=df)
  df$xr <- coefficients(m00)[[1]] + resid(m00)

  q = cor(df$xr, df$mort , method="kendall")
  r = nrow(df)
  res <- c(tau = q, p = pKendall(q,r, lower=F))
})
pR_CAP_stage <- data.frame(t(pR_CAP_stage))

padjR_CAP_stage <- pR_CAP_stage
padjR_CAP_stage$p <- p.adjust(padjR_CAP_stage$p, "bonferroni")
# padjR_CAP_stage[grep("mitoch", rownames(padjR_CAP_stage)),]

padjR_CAP_stage_sig <- padjR_CAP_stage[padjR_CAP_stage$p <= 0.05,]
padjR_CAP_stage_sig_UP <- padjR_CAP_stage_sig[padjR_CAP_stage_sig$tau > 0,]
padjR_CAP_stage_sig_DOWN <- padjR_CAP_stage_sig[padjR_CAP_stage_sig$tau < 0,]

# Estimate significance of association of ensemble gene noise with Sepsis
pheno_SPS$mort_num[pheno_SPS$mort == "NA"] <- -1

pR_SPS_stage <- apply(affyDat_noise_paths_SPS, 2, function(x) {
  df <- data.frame(x=x, age=pheno_SPS$age, d = pheno_SPS$d, mort = pheno_SPS$mort_num)

  # subtract the effect of age
  m00 <- lm(x ~ age, data=df)
  df$xr <- coefficients(m00)[[1]] + resid(m00)

  q = cor(df$xr, df$mort , method="kendall")
  r = nrow(df)
  res <- c(tau = q, p = pKendall(q,r, lower=F))
})
pR_SPS_stage <- data.frame(t(pR_SPS_stage))
padjR_SPS_stage <- pR_SPS_stage
padjR_SPS_stage$p <- p.adjust(padjR_SPS_stage$p, "bonferroni")
# padjR_SPS_stage[grep("mitoch", rownames(padjR_SPS_stage)),]

padjR_SPS_stage_sig <- padjR_SPS_stage[padjR_SPS_stage$p <= 0.05,]
padjR_SPS_stage_sig_UP <- padjR_SPS_stage_sig[padjR_SPS_stage_sig$tau > 0,]
padjR_SPS_stage_sig_DOWN <- padjR_SPS_stage_sig[padjR_SPS_stage_sig$tau < 0,]

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("CAP_Sepsis_Ensemble_noise_sig.RDat")
