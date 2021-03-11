library(venn)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("CAP_Sepsis_Ensemble_noise_sig.RDat")
load("H1N1_Ensemble_noise_sig.RDat")
load("COVID19_Ensemble_noise.RDat")

################################################################################

# Respiratory chain complex I (subcomplex I alpha) is redundant, although gene sets differ slightly
# We choose mitochondrial respiratory chain complex I assembly for the ensemble gene noise for the COVID19 patients
# as it has bigger association with COVID19 mild/sever stage
idx1 <- grep("Respiratory chain complex I (subcomplex I alpha)", rownames(COVID_ENS_noise_sig), fixed = T )
idx2 <- grep("mitochondrial respiratory chain complex I assembly", rownames(COVID_ENS_noise_sig), fixed = T )
COVID_ENS_noise_sig[idx1,] <- COVID_ENS_noise_sig[idx2,]
COVID_ENS_noise_sig <- COVID_ENS_noise_sig[-idx2,]
COVID_ENS_noise_sig$padj <- p.adjust(COVID_ENS_noise_sig$p, "bonferroni")


################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")

pdf("Figure2A.pdf")

venn(list( 'COVID-19' = rownames(COVID_ENS_noise_sig)[COVID_ENS_noise_sig$padj < 0.01],
  H1N1 = rownames(padjR_H1N1_stage_sig_UP),
                CAP = rownames(padjR_CAP_stage_sig_UP),
                sepsis = rownames(padjR_SPS_stage_sig_UP)), ilabels = TRUE, zcolor = "style"
)

dev.off()
