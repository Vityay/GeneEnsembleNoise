setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
load("CAP_Sepsis_Ensemble_noise.RDat")

var_affyDat_noise_paths <- apply(affyDat_noise_paths, 2, function(x) var(as.numeric(x)))
var_affyData_out <- apply(affyData_out, 1, function(x) var(as.numeric(x)))
z <- as.numeric(affyData_out[grep("GAPDH", rownames(affyData_out)), ])
var_affyData_out_GAPDH <- apply(affyData_out, 1, function(x) var(as.numeric(x)-z))

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures")
pdf("FigureS3.pdf", width = 3, height = 5)

boxplot(list(v_genes = var_affyData_out, v_noise = var_affyDat_noise_paths, 'v_genes:GAPDH' = var_affyData_out_GAPDH), outline=F, las=1)
legend("topleft", legend = c("p < 0.0001", "p < 0.0001"))

dev.off()

# t.test(var_affyData_out, var_affyDat_noise_paths)
# t.test(var_affyData_out, var_affyData_out_GAPDH)
