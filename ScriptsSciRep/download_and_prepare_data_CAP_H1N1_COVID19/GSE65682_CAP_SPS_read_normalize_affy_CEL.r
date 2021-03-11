# Script to process affymetrix gene expression dataset GSE65682
# "Genome-wide blood transcriptional profiling in critically ill patients - MARS consortium"
# CAP and sepsis
# PART 1

library(GEOquery)
library(oligo)
library(affycoretools)
library(hgu219.db)
library(pd.hg.u219)
library(mvoutlier)

setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")

# read 	Series Matrix File
# setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
# samples <- getGEO(file = "GSE65682_series_matrix.txt.gz")
# OR, from the url https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/GSE65682_series_matrix.txt.gz
samples <- getGEO("GSE65682")[[1]]

samples_pData <-  pData(phenoData(samples))

# read CEL files
# setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/GSE65682_CEL")
# files <- list.files()[grep("CEL.gz", list.files())]
# OR, download and read CEL files from the url
options(timeout = max(60*60, getOption("timeout")))
getGEOSuppFiles("GSE65682")
setwd("GSE65682/")
untar("GSE65682_RAW.tar")
files <- list.files()[grep("CEL.gz", list.files())]

GSM <- substring(files, 1,10)
files <- files[match(rownames(samples_pData), GSM)]
affyData <- read.celfiles(filenames = files)

# process (RMA-normalization) CEL files
affyData <- oligo::rma(affyData, background = T, normalize=T)

# convert probe IDs to gene symbols
affyData <- annotateEset(affyData, hgu219.db)
genes <- featureData(affyData)@data[,c(1,3)]

# gene expression matrix
affyData <- exprs(affyData)
colnames(affyData) <- substring(colnames(affyData), 1,10)
affyData <- merge(na.omit(genes), affyData, by.x=1, by.y='row.names')
affyData <- affyData[duplicated(affyData$SYMBOL) == FALSE,]
rownames(affyData) <- affyData$SYMBOL
affyData <- affyData[,-c(1:2)]
# sum(samples_pData$geo_accession != colnames(affyData))

rm(genes, files, GSM)

# A problem of arrays is that it is difficult to discriminate expressed genes from non-expressed genes.
# Genes with hybridization intensities close to backgroung are expected to be noisy accross all samples.
# Thus, we used mvoutlier sign2 algorithm to detect such genes and remove them.
out <- sign2(affyData, explvar = 0.95, qcrit = 0.9)
affyData_out <- affyData[out$wfinal01 == 0, ]
# Compare
# plot(density(unlist(affyData[out$wfinal01 == 1, ])), xlim = c(1,14))
# lines(density(unlist(affyData[out$wfinal01 == 0, ])), col=2)

rm(out)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("GSE65682_series_all.Rdat")

################################################################################
################################################################################
################################################################################

# PART 2
# We are interested only in samples with mortality information, thus we repeat
# the above analysis for the selected samples

library(GEOquery)
library(oligo)
library(affycoretools)
library(hgu219.db)
library(pd.hg.u219)
library(mvoutlier)


# read 	Series Matrix File
setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
samples <- getGEO(file = "GSE65682_series_matrix.txt.gz")
samples_pData <-  pData(phenoData(samples))

# Choose samples for which mortality_event_28days is available for CAP and Sepsis
samples_pData_subset <-
subset( samples_pData,
samples_pData$'icu_acquired_infection:ch1' == "healthy" |
samples_pData$'mortality_event_28days:ch1' != "NA"
)

# Simlify phenotypic information
pheno <- data.frame( sex = samples_pData_subset$'gender:ch1',
age = samples_pData_subset$'age:ch1',
cohort = samples_pData_subset$'endotype_cohort:ch1',
d = samples_pData_subset$'pneumonia diagnoses:ch1',
mort = samples_pData_subset$'mortality_event_28days:ch1',
time = as.numeric(samples_pData_subset$'time_to_event_28days:ch1'))

pheno$d <- as.character(pheno$d)
pheno$d[pheno$d == "cap"] <- "cap"
pheno$d[pheno$d == "hap"] <- "cap"
pheno$d[samples_pData_subset$'icu_acquired_infection:ch1' == "healthy"] <- "ctl"
pheno$d[pheno$d == "NA"] <- "other"
pheno$d <- factor(pheno$d, levels=c("ctl", "cap", "other"))
pheno$cens <- 1
pheno$cens[pheno$time == 28] <- 0
pheno$mort_num <- as.numeric(pheno$mort)
pheno$mort_num[is.na(pheno$mort_num)] <- 0
pheno$age <- as.numeric(as.character(pheno$age))


setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/GSE65682_CEL")
files <- list.files()
GSM <- substring(files, 1,10)
files <- files[match(rownames(samples_pData_subset), GSM)]

affyData <- read.celfiles(filenames = files)
affyData <- oligo::rma(affyData, background = T, normalize=T)

# convert probe IDs to gene symbols
affyData <- annotateEset(affyData, hgu219.db)
genes <- featureData(affyData)@data[,c(1,3)]

# gene expression matrix
affyData <- exprs(affyData)
colnames(affyData) <- substring(colnames(affyData), 1,10)
affyData <- merge(na.omit(genes), affyData, by.x=1, by.y='row.names')
affyData <- affyData[duplicated(affyData$SYMBOL) == FALSE,]
rownames(affyData) <- affyData$SYMBOL
affyData <- affyData[,-c(1:2)]
# sum(samples_pData_subset$geo_accession != colnames(affyData))
rm(genes, files, GSM)


# A problem of arrays is that it is difficult to discriminate expressed genes from non-expressed genes.
# Genes with hybridization intensities close to backgroung are expected to be noisy accross all samples.
# Thus, we used mvoutlier sign2 algorithm to detect such genes and remove them.
out <- sign2(affyData, explvar = 0.95, qcrit = 0.9)
affyData_out <- affyData[out$wfinal01 == 0, ]
# Compare
# plot(density(unlist(affyData[out$wfinal01 == 1, ])), xlim = c(1,14))
# lines(density(unlist(affyData[out$wfinal01 == 0, ])), col=2)

rm(out)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("GSE65682_series_all_mortality.Rdat")
