# Gamlss estimation of gene noise for each gene for COVID19

library(GEOquery)
library(edgeR)

setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")

# read Series Matrix File
# samples <- getGEO(file = "GSE157103_series_matrix.txt.gz")
# OR, from the url https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103
samples <- getGEO("GSE157103")[[1]]

# format	Series Matrix File
samples_pData <-  pData(phenoData(samples))
samples_pData <-  samples_pData[,c(
'description', 'Sex:ch1', 'age (years):ch1', 'disease state:ch1',
'hospital-free days post 45 day followup (days):ch1', 'hospital-free days post 45 day followup (days):ch1', 'icu:ch1',
'apacheii:ch1', 'charlson score:ch1', 'sofa:ch1',  'dm:ch1', 'ventilator-free days:ch1',
'crp (mg/l):ch1', 'ddimer (mg/l_feu):ch1', 'ferritin (ng/ml):ch1', 'fibrinogen:ch1', 'lactate (mmol/l):ch1',
'procalcitonin (ng/ml):ch1'
 )]

colnames(samples_pData) <- c(
"samples", "Sex", "age", "disease", "HFD45", "stage", "icu",
"apacheii", "charlson", "sofa", "dm", "VFD",
"crp", "ddimer", "ferritin", "fibrinogen", "lactate", "procalcitonin"
)
samples_pData$age[samples_pData$age == ">89"] <- "90"
samples_pData[,c(3,8:18)] <- apply(samples_pData[,c(3,8:18)], 2, as.numeric)

samples_pData$disease[samples_pData$disease == "COVID-19"] <- 1
samples_pData$disease[samples_pData$disease == "non-COVID-19"] <- 0
samples_pData$disease <- factor(samples_pData$disease)

samples_pData$HFD45 <- as.numeric(samples_pData$HFD45)
samples_pData$stage <- "severe"
samples_pData$stage[samples_pData$HFD45 >= 26] <- "mild"
samples_pData$stage[samples_pData$disease == 0] <- "non-covid"


################################################################################
################################################################################
################################################################################
# read and process RNA counts from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103

# read RNA counts from the file
# RNAcounts <- read.table("GSE157103_genes.ec.tsv.gz", header = T, comment.char = "", sep = "\t", row.names="X.symbol")
# OR, read RNA counts from the url
options(timeout = max(60*60, getOption("timeout")))
getGEOSuppFiles("GSE157103")
setwd("GSE157103/")
RNAcounts <- read.table("GSE157103_genes.ec.tsv.gz", header = T, comment.char = "", sep = "\t", row.names="X.symbol")

RNAcounts <- round(RNAcounts)
# sum(colnames(RNAcounts) != samples_pData$description)

# select only samples from COVID19 patients (non-COVID19 patients are also ICU patients, but cause of pathology is not available)
idx <- samples_pData$disease == "1"
samples_pData <- samples_pData[idx,]
RNAcounts <- RNAcounts[,idx]
rm(idx)

# For convinience put RNA counts to DGEList
DGE_RNAcounts <- DGEList(RNAcounts, samples=samples_pData )
DGE_RNAcounts <- calcNormFactors(DGE_RNAcounts)
DGE_RNAcounts$samples$offset <- log(DGE_RNAcounts$samples$lib.size*DGE_RNAcounts$samples$norm.factors)
DGE_RNAcounts$counts <- DGE_RNAcounts$counts[apply(DGE_RNAcounts$counts, 1, function(x) mean(x==0)) <= 0.05, ]

rm(RNAcounts)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("COVID19_RNAseq.Rdat")
################################################################################

# Prepare log2 CPM
library(edgeR)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("COVID19_RNAseq.Rdat")

idx <- apply(DGE_RNAcounts$counts, 1, function(x) all(x > 0) )
COVID_logCPM <- log2(sweep(DGE_RNAcounts$counts[idx,], 2, 1e-06*DGE_RNAcounts$samples$lib.size, "/"))
COVID_pheno <- DGE_RNAcounts$samples[,c("samples", "Sex", "age", "stage", "HFD45")]

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("COVID_logCPM.Rdat")
