# Script to process illumina gene expression dataset GSE21802
# "Hosts responses in critical disease caused by pandemic H1N1"
# H1N1 data

library(Biobase)
library(GEOquery)
library(illuminaHumanv2.db)
library(beadarray)
library(blima)
library(xtable)
library(data.table)
library(preprocessCore)

# convert illumina IDs to gene symbol
illuminaToSymbol = toTable(illuminaHumanv2SYMBOLREANNOTATED)

setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")

options(timeout = max(60*60, getOption("timeout")))
getGEOSuppFiles("GSE21802")
setwd("GSE21802/")

# Illumina detection p values will be used to remove genes with low expression/low detection
pVal <- read.table("GSE21802_non-normalized_data.txt.gz", header=T, sep="\t")
ILMN <- pVal$ID_REF
pVal <- pVal[,grep("Pval", colnames(pVal))]
rownames(pVal) <- ILMN
# FDR adjusted illumina detection p values (FDR < 0.1 - significant detection)
pVal <- apply(pVal, 2, function(x) p.adjust(x, "BH"))
pVal[pVal <= 0.1] <- 0
pVal[pVal > 0.1] <- 1
pVal <- abs(pVal-1)
pVal <- pVal[rowMeans(pVal) == 1,]
pVal <- rownames(pVal)

# Select genes with significant detection
illuminaToSymbol <- illuminaToSymbol[match(pVal, illuminaToSymbol$IlluminaID),]

# get GSE21802 sample and expression data
gset <- getGEO("GSE21802", destdir=".", GSEMatrix =TRUE, AnnotGPL = TRUE, getGPL = TRUE)
gset <- gset[[1]]

# get GSE21802 sample data
samples_H1N1 <-  pData(phenoData(gset))
samples_H1N1$'patient:ch1'[samples_H1N1$'patient:ch1' == "control"] <- paste0("CTL", which(samples_H1N1$'patient:ch1' == "control"))

# get GSE21802 expression data and select genes with significant detection
eset_H1N1 = log2(exprs(gset))
eset_H1N1 <- eset_H1N1[na.omit(match(illuminaToSymbol$IlluminaID, rownames(eset_H1N1))),]
illuminaToSymbol <- illuminaToSymbol[na.omit(match(rownames(eset_H1N1), illuminaToSymbol$IlluminaID)),]
# sum(rownames(eset_H1N1) != illuminaToSymbol$IlluminaID)

# gene expression matrix
df <- data.table(Gene = illuminaToSymbol$SymbolReannotated, eset_H1N1)
df <- df[ ,lapply(.SD, mean), by = Gene]
Genes <- df[,1][[1]]
df <- as.matrix(df[,-1])

# average technical replicates
idx <- split(c(1:nrow(samples_H1N1)), samples_H1N1$'patient:ch1')
idx2 <- sapply(idx, function(x) x[[1]])
df <- Map(function(i) if(length(i) > 1) { res <- rowMeans(df[,i]) } else { res <- df[,i] } , i = idx )
df <- do.call(cbind, df)
colnames(df) <- samples_H1N1$geo_accession[idx2]
rownames(df) <- Genes

# finalize
samples_H1N1 <- samples_H1N1[idx2,] # samples data.frame
eset_H1N1 <- df # expression data.frame
rm(gset, df, idx, idx2, illuminaToSymbol, Genes, pVal, ILMN)

eset_H1N1_m <- as.matrix(eset_H1N1) # expression data.frame to matrix
eset_H1N1_m <- normalize.quantiles(eset_H1N1_m) # normalize expression data
rownames(eset_H1N1_m) <- rownames(eset_H1N1)
colnames(eset_H1N1_m) <- colnames(eset_H1N1)

eset_H1N1 <- eset_H1N1_m; rm(eset_H1N1_m) # save normalized expression data to eset_H1N1

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("eset_H1N1.Rdat")
