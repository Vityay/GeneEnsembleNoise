library(WGCNA)
library("flashClust")
library("doParallel")
# options(stringsAsFactors = FALSE);

### CAP_sepsis ###
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
load("GSE65682_series_all_mortality.Rdat")
affyData_CAP <- affyData_out[,pheno$d == "cap"]
phenotype_CAP <- pheno$mort[pheno$d == "cap"]

datExpr = t(affyData_CAP)
GeneNames = rownames(affyData_CAP)

# Calculate a soft-treshold power
powers = c(c(1:10), seq(from = 12, to=30, by=2));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,
  corFnc = cor,corOptions = list(use = 'p'), networkType = "unsigned")

################################################################################
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
################################################################################
# Choose a soft-treshold power

softPower = 12;

#calculate the adjacency matrix
adj <- adjacency(datExpr,type = "unsigned", power = softPower);

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
# TOM <- TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
TOM <- TOMsimilarity(adj,  TOMType = "unsigned")
rownames(TOM) <- GeneNames
colnames(TOM) <- GeneNames
dissTOM <- 1-TOM

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 10;

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

datME = moduleEigengenes(datExpr,dynamicColors)$eigengenes

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep")
save.image("CAP_WGCNA.Rdat")
