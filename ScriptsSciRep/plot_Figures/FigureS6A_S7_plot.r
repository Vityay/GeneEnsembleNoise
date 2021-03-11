library(parallel)
library(caret)
library(MASS)
library(xgboost)

library(OptimalCutpoints)
library(pROC)
library(ROCR)

library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(gplots)

setwd("/media/yuri/Data/home1/GNI_data/infection/ModelsSciRep")
load('CAP_sepsis_WGCNA_XGB.Rdat')
XGB_model_bst <- xgb.load('CAP_sepsis_WGCNA_XGB.model')

dtrain0 <- xgb.DMatrix(data = as.matrix(training0[,-1]), label = as.numeric(training0$Y)-1)
dvalid <- xgb.DMatrix(data = as.matrix(validation[,-1]), label = as.numeric(validation$Y)-1)

# predict discovery and validation cohorts (probs)
Ytrain_prob <- predict(XGB_model_bst, dtrain0, ntree = XGB_model_ntree, outputmargin=TRUE)
Yvalid_prob <- predict(XGB_model_bst, dvalid,  ntree = XGB_model_ntree, outputmargin=TRUE)

################################################################################
# Find optimal cutpoint from the ROC curve maximizing the product of Sensitivity and Specificity
p_cut <- optimal.cutpoints(y~x, status=Ytrain0, tag.healthy = 0, data= data.frame(y=Ytrain_prob, x = as.numeric(Ytrain0) -1),
methods="MaxProdSpSe", ci.fit=TRUE, control = control.cutpoints(ci.SeSp = "AgrestiCoull", ci.PV = "AgrestiCoull") )
cutoff <- p_cut$MaxProdSpSe$Global$optimal.cutoff$cutoff[1]

# predict discovery and validation cohorts (classes)
p_train_class <- factor(as.numeric(ifelse(Ytrain_prob <= cutoff, 0, 1)), levels = c(0,1))
p_valid_class  <- factor(as.numeric(ifelse(Yvalid_prob  <= cutoff, 0, 1)), levels = c(0,1))

# Model summary bACC - balanced accuracy, Sensitivity, Specificity
conf_train <- confusionMatrix(p_train_class, Ytrain0, "0")
conf_valid  <- confusionMatrix(p_valid_class, Yvalid, "0")

conf_table <- data.frame( bACC = c( conf_train$byClass[[11]], conf_valid$byClass[[11]] ),
Sensitivity = c( conf_train$byClass[[1]], conf_valid$byClass[[1]] ),
Specificity = c( conf_train$byClass[[2]], conf_valid$byClass[[2]] )
)
conf_table <- t(conf_table)
colnames(conf_table) <- c("discovery", "validation")
conf_table <- round(conf_table, 3)
################################################################################
# Variable importance
var_imp <- xgb.importance(dimnames(dtrain0)[[2]], model = XGB_model_bst)
# var_imp <- var_imp[ order(as.numeric(gsub("x", "", var_imp$Feature))),]
# var_imp$Feature <- ids[as.numeric(gsub("x", "", var_imp$Feature))]
var_imp <- var_imp[order(var_imp$Gain, decreasing=F),]
var_imp <- as.data.frame(var_imp)
rownames(var_imp) <- gsub("ME", "", var_imp[,1])
var_imp <- var_imp[,-1]

colors <- rownames(var_imp)
colors[colors == "age"] <- "black"

################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures/")
pdf("FigureS7A.pdf", width=4, height=9)

k <- which(colors == "grey")+1
m <- length(colors)
par(mar=c(4,16,1,0.5))
barplot(var_imp[k:m,1], horiz = T, names.arg = rownames(var_imp)[k:m], col = col2hex(colors)[k:m],
      xlab = "Gain", main = "variable importance",las=1, cex.axis=0.7, cex.names = 0.7)

dev.off()
################################################################################

Train_roc <- roc(Ytrain0, Ytrain_prob)
Train_auc <- lapply(pROC::ci(Ytrain0, Ytrain_prob), function(x) round(x, 3))

Valid_roc <- roc(Yvalid, Yvalid_prob)
Valid_auc <- lapply(pROC::ci(Yvalid, Yvalid_prob), function(x) round(x, 3))

################################################################################
pdf("FigureS6A.pdf", width=6, height=6)

plot(Train_roc, lty = 1, lwd = 2, col=4, type="s")
plot(Valid_roc, add=TRUE, lty = 1, lwd = 2, col=2, type="s")
legend("bottomright", legend =
  c(paste("discovery cohort: AUC = ", Train_auc[[2]], " (", Train_auc[[1]], ", ", Train_auc[[3]], ")", sep=""),
    paste("validation cohort: AUC = ", Valid_auc[[2]], " (", Valid_auc[[1]], ", ", Valid_auc[[3]], ")", sep="")),
lty=1, lwd=2, col=c(4,2), bty="n")

dev.off()
################################################################################
# enrichment


setwd("/media/yuri/Data/home1/GNI_data/infection/DataSciRep/")
load("org.Hs.core_paths.Rdat")

names(kegg_path) <- paste0("kegg:", names(kegg_path))
names(corum_core) <- paste0("complex:", names(corum_core))
corum_core$`complex:mitochondrial respiratory chain complex I assembly` <- NULL
core_paths <- c(corum_core, kegg_path)

topME <- rev(colors)
topME <- topME[topME!="black"]

# For explanetion of gene set enrichment see
# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
N <- length(unlist(core_paths)[!duplicated(unlist(core_paths))])
m <- sapply(core_paths, length)
n <- N-m

topME_ens <- lapply(topME, function(module) {

  k <- sum(GeneNames[dynamicColors %in% module] %in% unlist(core_paths)[!duplicated(unlist(core_paths))])
  x <- sapply(core_paths, function(x) sum(x %in% GeneNames[dynamicColors %in% module]) )

  sig <- phyper(q = x -1, m=m, n=n, k=k, lower.tail=FALSE)
  sig <- p.adjust(sig, "fdr")
  res <- sig[sig < 1e-03]

  if(length(res) > 0) {
    res <- data.frame(module = module, ensemble = names(res), fdr = -log10(res))
    res <- res[order(res$fdr, decreasing = T),]
    rownames(res) <- NULL
  } else {res <- NULL}

  res

})

topME_ens <- do.call(rbind, topME_ens)
topME_ens <- topME_ens[!duplicated(topME_ens$ensemble),]
topME_ens <- topME_ens[nrow(topME_ens):1,]
topME_ens$fdr[topME_ens$fdr > -log10(1e-10)] <- -log10(1e-10)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures/")
pdf("FigureS7B.pdf", width=4, height=9)
  k <- max(which(topME_ens[,'module'] == "grey"))+1
  m <- length(topME_ens[,'module'])
  par(mar=c(4,16,1,0.5))
  barplot(topME_ens[k:m,'fdr'], horiz = T, names.arg = topME_ens[k:m,'ensemble'], col = col2hex(topME_ens[k:m,'module']),
        xlab = "-log10(fdr)", main = "WGCNA enrichment (fdr < 0.001)",las=1, cex.axis=0.7, cex.names = 0.7)
dev.off()
