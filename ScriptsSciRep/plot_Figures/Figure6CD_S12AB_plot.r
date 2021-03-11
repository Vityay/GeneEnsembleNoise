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

library(survival)
library(survminer)

setwd("/media/yuri/Data/home1/GNI_data/infection/ModelsSciRep")
load('COVID19_Ensemble.Rdat')
XGB_model_bst <- xgb.load('COVID19_Ensemble.model')
XGB_model_ntree <- XGB_model_bst$best_iteration

# create xgb.DMatrix
dtrain0 <- xgb.DMatrix(data = as.matrix(trn_x[,-1]), label = trn_y)
dvalid  <- xgb.DMatrix(data = as.matrix(val_x[,-1]), label = val_y)

# predict discovery and validation cohorts (probs)
Ytrain_prob <- predict(XGB_model_bst, dtrain0, ntree = XGB_model_ntree, outputmargin=TRUE)
Yvalid_prob <- predict(XGB_model_bst, val_x[,-1],  ntree = XGB_model_ntree, outputmargin=TRUE)

################################################################################
# Find optimal cutpoint from the ROC curve maximizing the product of Sensitivity and Specificity
p_cut <- optimal.cutpoints(y~x, status=trn_y, tag.healthy = 0, data= data.frame(y=Ytrain_prob, x = trn_y),
methods="MaxProdSpSe", ci.fit=TRUE, control = control.cutpoints(ci.SeSp = "AgrestiCoull", ci.PV = "AgrestiCoull") )
cutoff <- p_cut$MaxProdSpSe$Global$optimal.cutoff$cutoff[1]

# predict discovery and validation cohorts (classes)
p_train_class  <- factor(as.numeric(ifelse(Ytrain_prob <= cutoff, 0, 1)), levels = c(0,1))
p_valid_class  <- factor(as.numeric(ifelse(Yvalid_prob  <= cutoff, 0, 1)), levels = c(0,1))

# Model summary bACC - balanced accuracy, Sensitivity, Specificity
conf_train  <- confusionMatrix(p_train_class, factor(trn_y), "1")
conf_valid  <- confusionMatrix(p_valid_class, factor(val_y), "1")

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
var_imp <- var_imp[order(var_imp$Gain, decreasing=T),][c(1,3,5,7,8,10,11,13:15),]
var_imp <- as.data.frame(var_imp)
rownames(var_imp) <- var_imp[,1]
var_imp <- var_imp[,-1]

var_imp <- var_imp[order(var_imp$Gain, decreasing=F),]

################################################################################
# ROC analysis
Train_roc <- roc(trn_y, Ytrain_prob)
Train_auc <- lapply(pROC::ci(trn_y, Ytrain_prob), function(x) round(x, 3))

Valid_roc <- roc(val_y, Yvalid_prob)
Valid_auc <- lapply(pROC::ci(val_y, Yvalid_prob), function(x) round(x, 3))

All_roc <- roc(c(trn_y, val_y), c(Ytrain_prob, Yvalid_prob))
All_auc <- lapply(pROC::ci(c(trn_y, val_y), c(Ytrain_prob, Yvalid_prob)), function(x) round(x, 3))


################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures/")
pdf("Figure6D.pdf", width=6, height=3.5)

par(mar=c(4,16,1,0.5))
barplot(var_imp[,1], horiz = T, names.arg = rownames(var_imp),
      xlab = "Gain", main = "variable importance",las=1, cex.axis=0.7, cex.names = 0.7)

dev.off()
################################################################################
pdf("Figure6C.pdf", width=6, height=6)

plot(Train_roc, lty = 1, lwd = 2, col=4, type="s")
plot(Valid_roc, add=TRUE, lty = 1, lwd = 2, col=2, type="s")
legend("bottomright", legend =
  c(paste("discovery cohort: AUC = ", Train_auc[[2]], " (", Train_auc[[1]], ", ", Train_auc[[3]], ")", sep=""),
    paste("validation cohort: AUC = ", Valid_auc[[2]], " (", Valid_auc[[1]], ", ", Valid_auc[[3]], ")", sep="")),
lty=1, lwd=2, col=c(2,4), bty="n")

dev.off()
################################################################################


################################################################################

trn_res <- data.frame(trueY = trn_y, Yhat = p_train_class, Yp = Ytrain_prob,
  Acc = as.numeric(trn_y == p_train_class), hfd = COVID_pheno$HFD45[train_idx],
  cohort = "discovery")
val_res <- data.frame(trueY = val_y, Yhat = p_valid_class, Yp = Yvalid_prob,
  Acc = as.numeric(val_y == p_valid_class), hfd = COVID_pheno$HFD45[test_idx],
  cohort = "validation")

plot_res <- rbind(trn_res, val_res)
plot_res$shape <- paste0(plot_res$cohort, plot_res$Acc)

p1 <- ggplot(plot_res, aes(y=Yp, x = as.factor(trueY) )) + geom_boxplot(aes(color = as.factor(trueY)), outlier.size=-1) +
geom_hline(yintercept=cutoff) +
geom_jitter(aes(color = as.factor(trueY), shape = shape), size=3, width = 0.1) +
scale_shape_manual(values = c(3,2,4,17)) +
scale_color_manual(values = c("#4591cc", "#21a883")) +
xlab("") + ylab("model score") +
theme_bw() +
theme(text = element_text(size=12),
        axis.text.y = element_text(size=12))

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures/")
ggsave("FigureS12Aleft.pdf", p1, width=4, height=6)
