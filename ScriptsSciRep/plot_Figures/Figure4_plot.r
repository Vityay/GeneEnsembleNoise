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
load('CAP_Ensemble.Rdat')
XGB_model_bst <- xgb.load('CAP_Ensemble.model')

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
conf_train <- confusionMatrix(p_train_class, Ytrain0, "1")
conf_valid  <- confusionMatrix(p_valid_class, Yvalid, "1")

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
var_imp <- var_imp[ order(as.numeric(gsub("x", "", var_imp$Feature))),]
var_imp$Feature <- ids[as.numeric(gsub("x", "", var_imp$Feature))]
var_imp <- var_imp[order(var_imp$Gain, decreasing=F),]
var_imp <- var_imp[-grep("Pr.cl.", var_imp$Feature),]
var_imp <- as.data.frame(var_imp)
rownames(var_imp) <- var_imp[,1]
var_imp <- var_imp[,-1]

################################################################################
# ROC analysis
Train_roc <- roc(Ytrain0, Ytrain_prob)
Train_auc <- lapply(pROC::ci(Ytrain0, Ytrain_prob), function(x) round(x, 3))

Valid_roc <- roc(Yvalid, Yvalid_prob)
Valid_auc <- lapply(pROC::ci(Yvalid, Yvalid_prob), function(x) round(x, 3))

################################################################################
setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/Figures/")
pdf("Figure4D.pdf", width=6, height=3.5)

par(mar=c(4,16,1,0.5))
barplot(var_imp[,1], horiz = T, names.arg = rownames(var_imp),
      xlab = "Gain", main = "variable importance",las=1, cex.axis=0.7, cex.names = 0.7)

dev.off()
################################################################################
pdf("Figure4B.pdf", width=6, height=6)

plot(Train_roc, lty = 1, lwd = 2, col=4, type="s")
plot(Valid_roc, add=TRUE, lty = 1, lwd = 2, col=2, type="s")
legend("bottomright", legend =
  c(paste("discovery cohort: AUC = ", Train_auc[[2]], " (", Train_auc[[1]], ", ", Train_auc[[3]], ")", sep=""),
    paste("validation cohort: AUC = ", Valid_auc[[2]], " (", Valid_auc[[1]], ", ", Valid_auc[[3]], ")", sep="")),
lty=1, lwd=2, col=c(2,4), bty="n")

dev.off()
################################################################################
pdf("Figure4A.pdf", width=6, height=6)

par(mfrow=c(1,2))
y.lim = range(c(Ytrain_prob, Yvalid_prob))
boxplot(Ytrain_prob ~ Ytrain0, pch="",
        col="#00000022", border=c("#4591CC", "#21A883"), lty=1, las=1, ylab = "model score (discovery)", xlab = NA, ylim=y.lim)
stripchart(Ytrain_prob ~ Ytrain0, vertical=T, method="jitter", pch=c(15,16), col=c("#4591CCAA", "#21A883AA"), add=T)
abline(h = cutoff, lty=2)

boxplot(Yvalid_prob ~ Yvalid, pch="", col="#00000022", border=c("#4591CC", "#21A883"), lty=1, las=1,
        ylab = "model score (validation)", xlab = NA, ylim=y.lim)
stripchart(Yvalid_prob ~ Yvalid, vertical=T, method="jitter", pch=c(15,16), col=c("#4591CCAA", "#21A883AA"), add=T)
abline(h = cutoff, lty=2)

dev.off()
################################################################################
ggsave(filename="Table2_CAP_Ensemble_XGB.pdf", grid.arrange(tableGrob(conf_table)),
        width = 6, height = 4, units = "cm")
################################################################################

pheno_train <- pheno_CAP[inTraining,]
pheno_valid <- pheno_CAP[-inTraining,]

surv_df_train <- data.frame(time = pheno_train$time, cens = pheno_train$cens, pred_cl = p_train_class, true_cl = Ytrain0)
surv_train <- survfit(Surv(time, cens) ~ pred_cl, data=surv_df_train)

surv_df_valid <- data.frame(time = pheno_valid$time, cens = pheno_valid$cens, pred_cl = p_valid_class, true_cl = Yvalid)
surv_valid <- survfit(Surv(time, cens) ~ pred_cl, data=surv_df_valid)

p1 <- ggsurvplot(surv_train, data = surv_df_train, linetype = "strata", conf.int = TRUE, pval = TRUE, palette= c("#4591CC", "#21A883"))
p2 <- ggsurvplot(surv_valid, data = surv_df_valid, linetype = "strata", conf.int = TRUE, pval = TRUE, palette= c("#4591CC", "#21A883"))

pdf("Figure4C.pdf", paper="a4", width=3.5, height=3.7)
print(p1, newpage=F)
print(p2, newpage=T)
dev.off()
