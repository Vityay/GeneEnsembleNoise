library(parallel)
library(caret)
library(xgboost)
library(OptimalCutpoints)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
load("Sepsis_WGCNA.Rdat")

rm(adj, affyData, affyData_SPS, affyData_out, cex1, datExpr, dissTOM,
  dynamicMods, geneTree, minModuleSize, phenotype_SPS, powers,
  samples, samples_pData, samples_pData_subset, sft, softPower, TOM)


pheno <- pheno[pheno$d == "other",]

################################################################################
# Split data into discovery and validation as in GSE65682
inTraining <- which(pheno$cohort == "discovery")

X <- data.frame(datME, age=pheno$age )
Y <- factor(pheno$mort)

# Split data into discovery and validation as in GSE65682
Xtrain0 <- X[inTraining,]
Ytrain0 <- Y[inTraining]
training0 <- data.frame(Y = Ytrain0, Xtrain0)

Xvalid <- X[-inTraining,]
Yvalid <- Y[-inTraining]
validation <- data.frame(Y = Yvalid, Xvalid)

# create xgb.DMatrix
dtrain0 <- xgb.DMatrix(data = as.matrix(training0[,-1]), label = as.numeric(training0$Y)-1)
dvalid <- xgb.DMatrix(data = as.matrix(validation[,-1]), label = as.numeric(validation$Y)-1)

################################################################################
# Hyperparameters grid (see xgboost for details)
params <- expand.grid(
  eta = seq(0.2,0.3,0.02),
  max_depth = c(10:20),
  gamma = seq(0.5, 3, 0.05),
  colsample_bytree = seq(0.1, 1, 0.1),
  min_child_weight = 1,
  subsample = seq(0.1, 1, 0.1)
)

# Find optimal hyperparameters from the grid using 5X cross-validation
cv_xgb <- lapply(seq(nrow(params)), function(i) {
  p <- split(t(params[i,]),
    factor(names(params[i,]),
    levels=c("eta", "max_depth", "gamma", "colsample_bytree", "min_child_weight", "subsample") )
  )
  cv <- xgb.cv(data = dtrain0, params = p, nrounds = 100, nthread = 5, nfold = 5, eval_metric='auc',
                     early_stopping_rounds = 10, objective = "binary:logistic")
  res <- cv$evaluation_log[ cv$best_iteration, 'test_auc_mean'][[1]]
  res
})

# Choose optimal hyperparameters from the grid using 5X cross-validation
params <- params[which.max(unlist(cv_xgb)),]
params <- split(t(params),
  factor(names(params),
  levels=c("eta", "max_depth", "gamma", "colsample_bytree", "min_child_weight", "subsample") )
)

# Train XGB models
XGB_model <- mclapply(c(1:1000), function(i) {

  cat(i, "\r")

  # keep 25% of data from the discovery cohort which will be used for early early_stopping
  # see Feature Engineering and Selection: A Practical Approach for Predictive Models by Max Kuhn and Kjell Johnson
  # 3.4.5 "Validation Sets" for the discussion of this approach
  # http://www.feat.engineering/
  inTrain <- createDataPartition(Ytrain0, p = .75, list = TRUE)$Resample1
  inTest <- seq(Ytrain0)[-inTrain]

  dtrain <- slice(dtrain0, inTrain)
  dtest  <- slice(dtrain0, inTest)

  sink(file="/dev/null")
    # Train the model with optimal hyperparameters (see above) and
    # dtest - 25% of data from the discovery cohort used to determine early stopping
    xgb <- xgb.train(data = dtrain, params = params, nround=100, early_stopping_rounds = 10, nthread = 1,
           watchlist = list(train=dtest), eval_metric='auc', maximize=T)
  sink()

  # predict discovery and validation cohorts (probabilities)
  Ytrain0_prob <- predict(xgb, dtrain0,  ntree=xgb$best_iteration, outputmargin=TRUE)
  Yvalid_prob <- predict(xgb, dvalid,  ntree=xgb$best_iteration, outputmargin=TRUE)

  # Find optimal cutpoint from the ROC curve maximizing the product of Sensitivity and Specificity
  df0 <- data.frame(y=Ytrain0_prob, x = as.numeric(Ytrain0) -1)
  p_cut <- optimal.cutpoints(y~x, status=Ytrain0, tag.healthy = 0, data= df0, methods="MaxProdSpSe",
    ci.fit=FALSE, control = control.cutpoints(ci.SeSp = "AgrestiCoull", ci.PV = "AgrestiCoull") )
  cutoff <- p_cut$MaxProdSpSe$Global$optimal.cutoff$cutoff[1]

  # predict discovery and validation cohorts (classes)
  Ytrain0_class <- factor(as.numeric(ifelse(Ytrain0_prob <= cutoff, "0", "1")), levels = c(0,1))
  Yvalid_class  <- factor(as.numeric(ifelse(Yvalid_prob  <= cutoff, "0", "1")), levels = c(0,1))

  # ballanced accuracy for the discovery and validation cohorts
  bACC_train  <- confusionMatrix(Ytrain0_class, factor(as.numeric(Ytrain0) -1), "1")$byClass[[11]]
  bACC_valid  <- confusionMatrix(Yvalid_class, factor(as.numeric(Yvalid) -1), "1")$byClass[[11]]

  res <- list(model = xgb, auc_train=bACC_train, auc_valid = bACC_valid, nround = xgb$best_iteration, params = params)

  res

}, mc.cores=6 )

# Choose the model which generalize best the validation cohort
idx <- which.max(sapply(XGB_model, function(x) x$auc_valid))
XGB_model_bst <- XGB_model[[idx]]$model
# keep ntree of the selected model
XGB_model_ntree <- XGB_model[[idx]]$nround
# XGB_model_params <- XGB_model[[idx]]$params

rm(idx,cv_xgb, XGB_model, params, dtrain0, dvalid)

# Save results
setwd("/media/yuri/Data/home1/GNI_data/infection/ModelsSciRep")
save.image('Sepsis_WGCNA_XGB.Rdat')
xgb.save(XGB_model_bst, 'Sepsis_WGCNA_XGB.model')
