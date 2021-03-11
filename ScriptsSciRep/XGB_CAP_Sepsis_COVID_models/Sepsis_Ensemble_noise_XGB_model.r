library(parallel)
library(caret)
library(xgboost)
library(OptimalCutpoints)

setwd("/media/yuri/Data/home1/GNI_data/infection/ResultsSciRep/")
load("CAP_Sepsis_Ensemble_noise_model.RDat")

rm(affyDat_noise_paths, affyDat_noise_paths_CAP, pheno, pheno_CAP)
################################################################################

inTraining <- which(pheno_SPS$cohort == "discovery")

features <- sapply(seq(ncol(affyDat_noise_paths_SPS)), function(i)
    t.test(affyDat_noise_paths_SPS[inTraining,i] ~ pheno_SPS$mort[inTraining])$p.value )


X <- data.frame(affyDat_noise_paths_SPS[, features <= 0.01], pheno=pheno_SPS$age )
ids <- colnames(X)
colnames(X) <- paste0("x", seq(ncol(X)))
Y <- factor(pheno_SPS$mort)

Xtrain0 <- X[inTraining,]
Ytrain0 <- Y[inTraining]
training0 <- data.frame(Y = Ytrain0, Xtrain0)

Xvalid <- X[-inTraining,]
Yvalid <- Y[-inTraining]
validation <- data.frame(Y = Yvalid, Xvalid)

dtrain0 <- xgb.DMatrix(data = as.matrix(training0[,-1]), label = as.numeric(training0$Y)-1)
dvalid <- xgb.DMatrix(data = as.matrix(validation[,-1]), label = as.numeric(validation$Y)-1)

################################################################################

params <- expand.grid(
  eta = seq(0.15,0.25,0.01),
  max_depth = c(5:20),
  gamma = seq(0.005, 0.01, 0.001),
  colsample_bytree = seq(0.1, 1, 0.1),
  min_child_weight = 1,
  subsample = seq(0.1, 1, 0.1)
)

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

params <- params[which.max(unlist(cv_xgb)),]
params <- split(t(params),
  factor(names(params),
  levels=c("eta", "max_depth", "gamma", "colsample_bytree", "min_child_weight", "subsample") )
)

XGB_model <- mclapply(c(1:1000), function(i) {

  cat(i, "\r")

  inTrain <- createDataPartition(Ytrain0, p = .75, list = TRUE)$Resample1
  inTest <- seq(Ytrain0)[-inTrain]

  dtrain <- slice(dtrain0, inTrain)
  dtest  <- slice(dtrain0, inTest)

  sink(file="/dev/null")

    xgb <- xgb.train(data = dtrain, params = params, nround=100, early_stopping_rounds = 10, nthread = 1,
           watchlist = list(train=dtest), eval_metric='auc', maximize=T)
  sink()


  Ytrain0_prob <- predict(xgb, dtrain0,  ntree=xgb$best_iteration, outputmargin=TRUE)
  Yvalid_prob <- predict(xgb, dvalid,  ntree=xgb$best_iteration, outputmargin=TRUE)

  df0 <- data.frame(y=Ytrain0_prob, x = as.numeric(Ytrain0) -1)
  p_cut <- optimal.cutpoints(y~x, status=Ytrain0, tag.healthy = 0, data= df0, methods="MaxProdSpSe",
    ci.fit=FALSE, control = control.cutpoints(ci.SeSp = "AgrestiCoull", ci.PV = "AgrestiCoull") )
  cutoff <- p_cut$MaxProdSpSe$Global$optimal.cutoff$cutoff[1]

  Ytrain0_class <- factor(as.numeric(ifelse(Ytrain0_prob <= cutoff, "0", "1")), levels = c(0,1))
  Yvalid_class  <- factor(as.numeric(ifelse(Yvalid_prob  <= cutoff, "0", "1")), levels = c(0,1))

  bACC_train  <- confusionMatrix(Ytrain0_class, factor(as.numeric(Ytrain0) -1), "1")$byClass[[11]]
  bACC_valid  <- confusionMatrix(Yvalid_class, factor(as.numeric(Yvalid) -1), "1")$byClass[[11]]

  res <- list(model = xgb, auc_train=bACC_train, auc_valid = bACC_valid, nround = xgb$best_iteration, params = params)

  res

}, mc.cores=6 )

idx <- which.max(sapply(XGB_model, function(x) x$auc_valid))
XGB_model_bst <- XGB_model[[idx]]$model
XGB_model_ntree <- XGB_model[[idx]]$nround
# XGB_model_params <- XGB_model[[idx]]$params

rm(idx,cv_xgb, XGB_model, params, dtrain0, dvalid)

setwd("/media/yuri/Data/home1/GNI_data/infection/ModelsSciRep")
save.image('Sepsis_Ensemble.Rdat')
xgb.save(XGB_model_bst, 'Sepsis_Ensemble.model')
