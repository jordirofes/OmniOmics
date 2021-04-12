setGeneric("mlFit", function(dt, ...){
    standardGeneric("mlFit")
})
setMethod("mlFit", definition = function(dt, groupvar, metric = "Accuracy", cvmethod = "loo", kfolds = 10, ntimes = 5, preproc = c("center", "scale"), tlength = 20){
    dt2 <- t(extractData(dt))
    grp <- extractPhenoData(dt)[[groupvar]]
    if(cvmethod == "loo"){
        kfolds <- length(grp)
    }
    folds <- createMultiFolds(grp, k = kfolds, times = ntimes)
    control <- trainControl("repeatedcv", index = folds,
                            selectionFunction = "oneSE")
    mod <- train(dt2, grp, method = method, metric = metric,
                tuneLength = tlength, trControl = control,
                preProcess = preproc)

    # rf_mod <- train(dt2, grp,method = "rf", metric = "Accuracy",
    #                 tuneLength = 20, trControl = control,
    #                 preProcess = c("zv", "center", "scale"))
    # svm_mod <- train(dt2, grp,method = "svmLinear", metric = "Accuracy",
    #                 tuneLength = 20, trControl = control,
    #                 preProcess = c("zv", "center", "scale"))
})

mlPredictCM <- function(mlmod, newdt, groupvar, posclass){
    newdt2 <- extractData(newdt)
    pred_dt <- predict(mlmod, newdt)
    newdt_labs <- extractPhenoData(newdt)
    cm <- confusionMatrix(table(pred_dt, newdt_labs), positive = posclass)
    return(cm)
}

mlPredictROC <- function(mlmod, newdt, groupvar, posclass){
    newdt2 <- extractData(newdt)
    pred_dt <- predict(mlmod, newdt)
    newdt_labs <- extractPhenoData(newdt)
    dt_class <- unique(newdt_labs)
    idx <- which(dt_class == posclass)
    idx2 <- which(dt_class != posclass)
    pred <- prediction(mlmod, labels = newdt_labs, label.ordering = dt_class[c(idx2, idx)])
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    perf2 <- performance(pred2, measure = "auc")
    auc <- unlist(perf2@y.values)
    plot(perf, main = paste("ROC curve \n AUC = ", round(auc,2)))
}
