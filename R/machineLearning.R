#'@title Caret model fitting
#'@author Jordi Rofes Herrera
#'@description Trains selected machine learning model from features
#'@param dt A SummarizedExperiment or ExpressionSet object
#'@param filenum An optional numeric vector indicating the files to plot
#'@param groupvar A numeric or string indicating the variable from the phenodata
#'to use for grouping (must be a factor of length two)..
#'@param metric A string with the metric to compute for each training step (see ?train for more details)
#'@param cvmethod A string indicating de crossvalidation method "loo" for leave one out.
#'@param kfold Number of folds for k-fold crossvalidation
#'@param ntimes Number of partition to create
#'@param preproc A string vector with pre-processing options. See ?train for more details.
#'@param tlenght The granularity of the training parameters
#'@return A train object
#'@export
setGeneric("mlFit", function(dt, groupvar, method, metric = "Accuracy",
                                cvmethod = "loo", kfolds = 10, ntimes = 5,
                                preproc = c("center", "scale"), tlength = 20,
                                ...){
    standardGeneric("mlFit")
})
#'@export
setMethod("mlFit",
            definition = function(dt, groupvar, method, metric = "Accuracy",
                                cvmethod = "loo", kfolds = 10, ntimes = 5,
                                preproc = c("center", "scale"), tlength = 20){
    dt2 <- as.data.frame(t(extractData(dt)))
    grp <- extractPhenoData(dt)[[groupvar]]
    if(cvmethod == "loo"){
        kfolds <- length(grp)
    }
    folds <- createMultiFolds(grp, k = kfolds, times = ntimes)
    control <- trainControl("repeatedcv", index = folds,
                            selectionFunction = "oneSE")
    mod <- train(dt2, grp, method = method, metric = metric,
                tuneLength = tlength, trControl = control,
                preProcess = preproc, )
    # rf_mod <- train(dt2, grp,method = "rf", metric = "Accuracy",
    #                 tuneLength = 20, trControl = control,
    #                 preProcess = c("zv", "center", "scale"))
    # svm_mod <- train(dt2, grp,method = "svmLinear", metric = "Accuracy",
    #                 tuneLength = 20, trControl = control,
    #                 preProcess = c("zv", "center", "scale"))
})
#'@title Caret confusion matrix
#'@author Jordi Rofes Herrera
#'@description Calculates the caret confusion matrix with the machine learning model and new data.
#'@param newdt A SummarizedExperiment or ExpressionSet object with the test data
#'@param mlmod A machine learning model with a predict method.
#'@param prepro_obj An optional pre-processing object created by the preProcess() function. Not required with a training object.
#'@param groupvar A numeric or string indicating the variable from the phenodata
#'to use for grouping (must be a factor of length two).
#'@param posclass A string indicating the positive class from the groupvar
#'@return A confusion matrix
#'@export
mlPredictCM <- function(mlmod, newdt, prepro_obj, groupvar, posclass){
    newdt2 <- t(extractData(newdt))
    if(!missing(prepro_obj)){
        newdt2 <- predict(prepro_obj, newdt2)
    }
    pred_dt <- predict(mlmod, newdt2)
    newdt_labs <- extractPhenoData(newdt)[[groupvar]]
    cm <- confusionMatrix(table(pred_dt, newdt_labs), positive = posclass)
    return(cm)
}
#'@title Receiver Operand Curve (ROC)
#'@author Jordi Rofes Herrera
#'@description Calculates the ROC curve with the machine learning model and new data.
#'@param newdt A SummarizedExperiment or ExpressionSet object with the test data
#'@param mlmod A machine learning model with a predict method with "prob" option.
#'@param groupvar A numeric or string indicating the variable from the phenodata
#'to use for grouping (must be a factor of length two).
#'@param posclass A string indicating the positive class from the groupvar
#'@return A ROC curve plot
#'@export
mlPredictROC <- function(mlmod, newdt, prepro_obj, groupvar, posclass){
    newdt2 <- t(extractData(newdt))
    pred_dt <- predict(mlmod, newdt2, type = "prob")
    pred_dt <- pred_dt[,which(colnames(pred_dt) == posclass)]
    if(!missing(prepro_obj)){
        newdt2 <- predict(prepro_obj, newdt2)
    }
    newdt_labs <- extractPhenoData(newdt)[[groupvar]]
    dt_class <- unique(newdt_labs)
    idx <- which(dt_class == posclass)
    idx2 <- which(dt_class != posclass)
    pred <- prediction(pred_dt, labels = newdt_labs,
                        label.ordering = dt_class[c(idx2, idx)])
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    perf2 <- performance(pred, measure = "auc")
    auc <- unlist(perf2@y.values)
    plot(perf, main = paste("ROC curve \nAUC = ", round(auc,2)))
}
