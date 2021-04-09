

normalizeData <- function(train_dt, test_dt, method){


}

setGeneric("mlFit", function(dt, ...){
    standardGeneric("mlFit")
})
setMethod("mlFit", signature = "ExpressionSet", function(dt, groupvar, cvmethod = "loo"){
    dt <- t(exprs(dt))
    grp <- pData(dt)[[groupvar]]

    if(nrow(dt) <= 10){
        cvmethod <- "loo"
    }
    if(cvmethod == "loo"){
        dt_indx <- lapply(1:nrow(dt), function(x){
            train_dt <- dt[-x,]
            test_dt <- dt[x,]
            train_grp <- grp[-x]
            test_dt <- grp[x]
            pl_mod <- plsda(train_dt, train_grp)
        })
    }

    pldsa_mod <- plsda(dt, grp)

})
setMethod("mlFit", signature = "SummarizedExperiment", function(dt){

})

mlFitfun <- function(a){


}
