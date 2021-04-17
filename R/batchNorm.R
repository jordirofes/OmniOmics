#'@export
featurebatchQc <- function(features, groupvar, injectionvar, qcname,
                            logscale = FALSE, interactive = TRUE){
    ydt <- apply(assay(features), 2, function(x){
        sum(x, na.rm = TRUE)
    })
    if(!missing(injectionvar)){
        xdt <- colData(features)[[injectionorder]]
    } else{
        xdt <- seq_len(ncol(features))
    }
    if(logscale){
        logscale <- "log10"
    } else{
        logscale <- NA
    }
    dt <- list(xdt, ydt)
    p <- ggStandardPlot(dt = dt, groups = colData(features)[[groupvar]],
                    plottype = "scatter",
                    ptitle = "Total sum of feature intensity for each sample",
                    xlab = "Injection Order", ylab = "Total sum of features",
                    ytrans = logscale)
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@export
featureBatchPVCA <- function(features, phenovars, threshold){
    covariates <- colnames(pData(features))[phenovars]
    pvcadt <- pvcaBatchAssess(features, covariates, threshold)
    p <- ggStandardPlot(dt = pvcadt$dat, plabs = factor(pvcadt$label,
                                                        levels =  pvcadt$label,
                                                        ordered = TRUE),
                    plottype = "column",
                    ptitle = "Princpal variation component analysis (PVCA)",
                    angle = 0, xlab = "Variables", ylab = "Variance") +
        theme(axis.text.x = element_text(hjust = 0.5)) +
        geom_text(aes(x = 1:length(pvcadt$dat), y = pvcadt$dat + 0.01,
                        label = as.character(round(pvcadt$dat, 2))))
    return(p)
}
#'@export
batchNormalization <- function(features, method, injectionorder, batchnum, groups,
                                qcname, covariate, covariate2 = NULL){
    if(method == "qcnorm"){
        if(length(injectionorder) == 1){
            injectionorder <- colData(features)[[injectionorder]]
        }
        batchnum <- colData(features)[[batchnum]]
        groups <- colData(features)[[groups]]
            return(QCRSC(df = features, order = injectionorder,
                batch = batchnum, classes = groups, qc_label = qcname, spar = 0,
                minQC = 4))
    } else if(method == "covnorm"){
        feature_dt <- extractData(features)
        covariate <- extractPhenoData(features)[[covariate]]
        batch_dt <- removeBatchEffect(feature_dt, batch = covariate,
                            batch2 = covariate2)
        exprs(features) <- batch_dt
        return(features)
    } else{
        stop("Enter a valid normalization method")
    }
}
# setGeneric("batchNorma", function(object, method){
#     standardGeneric("batchNorma")
# })
# setMethod("batchNorma", "SummarizedExperiment",
#             function(features, samporder, batchnum, groups){
#     return(normaliseIllumina(object, method = method))
# })
# setMethod("batchNorma", "GeneFeatureSet", function(object){
#     return(rma(object))
# })
