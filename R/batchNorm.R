featurebatchQc <- function(features, groups, injectionorder, qcname,
                            logscale = NA, interactive = TRUE){
    ydt <- apply(assay(features), 2, function(x){
        sum(x, na.rm = TRUE)
    })
    dt <- list(colData(features)[[injectionorder]], ydt)
    p <- ggStandardPlot(dt = dt, groups = colData(features)[[groups]],
                    plottype = "scatter",
                    ptitle = "Feature total sum for each sample",
                    xlab = "Injection Order", ylab = "Total sum of features",
                    smoothing = TRUE, ytrans = logscale)
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}

featureBatchPVCA <- function(features, phenovars, threshold){
    covariates <- pData(rawData)[phenovars]
    pvcadt <- pvcaBatchAssess(features, phenovars, threshold)
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

batchNormalization <- function(features, norm, injectionorder, batchnum, groups,
                                qcname, covariate, covariate2){
    if(norm == "qcnorm"){
            return(QCRSC(df = features, order = injectionorder,
                batch = batchnum, classes = groups, qc_label = qcname))
    } else if(norm == "covnorm"){
        return(removeBatchEffect(log10(assay(features)), batch = covariate,
                            batch2 =  covariate2))
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
