#'@title Feature sum intensity plot
#'@author Jordi Rofes Herrera
#'@description Creates a scatter plot with the sum of intensities of all features
#' for each sample ordered by injection order to visualize batch effect
#'@param features A SummarizedExperiment with features
#'@param groupvar A numeric indicating the column to use as grouping factor or a
#'  character indicating it's name
#'@param injectionvar An optional string or numeric indicating the variable
#'in the phenodata used to indicate the injection order, no value will default
#'to the file order
#'@param interactive A boolean indicating if the plot will be converted to an
#'  interactive `ggplotly()`
#'@param logscale log10 intensity transformation
#'@return Returns a ggplot object or a ggplotly
#'@export
featurebatchQc <- function(features, groupvar, injectionvar,
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
#'@title PVCA batch effect plot
#'@author Jordi Rofes Herrera
#'@description Creates a column plot with the variance of each covariate in the
#'principal variation component analysis
#'@param features An ExpressionSet object
#'@param phenovar A numeric or string vector indicating the columns of the pehnodata to include
#'in the pvca analysis
#'@param threshold A precentile value for the minimum variance the selected
#'principal components need to explain
#'@return Returns a ggplot object or a ggplotly
#'@export
featureBatchPVCA <- function(features, phenovars, threshold){
    covariates <- colnames(extractPhenoData(features))[phenovars]
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
#'@title Batch Normalization
#'@author Jordi Rofes Herrera
#'@description Makes the batch normalization for injection order (Metabolomics)
#'using the pmp package and covariates (Transcriptomics) using the pvca package.
#'@param features An ExpressionSet or SummarizedExperiment object
#'@param method A string for the normalization method: qcnorm for
#'injection order and batch normalization and covnorm for covariate normalization
#'@param injectionorder A string or numerical indicating the variable in the phenodata
#'with the injection order or a vector with the injection order.
#'@param batchnum A string or numerical indicating the variable in the phenodata
#'with the batch number
#'@param groups A string or numerical indicating the variable in the phenodata
#'with the grouping variable (QC/SAMPLE/...)
#'@param qcname A string with the name of the QC samples in the groups variable
#'@param covariate A string or numerical indicating the variable in the phenodata
#'with the first covariate to normalize
#'@param covariate2 An optional string or numerical indicating the variable in the phenodata
#'with the second covariate to normalize
#'@return Returns a covariate/batch corrected ExpressionSet or SummarizedExperiment
#'@export
batchNormalization <- function(features, method, injectionorder, batchnum, groups,
                                qcname, covariate, covariate2 = NULL){
    
    if(method == "qcnorm"){
        if(length(injectionorder) == 1){
            injectionorder <- colData(features)[[injectionorder]]
        }
        if(length(batchnum) == 1){
            batchnum <- colData(features)[[batchnum]]
        }
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
