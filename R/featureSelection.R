# Variance plot
#'@export
setGeneric("featureVarPlot", function(features, ...){
    standardGeneric("featureVarPlot")
})
#'@export
setMethod("featureVarPlot", "ExpressionSet",
            function(features, varfun, interactive){
    feat_sd <- sort(apply(exprs(features), 1, varfun))
    feat_index <- 1:nrow(features)
    p <- varPlot(feat_index = feat_index, feat_sd = feat_sd,
                    interactive = interactive)
    return(p)
})
#'@export
setMethod("featureVarPlot", "SummarizedExperiment",
            function(features, varfun, interactive){
    feat_sd <- sort(apply(assay(features), 1, varfun))
    feat_index <- 1:nrow(features)
    p <- varPlot(feat_index = feat_index, feat_sd = feat_sd,
                    interactive = interactive)
    return(p)
})
varPlot <- function(feat_index, feat_sd, interactive = TRUE){
    p <- ggStandardPlot(dt = list(feat_index, feat_sd), plottype = "scatter",
                    ptitle = "Distribution of feature variability",
                    xlab = "Feature Index", ylab = "Standard deviation",
                    angle = 0, groups = NULL)
    p <- p + geom_vline(xintercept = length(feat_index)*0.90) +
        geom_vline(xintercept = length(feat_index)*0.95)
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}

#'@export
geneFeatureFilter <- function(features, entrez, varfilt, varcutoff, var.func){
    return(nsFilter(features, require.entrez = entrez, remove.dupEntrez = TRUE,
                var.filter = varfilt, var.func = var.func, var.cutoff = varcutoff,
                filterByQuantile = TRUE))
}
#'@export
metabFeatureFilter <- function(features, intensity_thr, ism0, hasan){
    feature_info <- data.frame()
    if(!missing(intensity_thr)){
        feature_info$intensity <- apply(assay(features), 1,
                                        function(x){any(x >= thint)})
    }
    if(ism0){
        feature_info$ism0 <- grepl("M0", rowData(features)$isotope)
    }
    if(hasan){
        feature_info$hasan <- !is.na(rpwData(features)$annotation)
    }
    filt_indx <- apply(feature_info, 1, function(x){all(x)})
    return(features[filt_indx,])
}


groupFeatureFilter <- function(features, phenovar){

}

#'@export
# Machine learning feature selection
featureSign <- function(features, phenoVar, boot){
    mod <- biosigner::biosign(features, phenoVar)
    return(mod)
}
#'@export
featureSelection <- function(biosigndata, model = 1, scoremin = "A", minmod){
    lev <- c("E", "B", "A", "S")
    scores <- as.data.frame(biosigndata@tierMC) %>%
        mutate(plsda = factor(plsda, levels = lev, ordered = TRUE)) %>%
        mutate(randomforest = factor(randomforest, levels = lev, ordered = TRUE)) %>%
        mutate(svm = factor(svm, levels = lev, ordered = TRUE)) %>%
        mutate(highscore = plsda >= scoremin | randomforest >= scoremin | svm >= scoremin)
    return(biosigndata@eset[scores$highscore,])
}

