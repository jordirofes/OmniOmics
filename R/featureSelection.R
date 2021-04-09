# Variance plot
#'@export
setGeneric("featureVarPlot", function(features, varfun, interactive = TRUE){
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
# Gene feature filtering
#'@export
geneFeatureFilter <- function(features, entrez, rem.dupEntrez, varfilt, varcutoff, var.func){
    return(nsFilter(features, require.entrez = entrez, remove.dupEntrez = rem.dupEntrez,
                var.filter = varfilt, var.func = var.func, var.cutoff = varcutoff,
                filterByQuantile = TRUE))
}
# Metabolomic basic feature filtering
#'@export
metabFeatureFilter <- function(features, phenovar, blank.int.ratio,
                            blank.int.ratio.thr = 2, blancname = "blank",
                            samplename = "sample", varfilter, varfun, varthr,
                            varquant, intensity_thr, ism0, hasan){
    feat_dt <- assay(features)
    feature_info <- data.frame()
    if(!missing(blank.int.ratio)){
        feature_info$blanc_int <- apply(feat_dt, 1, function(x){
            b <- mean(x[colData(features)[[phenovar]] == samplename])
            s <- mean(x[colData(features)[[phenovar]] == blankname])
            s/b > blank.int.ratio.thr
        })
    }
    if(!missing(intensity_thr)){
        feature_info$intensity <- apply(feat_dt, 1,
                                        function(x){any(x >= thint)})
    }
    if(varfilter){
        feature_info$varfilt <- apply(feat_dt[,colData(features)[[phenovar]] == samplename], 1, varfun)
        if(varquant){
            thr <- sort(feature_info$varfilt)[round(length(feature_info$varfilt) * varthr)]
        } else{
            thr <- varthr
        }
        feature_info$varfilt <- feature_info$varfilt >= thr
    }
    if(ism0){
        feature_info$ism0 <- grepl("M0", rowData(features)$isotope)
    }
    if(hasan){
        feature_info$hasan <- !is.na(rowData(features)$annotation)
    }
    filt_indx <- apply(feature_info, 1, function(x){all(x)})
    return(features[filt_indx,])
}

# Limma gene filtering
model.mat <- function(features, phenovar){
    modmat <- model.matrix(~ 0 + pData(features)[[phenovar]])
    colnames(modmat) <- sub(pattern = "pData(features)[[phenovar]]",
                            replacement = "", x = colnames(modmat),
                            fixed = TRUE)
    return(modmat)
}

groupFeatureComp <- function(features, modelMatrix, contrastMat, adjpval = "fdr"){
    fitmod <- lmFit(features, modelMatrix)
    fitdt <- contrasts.fit(fitmod, contrastMat)
    fitdt <- eBayes(fitdt)
    # fittables <- lapply(seq_along(rownames(fitdt$cov.coefficients)), function(x){
    #     a <- topTable(fitdt, coef = rownames(fitdt$cov.coefficients)[x], number = nrow(fitdt), adjust.method = adjpval)
    #     l <- rownames(fitdt$cov.coefficients)[x]
    # })
    fittables <- vector(length = ncol(fitdt$cov.coefficients))
    for(i in seq_along(rownames(fitdt$cov.coefficients))){
        a <- topTable(fitdt, coef = rownames(fitdt$cov.coefficients)[i], number = nrow(fitdt), adjust.method = adjpval)
        a$compname <- rep(rownames(fitdt$cov.coefficients)[i], nrow(a))
        fittables[i] <- list(a)
    }
    return(fittables)
}
groupFeatureVolcano <- function(comptable, interactive = TRUE){
    compname <- comptable$compname[1]
    dt <- comptable[,c("logFC", "adj.P.Val")]
    dt$adj.P.Val <- -log10(dt$adj.P.Val)
    dt <- as.list.data.frame(dt)
    p <- ggStandardPlot(dt = dt, groups = NULL, plottype = "scatter",
                        ptitle = paste("Volcano plot for:", compname),
                        xlab = "Fold-Change (Log2)", ylab = "p-values (-log10)",
                        angle = 0) +
        geom_vline(xintercept = -1) + geom_vline(xintercept = 1) +
        geom_hline(yintercept = -log10(0.05))
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}

groupFeatureFilter <- function(features, comptable, pvalthr, logFCthr){
    tokeep <- which(abs(comptable[["logFC"]]) >= logFCthr & comptable[["adj.P.Val"]] <= pvalthr)
    return(features[tokeep,])
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

