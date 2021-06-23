# Variance plot
#'@title Variance Plot
#'@author Jordi Rofes Herrera
#'@description Creates a scatter plot with the ordered feature values after
#'evaluating the given function for each feature usually used to visualize the
#'distribution of IQR, standard deviation or coefficient of variation (CV) for each feature.
#'@param features A SummarizedExperiment or ExpressionSet
#'@param varfun A function to calculate for each feature, usually IQR / SD or function(x){sd(x)/mean(x)}
#'@param interactive A boolean indicating if the plot will be converted to an
#'  interactive `ggplotly()`
#'@return Returns a ggplot object or a ggplotly
#'@export
setGeneric("featureVarPlot", function(features, varfun, interactive = TRUE){
    standardGeneric("featureVarPlot")
})
#'@export
setMethod("featureVarPlot", definition =  function(features, varfun, interactive){
    feat_dt <- extractData(features)
    feat_sd <- sort(apply(feat_dt, 1, varfun, na.rm = TRUE))
    feat_index <- 1:nrow(features)
    p <- varPlot(feat_index = feat_index, feat_sd = feat_sd,
                    interactive = interactive)
    return(p)
})
#'@export
varPlot <- function(feat_index, feat_sd, interactive = TRUE){
    p <- ggStandardPlot(dt = list(feat_index, feat_sd), plottype = "scatter",
                    ptitle = "Distribution of feature variability",
                    xlab = "Feature Index", ylab = "Variation",
                    angle = 0, groups = NULL)
    p <- p + geom_vline(xintercept = length(feat_index)*0.90) +
        geom_vline(xintercept = length(feat_index)*0.95)
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
# Gene feature filtering
#'@title Gene Feature Filter
#'@author Jordi Rofes Herrera
#'@description Wrapper for the nsFilter() function. Use ?nsFilter for param information
#'@return Returns filtered gene features
#'@export
geneFeatureFilter <- function(features, entrez, rem.dupEntrez, varfilt, varcutoff, var.func){
    return(nsFilter(features, require.entrez = entrez, remove.dupEntrez = rem.dupEntrez,
                var.filter = varfilt, var.func = var.func, var.cutoff = varcutoff,
                filterByQuantile = TRUE)$eset)
}
# Metabolomic basic feature filtering
#'@title Metablomics feature filter
#'@author Jordi Rofes Herrera
#'@description Filters metabolomic features with a variety of criteria using different functions
#'@param features A SummarizedExperiment
#'@param groupvar A string or numeric indicating the variable in the phenodata
#'with the group variable (sample/qc/blank)
#'@param blankfilt Boolean indicating blank filtering
#'@param blankFoldChange A numeric indicating the fold change to filter features with high intensity in blank
#'@param blankname A string indicating the blank name in the groupvar
#'@param samplename A string indicating the sample name in the groupvar
#'@param cvqcfilt A boolean indicating the filter of features with high RDS in QC samples
#'@param cvqc_thr A numeric indicating the maximum percentage of RDS in QC samples
#'@param qcname A string indicating the qcname in the groupvar variable
#'@param nafilter A boolean indicating the filter of features with high percentage of missing values
#'@param naratioThr A numeric indicating the fraction of NA values allowed in the features
#'@param naratioMethod A string indicating the method for calculating the NA ratio: can be "QC" for within QC samples,
#'"within" for within each sample class in groupvar and "across" for across all samples.
#'@param varfilter A boolean indicating filtering depending on variability criteria
#'@param varfun A function to calculate for each feature, usually IQR / SD or function(x){sd(x)/mean(x)}
#'@param varthr The minimum threshold of variability to keep the features
#'@param varquant A boolean indicating if the threshold is given as a precentile or an absolute value (TRUE/FALSE)
#'@param intfilter A boolean indicating to filter by sample intensity
#'@param intensitythr A value indicating the minimum intensity of atleast one sample in each feature
#'@param ism0 A boolean indicating if filtering of features with non M0 annotations
#'(Requires annotated SummarizedExperiment as done with the pre-processing functions)
#'@param hasan A boolean indicating to filter features with no annotation
#'(Requires annotated SummarizedExperiment as done with the pre-processing functions)
#'@param samplfilter A boolean indicating to filter samples with high missing values
#'@param maxmv A numeric indicating the fraction of maximum missing values in each sample
#'@param filtername An optional name indicating the name of a variable in groupvar to eliminate all it's samples.
#'Useful to eliminate QC samples after using them for all pre-processing.
#'@param prepro A boolean indicating to do a pre-processing of all features
#'@param preprofuns A vector of strings indicating all the functions applied to the features. Available options are:
#'"pqn" point quotien normalization, "sum" normalization to the sum and "mvImp" multivariate imputation.
#'@param mvimpmethod A string indicating the multivariate imputation method: "knn" or "rf".
#'@return Returns a SummarizedExperiment with filtered features
#'@export
metabFeatureFilter <- function(features, groupvar, blankfilt = FALSE,
                            blankFoldChange = 2, blankname = "blank",
                            samplename, cvqcfilt = FALSE,
                            cvqc_thr = 30, qcname = "QC", nafilter = FALSE,
                            naratioThr, naratioMethod, varfilter = FALSE,
                            varfun, varthr, varquant, intfilter = FALSE,
                            intensitythr, ism0 = FALSE, hasan = FALSE,
                            sampfilter = FALSE, maxmv, filtername, prepro,
                            preprofuns, mvimpmethod){
    dt_groups <- extractPhenoData(features)[[groupvar]]

    if(blankfilt){
        features <- filter_peaks_by_blank(features, blankFoldChange, dt_groups,
                                        blankname, remove_peaks = TRUE,
                                        remove_samples = TRUE)
        dt_groups <- extractPhenoData(features)[[groupvar]]
    }
    if(intfilter){
        features <- intFilter(features, intensitythr)
    }
    if(cvqcfilt){
        features <- filter_peaks_by_rsd(features, classes = dt_groups,
                                        max_rsd = cvqc_thr, qc_label = qcname,
                                        remove_peaks = TRUE)
    }
    if(nafilter){
        features <- filter_peaks_by_fraction(features, min_frac = naratioThr,
                                            dt_groups, method = naratioMethod,
                                            qc_label = qcname,
                                            remove_peaks = TRUE)
    }
    if(varfilter){
        features <- varfunFilter(features, varfun, varquant, varthr, dt_groups,
                                samplename)
    }
    if(ism0){
        features <- ismoFilter(features)
    }
    if(hasan){
        features <- anotFilter(features)
    }
    if(prepro){
        features <- prePro(features, prefuns = preprofuns, method = mvimpmethod,
                        dt_groups, qcname, blankname)
    }
    if(sampfilter){
        features <- sampleFilter(features, groupvar = dt_groups,
                                groupname = filtername, maxmv = maxmv)
    }

    return(features)
}
#
# blankFilter <- function(features, groupvar, samplename, blankname){
#     feat_dt <- extractData(features)
#     blank_int <- apply(feat_dt, 1, function(x){
#             s <- mean(x[colData(feat_dt)[[groupvar]] == samplename], na.rm = TRUE)
#             b <- mean(x[colData(feat_dt)[[groupvar]] == blankname], na.rm = TRUE)
#             if(is.nan(b)){b <- 0}
#             s/b > blank.int.ratio.thr
#         })
#     return(features[blank_int,])
# }
#
# cvFilterQC <- function(features, groupvar, samplename, qcname, qctimes = 2){
#     feat_dt <- extractData(features)
#     cv_samp_qc <- apply(feat_dt, 1, function(x){
#             samp_indx <- colData(feat_dt)[[groupvar]] == samplename
#             qc_indx <- colData(feat_dt)[[groupvar]] == qcname
#             s <- sd(x[samp_indx], na.rm = TRUE)/mean(x[samp_indx], na.rm = TRUE)
#             qc <- sd(x[qc_indx], na.rm = TRUE)/mean(x[qc_indx], na.rm = TRUE)
#             s > (qctimes*qc)
#         })
#     return(features[cv_samp_qc,])
# }
#'@export
varfunFilter <- function(features, varfun, varquant, varthr, groupvar, samplename){
        feat_dt <- extractData(features)
        if(!missing(samplename)){
            feat_dt <- feat_dt[,groupvar == samplename]
        }
        varfilt <- apply(feat_dt, 1, varfun)
        if(varquant){
            thr <- sort(varfilt)[round(length(varfilt) * varthr)]
        } else{
            thr <- varthr
        }
        varfilt <- varfilt >= thr
        return(features[varfilt,])
}
#'@export
anotFilter <- function(features){
    feat_dt <- extractData(features)
    hasan <- !is.na(rowData(feat_dt)$annotation)
    return(features[hasan,])
}
#'@export
ismoFilter <- function(features){
    feat_dt <- extractData(features)
    ism0 <- grepl("M0", rowData(feat_dt)$isotope)
    return(features[ism0,])
}
#'@export
intFilter <- function(features, intthr){
    feat_dt <- extractData(features)
    intensity <- apply(feat_dt, 1, function(x){any(x >= intthr)})
    features[intensity,]
}

#'@export
# Limma gene filtering
#'@title Model matrix
#'@author Jordi Rofes Herrera
#'@description Creates a model matrix from an ExpressionSet
#'@param features An ExpressionSet object
#'@param phenovar A string or numeric indicating the group variable in the phenodata
#'to create the model matrix with
#'@return Returns filtered gene features
#'@export
model.mat <- function(features, phenovar){
    modmat <- model.matrix(~ 0 + extractPhenoData(features)[[phenovar]])
    colnames(modmat) <- sub(pattern = "extractPhenoData(features)[[phenovar]]",
                            replacement = "", x = colnames(modmat),
                            fixed = TRUE)
    return(modmat)
}
#'@title Limma group comparison
#'@author Jordi Rofes Herrera
#'@description Calculates comparison tables with the Limma package for each group pair
#'@param features An ExpressionSet object
#'@param modelMatrix A model matrix
#'@param contrastMat A contrast matrix
#'@param adjpval A string indicating the p-value adjust methodology "fdr" or "bonferroni"
#'to create the model matrix with
#'@return A list of comparisson tables for each pair of groups
#'@export
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
#'@title Annotate Data
#'@author Jordi Rofes Herrera
#'@description Annotates an ExpressionSet, SummarizedExperiment or comparison table.
#'@param features An ExpressionSet, SummarizedExperiment
#'@param tableList A list of comparison tables
#'@param anotpackage An annotation package for the corresponding microarray
#'@param metabList A list of metabolites for each feature
#'@return Returns an annotated ExpressionSet, SummarizedExperiment or comparison tables
#'@export
setGeneric("annotateData", function(features, tableList, anotpackage, metabList){
    standardGeneric("annotateData")
})
#'@export
setMethod("annotateData", c("ANY", "list", "character"), function(features, tableList, anotpackage){
    annotated_tables <- lapply(seq_along(tableList), function(x){
        annot <- annotateTable(tableList[[x]], anotpackage)
        rownames(tableList[[x]]) <- annot$SYMBOL
        cbind(tableList[[x]], annot)
    })
    return(annotated_tables)
})
#'@export
setMethod("annotateData", "ExpressionSet", function(features, tableList, anotpackage){
    features <- features[!is.na(rownames(features)),]
    annotation <- annotateTable(features, anotpackage)
    features <- features[!is.na(annotation$SYMBOL),]
    rownames(features) <- annotation$SYMBOL[!is.na(annotation$SYMBOL)]
    return(features)
})
#'@export
setMethod("annotateData", "data.frame", function(features, tableList, anotpackage){
    annotation <- annotateTable(features, anotpackage)
    features <- features[!is.na(annotation$SYMBOL),]
    rownames(features) <- annotation$SYMBOL[!is.na(annotation$SYMBOL)]
    return(features)
})
#'@export
setMethod("annotateData", "SummarizedExperiment", function(features, metabList){
    rownames(features) <- metabList
    return(features)
})
#'@export
annotateTable <- function(feat_dt, anotpackage){
    genes <- rownames(feat_dt)
    packData <- eval(parse(text = anotpackage))
    genes_anotat <- AnnotationDbi::select(packData, genes, c("SYMBOL", "ENTREZID", "GENENAME"))
    return(genes_anotat)

}
#'@title Volcano plot
#'@author Jordi Rofes Herrera
#'@description Creates a volcano plot for a comparison table
#'@param comptable A comparison table
#'@param adj.pvalue Boolean indicating if the adjusted p-values will be used
#'@param interactive Boolean indicating if a plotly will be output instead of a ggplot
#'@param jitterseed A seed for the random jitter to avoid names overlaps
#'@return Returns filtered gene features
#'@export
groupFeatureVolcano <- function(comptable, adj.pvalue = TRUE, interactive = TRUE, jitterseed = 123){
    compname <- comptable$compname[1]

    if(adj.pvalue){
        dt <- comptable[,c("logFC", "adj.P.Val")]
        dt_thr <- dt$adj.P.Val <= 0.05 & abs(dt$logFC) > 1
        groups <- ifelse(dt_thr, "royalblue4", "grey40")
        dt$adj.P.Val <- -log10(dt$adj.P.Val)
    } else{
        dt <- comptable[,c("logFC", "P.Value")]
        dt_thr <- dt$P.Value <= 0.05 & abs(dt$logFC) > 1
        groups <- ifelse(dt_thr, "royalblue4", "grey40")
        dt$P.Value <- -log10(dt$P.Value)
    }
    dt_labs <- rownames(comptable)
    dt <- as.list.data.frame(dt)
    p <- ggStandardPlot(dt = dt, groups = NULL, plottype = "scatter",
                        ptitle = paste("Volcano plot for:", compname),
                        xlab = "Fold-Change (Log2)", ylab = "p-values (-log10)",
                        angle = 0) +
        geom_vline(xintercept = -1, color = "lightgrey") + geom_vline(xintercept = 1, color = "lightgrey") +
        geom_hline(yintercept = -log10(0.05),color = "lightgrey") +
        ggplot2::annotate(geom = "point",x = dt[[1]], y = dt[[2]], colour = groups) +
        geom_text(aes(x = dt[[1]][dt_thr], y = dt[[2]][dt_thr], label = dt_labs[dt_thr]),
                  position = position_jitter(width = 0.1, height = 0.1, seed = jitterseed))
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@export
groupFeatureFilter <- function(features, comptable, pvalthr = 0.05, logFCthr = 1, padjusted){
    if(padjusted){
        tokeep <- which(abs(comptable[["logFC"]]) >= logFCthr & comptable[["adj.P.Val"]] <= pvalthr)
    } else{
        tokeep <- which(abs(comptable[["logFC"]]) >= logFCthr & comptable[["P.Value"]] <= pvalthr)
    }
    return(features[tokeep,])
}
#'@export
sampleFilter <- function(features, groupvar, groupname, maxmv){
    if(!missing(groupname)){
        features <- features[,groupvar != groupname]
        groupvar <- groupvar[groupvar != groupname]
    }

    if(!missing(maxmv)){
        features <- filter_samples_by_mv(features, max_perc_mv = maxmv, classes = groupvar, remove_samples = TRUE)
    }

    return(features)
}

#'@title Feature Sign
#'@author Jordi Rofes Herrera
#'@description Uses the biosign function from the biosigner package to create a biosign object for the given features
#'@param features An ExpressionSet or SummarizedExperiment object
#'@param groupvar A string or numeric indicating the group variable in the phenodata for classification (two factor only).
#'@return Returns a biosign object
#'@export

# Machine learning feature selection
featureSign <- function(features, groupvar, boot){
    feat_dt <- t(extractData(features))
    groupvar <- extractPhenoData(features)[[groupvar]]
    mod <- biosigner::biosign(x = feat_dt, y = groupvar)
    return(mod)
}

#'@title Feature Selection
#'@author Jordi Rofes Herrera
#'@description Filters features depending on it's score in a biosign model
#'@param features An ExpressionSet or SummarizedExperiment object
#'@param biosigndata A biosign object
#'@param model A numeric/string vector indicating the models from which the scores will be used.
#'@param scoremin The minimum score to select the features "A", "S", "B".
#'@return Returns a filtered ExpressionSet or SummarizedExperiment object.
#'@export
featureSelection <- function(features, biosigndata, model = 1, scoremin = "A"){
    lev <- c("E", "B", "A", "S")
    scores <- as.data.frame(biosigndata@tierMC) %>%
        mutate(plsda = factor(plsda, levels = lev, ordered = TRUE)) %>%
        mutate(randomforest = factor(randomforest, levels = lev, ordered = TRUE)) %>%
        mutate(svm = factor(svm, levels = lev, ordered = TRUE))
    highscore <- apply(as.data.frame(scores[,model]) >= scoremin, 1, any)
    return(features[highscore,])
}
#'@export
setGeneric("extractData", function(dt){
    standardGeneric("extractData")
})
#'@export
setMethod("extractData", "ExpressionSet", function(dt){
    return(exprs(dt))
})
#'@export
setMethod("extractData", "GeneFeatureSet", function(dt){
    return(exprs(dt))
})
#'@export
setMethod("extractData", "SummarizedExperiment", function(dt){
    return(assay(dt))
})
#'@export
setGeneric("extractPhenoData", function(dt){
    standardGeneric("extractPhenoData")
})
#'@export
setMethod("extractPhenoData", "ExpressionSet", function(dt){
    return(pData(dt))
})
#'@export
setMethod("extractPhenoData", "GeneFeatureSet", function(dt){
    return(pData(dt))
})
#'@export
setMethod("extractPhenoData", "SummarizedExperiment", function(dt){
    return(colData(dt))
})
#'@export
setMethod("extractPhenoData", "OnDiskMSnExp", function(dt){
    return(phenoData(dt)@data)
})
#'@export
setMethod("extractPhenoData", "XCMSnExp", function(dt){
    return(phenoData(dt)@data)
})
