# Variance plot
#'@export
setGeneric("featureVarPlot", function(features, varfun, interactive = TRUE){
    standardGeneric("featureVarPlot")
})
#'@export
setMethod("featureVarPlot", definition =  function(features, varfun, interactive){
    feat_dt <- extractData(features)
    feat_sd <- sort(apply(feat_dt, 1, varfun))
    feat_index <- 1:nrow(features)
    p <- varPlot(feat_index = feat_index, feat_sd = feat_sd,
                    interactive = interactive)
    return(p)
})
#'@export
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
metabFeatureFilter <- function(features, groupvar, blankfilt = FALSE,
                            blankFoldChange = 2, blankname = "blank",
                            samplename, cvqcfilt = FALSE,
                            cvqc_thr = 30, qcname = "QC", nafilter = FALSE,
                            naratioThr, naratioMethod, varfilter = FALSE,
                            varfun, varthr, varquant, intensitythr,
                            ism0 = FALSE, hasan = FALSE, sampfilter = FALSE,
                            maxmv, filtername, ...){
    dt_groups <- extractPhenoData(features)[[groupvar]]
    if(blankfilt){
        features <- filter_peaks_by_blank(features, blankFoldChange, dt_groups,
                                        blankname, qcname, remove_peaks = TRUE,
                                        remove_samples = TRUE)
        dt_groups <- extractPhenoData(features)[[groupvar]]
    }
    if(!missing(intensitythr)){
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
    if(sampfilter){
        features <- sampleFilter(features, groupvar = dt_groups,
                                groupname = filtername, maxmv = maxmv)
    }
    return(features)
}
#' #'@export
#' blankFilter <- function(features, groupvar, samplename, blankname){
#'     feat_dt <- extractData(features)
#'     blank_int <- apply(feat_dt, 1, function(x){
#'             s <- mean(x[colData(feat_dt)[[groupvar]] == samplename], na.rm = TRUE)
#'             b <- mean(x[colData(feat_dt)[[groupvar]] == blankname], na.rm = TRUE)
#'             if(is.nan(b)){b <- 0}
#'             s/b > blank.int.ratio.thr
#'         })
#'     return(features[blank_int,])
#' }
#' #'@export
#' cvFilterQC <- function(features, groupvar, samplename, qcname, qctimes = 2){
#'     feat_dt <- extractData(features)
#'     cv_samp_qc <- apply(feat_dt, 1, function(x){
#'             samp_indx <- colData(feat_dt)[[groupvar]] == samplename
#'             qc_indx <- colData(feat_dt)[[groupvar]] == qcname
#'             s <- sd(x[samp_indx], na.rm = TRUE)/mean(x[samp_indx], na.rm = TRUE)
#'             qc <- sd(x[qc_indx], na.rm = TRUE)/mean(x[qc_indx], na.rm = TRUE)
#'             s > (qctimes*qc)
#'         })
#'     return(features[cv_samp_qc,])
#' }
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
model.mat <- function(features, phenovar){
    modmat <- model.matrix(~ 0 + pData(features)[[phenovar]])
    colnames(modmat) <- sub(pattern = "pData(features)[[phenovar]]",
                            replacement = "", x = colnames(modmat),
                            fixed = TRUE)
    return(modmat)
}
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
    annotation <- annotateTable(features, anotpackage)
    rownames(features) <- annotation$SYMBOL
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
    if(!missing(maxmv)){
        features <- filter_samples_by_mv(features, max_perc_mv = maxmv, classes = groupvar, remove_samples = TRUE)
    }
    if(!missing(groupname)){
        features <- features[,groupvar != groupname]
    }
    return(features)
}

#'@export
# Machine learning feature selection
featureSign <- function(features, groupvar, boot){
    feat_dt <- t(extractData(features))
    groupvar <- extractPhenoData(features)[[groupvar]]
    mod <- biosigner::biosign(x = feat_dt, y = groupvar)
    return(mod)
}

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
