featureVarPlot <- function(features, interactive = TRUE){
    feat_sd <- sort(apply(exprs(features), 1, sd))
    feat_index <- 1:nrow(features)
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

geneFeatureFilter <- function(features, entrez, varfilt, varcutoff){
    return(nsFilter(features, require.entrez = entrez, remove.dupEntrez = TRUE,
                var.filter = varfilt, var.func = IQR, var.cutoff = varcutoff,
                filterByQuantile = TRUE))
}


metabFeatureFilter <- function(features, thint, snr, ism0, hasan){
    rowData(features) %>% mutate(idx_int = assay(features))

    data.frame(int = apply(assay(features), 1, function(x){any(x >= thint)}, snr))
    return(features[filt_indx,])
}

# Machine learning feature selection
featureSign <- function(features, phenoVar, boot){
    mod <- biosigner::biosign(features, phenoVar)
    return(mod)
}

featureSelection <- function(biosigndata, model = 1, scoremin = "A", minmod){
    lev <- c("E", "B", "A", "S")
    scores <- as.data.frame(biosigndata@tierMC) %>%
        mutate(plsda = factor(plsda, levels = lev, ordered = TRUE)) %>%
        mutate(randomforest = factor(randomforest, levels = lev, ordered = TRUE)) %>%
        mutate(svm = factor(svm, levels = lev, ordered = TRUE)) %>%
        mutate(highscore = plsda >= scoremin | randomforest >= scoremin | svm >= scoremin)
    return(biosigndata@eset[scores$highscore,])
}

