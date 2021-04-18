#'@export
setGeneric("groupComp", function(features, groupvar, elimlab = "none", test,
                                adj.method = "fdr", paired, var.equal = FALSE){
    standardGeneric("groupComp")
})
#'@export
setMethod("groupComp", definition = function(features, groupvar,
                                            elimlab = "none", test, adj.method,
                                            paired, var.equal = FALSE){
    group_tables <- featureGroupTables(features, groupvar, elimlab)
    groupvar <- extractPhenoData(features)[[groupvar]]
    comp_tables <- list()
    for(i in seq_along(group_tables)){
        for(j in seq_along(group_tables)){
            if(j < i){
                comp_tables <- c(comp_tables,
                                list(matComp(mat1 = group_tables[[i]],
                                    mat2 = group_tables[[j]], test = test,
                                    adj.method = adj.method, paired = paired,
                                    gr = c(unique(groupvar)[i],
                                            unique(groupvar)[j]),
                                    var.equal = var.equal)))
            }
        }
    }
    return(comp_tables)
})
#'@export
matComp <- function(mat1, mat2, test, adj.method = "fdr", paired = FALSE,
                    var.equal, gr){
    mat1 <- mat1
    mat2 <- mat2
    cm <- data.frame(row.names = colnames(mat1))
    cm$logFC <- sapply(1:ncol(mat1), function(x){
        m1 <- mean(mat1[,x], na.rm = TRUE)
        m2 <- mean(mat2[,x], na.rm = TRUE)
        log2(m2/m1)
    })
    if(test == "t-test"){
        tt <- lapply(1:ncol(mat1), function(x){
            t <- t.test(mat1[,x], mat2[,x], paired = paired,
                        var.equal = var.equal)
            return(c(t$statistic, t$p.value))
        })
    } else if(test == "wilcox"){
        tt <- lapply(1:ncol(mat1), function(x){
            t <- wilcox.test(mat1[,x], mat2[,x])
            return(c(t$statistic, t$p.value))
        })
    }
    tt <- do.call(rbind, tt)
    cm$t.stad <- tt[,1]
    cm$P.Value <- tt[,2]
    cm$adj.P.Val <- p.adjust(cm$P.Value, method = adj.method)
    cm$compname <- rep(paste(gr[1], "VS", gr[2]), ncol(mat1))
    return(cm)
}
#'@export
featureGroupTables <- function(features, groupvar, elimlab){
    groupvar <- extractPhenoData(features)[[groupvar]]
    feature_dt <- extractData(features)
    group_tables <- lapply(unique(groupvar)[unique(groupvar) != elimlab],
                            function(x){
        group_table <- t(feature_dt[,groupvar == x])
        return(group_table)
    })
    return(group_tables)
}

#'@export
setGeneric("multipca", function(features, prefuns, method, groupvar, qcname, scale = TRUE, ...){
    standardGeneric("multipca")
})
#'@export
setMethod("multipca", definition = function(features, prefuns, method, groupvar,
                                            qcname, scale = TRUE){
    if(!missing(prefuns)){
        groupvar <- extractPhenoData(features)
        features <- prePro(features, prefuns, method, groupvar, qcname)
    }
    dt <- extractData(features)
    pcdt <- prcomp(t(dt), scale = scale)
    return(pcdt)
})
#'@export
pcaPlot <- function(pcdt, object, pc = c(1,2), groupvar, interactive = TRUE,
                    ellipse = TRUE, ellipse.type = "t", level = 0.95,
                    bygroup = TRUE){
    if(length(pc) != 2){
        pc <- c(1,2)
    }
    plabs <- colnames(object)
    groups <- extractPhenoData(object)[[groupvar]]
    pcvar <- round(pcdt$sdev^2/sum(pcdt$sdev^2)*100, 2)
    pcvar <- pcvar[pc]
    pcdt <- pcdt$x[,pc]
    pcdt <- as.list.data.frame(as.data.frame(pcdt))
    p <- ggStandardPlot(dt = pcdt, plabs = plabs, groups = groups,
                        plottype = "scatter", ptitle = "PCA plot",
                        xlab = paste0("PC ", pc[1], " (", pcvar[pc[1]], "%)"),
                        ylab = paste0("PC ", pc[2], " (", pcvar[pc[2]], "%)"),
                        origin_lines = TRUE)
    if(ellipse){
        p <- p + stat_ellipse(aes(x = pcdt[[1]], y = pcdt[[2]]),
                            colour = "grey7", type = ellipse.type, linetype = 1,
                            level = level)
    }
    if(bygroup){
        if(any(table(plabs) < 5)){
            message("Too few points to calculate group ellipse")
        } else{
            p <- p + stat_ellipse(aes(x = pcdt[[1]], y = pcdt[[2]],
                                color = plabs), type = ellipse.type,
                                linetype = 1, level = level)
        }
    }
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@export
loadingsPlot <- function(pcdt, pc = 1, interactive = TRUE){
    dt <- data.frame(Features = rownames(pcdt$rotation),
                    Rotation = pcdt$rotation[,pc])
    p <- ggdotchart(dt, x = "Features", y = "Rotation", sorting = "ascending",
            add = "segments", ggtheme = theme_pubr()) +
        ggtitle(paste("PCA loading plot from pc", pc)) +
        theme(axis.text.x = element_text(size = 9),
            plot.title = element_text(hjust = 0.5))

    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@export
heatmapPlot <- function(features, groupvar, annotated){
    feature_dt <- extractData(features)
    group_dt <- extractPhenoData(features)[[groupvar]]
    p <- heatmap.2(feature_dt, Rowv = FALSE, Colv = TRUE,
            main = "Feature heatmap",
            scale = "row", sepcolor = "white",
            sepwidth = c(0.05,0.05), cexRow = 0.5, cexCol = 0.9, key = TRUE,
            keysize = 1.5, density.info = "histogram",
            tracecol = NULL, srtCol = 30, labCol = group_dt)
    return(p)
}


#'@export
setGeneric("compPlot", function(features, subset, groupvar, ...){
    standardGeneric("compPlot")
})
#'@export
setMethod("compPlot", definition = function(features, subset, groupvar, pal, ...){
    if(!missing(subset)){
        features <- features[subset,]
    }
    inte <- extractData(features)
    gr <- extractPhenoData(features)[[groupvar]]
    feat_dt <- lapply(rownames(inte), function(x){
        data.frame(intensity = inte[x,], group = gr)
    })
    if(missing(pal)){
        pal <- brewer.pal(n = length(unique(gr)), name = "Set1")
    }
    xlabs <- rownames(features)
    p <- lapply(seq_along(feat_dt), function(x){
        featureCompPlot(feat_dt[[x]], pal, xlabs[x], ...)
    })
    return(p)
})
#'@export
featureCompPlot <- function(dt, pal, xlab, ...){
    p <- ggpubr::ggboxplot(dt, x = "group", y = "intensity", shape = "group",
                            add = "jitter", color = "group", palette = pal,
                            xlab = xlab, ylab = "Intensity")
    return(p)
}
#'@export
prePro <- function(features, prefuns, method, groupvar, qcname, blankname, ...){
    groupvar <- extractPhenoData(features)[[groupvar]]
    for(i in 1:length(prefuns)){
        if(prefuns[i] == "blankfill"){
            features <- blankFill(features, groupvar, blankname)
        }
        if(prefuns[i] == "mvImp"){
            features <- mv_imputation(features, method = method)
        }
        if(prefuns[i] == "pqn"){
            features <- pqn_normalisation(features, groupvar, qcname)
        }
        if(prefuns[i] == "glog"){
            features <- glog_transformation(features, groupvar, qcname)
        }
        if(prefuns[i] == "sum"){
            features <- normalise_to_sum(features)
        }
    }
    return(features)
}
#'@export
blankFill <- function(features, groupvar, blankname){
    assay(features)[,groupvar == blankname][is.na(assay(features)[,groupvar == blankname])] <- 0
    return(features)
}
