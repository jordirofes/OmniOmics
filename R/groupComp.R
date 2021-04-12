#'@export
setGeneric("groupComp", function(features, groupvar, elimlab, test, adj.method,
                                paired, var.equal = FALSE){
    standardGeneric("groupComp")
})
#'@export
setMethod("groupComp", "ExpressionSet", function(features, groupvar, elimlab,
                                                test, adj.method, paired,
                                                var.equal = FALSE){
    groupvar <- pData(features)[[groupvar]]
    group_tables <- lapply(unique(groupvar)[unique(groupvar) != elimlab], function(x){
        return(t(exprs(features)[,groupvar == x]))
    })
    comp_tables <- list()
    for(i in seq_along(group_tables)){
        for(j in seq_along(group_tables)){
            if(j < i){
                comp_tables <- c(comp_tables, list(matComp(mat1 = group_tables[[i]],
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
setMethod("groupComp", "SummarizedExperiment", function(features, groupvar, elimlab, test, adj.method){
    groupvar <- colData(features)[[groupvar]]
    group_tables <- lapply(unique(groupvar)[unique(groupvar) != elimlab], function(x){
        return(t(assay(features)[,groupvar == x]))
    })
    comp_tables <- list()
    for(i in seq_along(group_tables)){
        for(j in seq_along(group_tables)){
            if(j < i){
                comp_tables <- c(comp_tables, list(matComp(mat1 = group_tables[[i]],
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
matComp <- function(mat1, mat2, test, adj.method = "fdr", paired = FALSE, var.equal,
                    gr){
    cm <- data.frame(row.names = colnames(mat1))
    cm$logFC <- sapply(1:ncol(mat1), function(x){
        m1 <- mean(mat1[,x])
        m2 <- mean(mat2[,x])
        log2(m2/m1)
    })
    tt <- lapply(1:ncol(mat1), function(x){
        t <- t.test(mat1[,x], mat2[,x], paired = paired, var.equal = var.equal)
        return(c(t$statistic, t$p.value))
    })
    tt <- do.call(rbind, tt)
    cm$t.stad <- tt[,1]
    cm$p.value <- tt[,2]
    cm$adj.P.Val <- p.adjust(cm$p.value, method = adj.method)
    cm$comp.group <- rep(paste(gr[1], "VS", gr[2]), ncol(mat1))
    return(cm)
}
#'@export
setGeneric("multipca", function(object, scale = TRUE){
    standardGeneric("multipca")
})
#'@export
setMethod("multipca", definition = function(object, scale = TRUE){
    dt <- extractData(object)
    pcdt <- prcomp(t(dt), scale = scale)
    dt[which(is.na(dt))] <- 0
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
    groups <- phenoData(object)[[groupvar]]
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
heatmapPlot <- function(features){
    p <- heatmap.2(features, Rowv = FALSE, Colv = FALSE,
            main = "Feature heatmap",
            scale = "row", sepcolor = "white",
            sepwidth = c(0.05,0.05), cexRow = 0.5, cexCol = 0.9, key = TRUE,
            keysize = 1.5, density.info = "histogram",
            tracecol = NULL, dendrogram = "none", srtCol = 30)
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
    gr <- pData(features)[[groupvar]]
    feat_dt <- lapply(rownames(inte), function(x){
        data.frame(intensity = inte[x,], group = gr)
    })
    if(missing(pal)){
        pal <- brewer.pal(n = length(unique(gr)), name = "Set1")
    }

    p <- lapply(feat_dt, function(x){
        featureCompPlot(x, pal, ...)
    })
    return(p)
})
#'@export
featureCompPlot <- function(dt, pal, ...){
    p <- ggpubr::ggboxplot(dt, x = "group", y = "intensity", shape = "group", add = "jitter", color = "group", palette = pal)
    return(p)
}
