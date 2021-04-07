setGeneric("groupComp", function(features, groupvar, elimlab, test, adj.method,
                                paired, var.equal = FALSE){
    standardGeneric("groupComp")
})
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

matComp <- function(mat1, mat2, test, adj.method = "fdr", paired = FALSE, var.equal,
                    gr){
    cm <- data.frame(row.names = colnames(mat1))
    cm$fold.change <- sapply(1:ncol(mat1), function(x){
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
    cm$adj.p.val <- cm %>% mutate(p.adjust(p.value, method = adj.method))
    cm$comp.group <- rep(paste(gr[1], "VS", gr[2]), nrow(mat1))
    return(cm)
}

setGeneric("multipca", function(object, scale = TRUE){
    standardGeneric("multipca")
})
setMethod("multipca", c("GeneFeatureSet","ExpressionSet"),
            function(object, scale = TRUE){

    dt <- exprs(object)
    pcdt <- prcomp(t(dt), scale = scale)
    return(pcdt)
})
setMethod("multipca", c("SummarizedExperiment"),
            function(object, scale = TRUE){
    dt <- assay(object)
    dt[which(is.na(dt))] <- 0
    pcdt <- prcomp(t(dt), scale = scale)
    return(pcdt)
})
pcaPlot <- function(pcdt, dt, pc = c(1,2), groupvar, interactive = TRUE){
    if(length(pc) != 2){
        pc <- c(1,2)
    }
    plabs <- colnames(dt)
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
    if(interactive){
        return(ggplotly(p))
    }
}


heatmapPlot <- function(features){
    p <- heatmap.2(features, Rowv = FALSE, Colv = FALSE,
            main = "Feature heatmap",
            scale = "row", sepcolor = "white",
            sepwidth = c(0.05,0.05), cexRow = 0.5, cexCol = 0.9, key = TRUE,
            keysize = 1.5, density.info = "histogram",
            tracecol = NULL, dendrogram = "none", srtCol = 30)
    return(p)
}
