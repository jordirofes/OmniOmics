# Creates a ggplot2 chromatogram or EIC from an XCMS object
ggChromPlot <- function(object, filenum, mz = NA, ppm = NA, pheno_var = 1,
                        chromtype = "max", logscale = TRUE, interactive = TRUE){
    if(any(is.null(mz))){
        mz <- c(-Inf, Inf)
    } else{
        if(any(is.na(mz))){
            mz <- c(-Inf, Inf)
        } else if(length(mz) != 2 & !is.na(ppm)){
            err <- (ppm/(10^6))*mz
            mz <- c(mz - err, mz + err)
        }
    }
    if(logscale){
        logscale <- "log10"
    }
    # Colors by sample names
    groups <- sub(object@phenoData@data[,pheno_var][filenum], pattern =
                                "\\.mzX?ML",replacement =  "")
    if(!any(is.na(filenum))){
        object <- filterFile(object, filenum)
    }
    # Chromatogram creation
    chrom_dt <- chromatogram(object, aggregationFun = chromtype, mz = mz)
    if(any(is.infinite(mz))){
        mz <- round(mz(chrom_dt[,1]), 4)
    }
    dt <- lapply(seq_along(chrom_dt), function(x){
        cbind(chrom_dt[[x]]@rtime,chrom_dt[[x]]@intensity)})
    p <- ggStandardPlot(dt = dt, groups = groups,
                        plottype = "line", ytrans = logscale,
                        ptitle = paste(mz[1], "-", mz[2], "mz",
                                        collapse = ""), xlab = "rt",
                        ylab = "Intensity", angle = 0)
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}

# Creates a violin/boxplot TIC plot from an XCMS object
ggTicQuality <- function(object, filenum = NA, pheno_var = 2, logscale = TRUE,
                        violin, interactive = TRUE, injection_order = TRUE,
                        pheno_order = 1){
    if(injection_order){
        plot_order <- findOrder(object, pheno_order)
    } else{
        plot_order <- seq_along(dim(raw_dt[[1]])[1])
    }
    # Data preparation
    if(any(is.na(filenum))){
        filenum <- seq_along(object@phenoData@data[,pheno_var])
    }
    if(logscale){
        logscale <- "log10"
    }
    # Creates de data
    bpdt <- split(tic(filterFile(object, filenum)),
                  f = fromFile(filterFile(object, filenum)))
    plabs <- sub(pattern = "\\.[mM][zZ]?[xX][mM][lL]$", replacement = "",
                x =  basename(fileNames(object)))[filenum]
    groups <- object@phenoData@data[,pheno_var][filenum]
    plottype <- ifelse(violin == TRUE, "violin", "boxplot")
    p <- ggStandardPlot(dt = bpdt, plabs = plabs, groups = groups,
                        plottype = plottype, ytrans = logscale,
                        ptitle = "Total Ion Count plot", xlab = "Sample",
                        ylab = "Intensity", plotorder = plot_order)
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}

ggStandardPlot <- function(dt, plabs, groups, plottype, ptitle, ytrans = NA,
                            xtrans = NA, xlab, ylab, angle = 90,
                            plotorder = dt, origin_lines = FALSE,
                            smoothing = FALSE){
    # plabs <- factor(plabs[plotorder], ordered = TRUE)
    # Creates the base plot
    chrom_base <- ggplot() + theme_classic() + xlab(xlab) +
    ylab(ylab) + ggtitle(ptitle) +
    theme(axis.line = element_line(colour = "black", size = 0.5,
                                    linetype = "solid"),
                                    legend.title = element_blank(),
                                    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = angle, vjust = 0.5, hjust = 1))


    # Creates a different type of plot for a given list of sample points
    # depending on the given plottype
    if(plottype == "violin"){
        chrom_plot <- lapply(seq_along(plotorder), function(x){
            geom_violin(aes(x = plabs[x], y = dt[[x]], fill = groups[x]))
        })
    } else if(plottype == "boxplot"){
        chrom_plot <- lapply(seq_along(plotorder), function(x){
            geom_boxplot(aes(x = plabs[plotorder[x]], y = dt[[plotorder[x]]],
                                fill = groups[x]))
        })
    } else if(plottype == "scatter"){
        chrom_plot <- list(geom_point(aes(x = dt[[1]], y = dt[[2]],
                                    color = groups)))

    } else if(plottype == "density"){
        chrom_plot <- lapply(seq_along(plotorder), function(x){
            geom_density(aes(x = dt[[x]], color = groups[x]))
        })
    } else if(plottype == "line"){
        chrom_plot <- lapply(seq_along(plotorder), function(x){
            geom_line(aes(x = dt[[x]][,1], y = dt[[x]][,2], color = groups[x]))
        })
    } else if(plottype == "column"){
        chrom_plot <- list(geom_col(aes(x = plabs, y = dt, fill = plabs)))
    }
    # Unites the plots into a single plot
    for(i in seq_along(chrom_plot)){
        chrom_base <- chrom_base + chrom_plot[[i]]
    }
    # Applies selected transformations to selected axis
    if(!is.na(ytrans)){
        chrom_base <- chrom_base + scale_y_continuous(trans = ytrans) +
        ylab(paste0("Intensity ", "(", ytrans, ")"))
    }
    if(!is.na(xtrans)){
        chrom_base <- chrom_base + scale_x_continuous(trans = xtrans) +
        xlab(paste0("Intensity ", "(", xtrans, ")"))
    }
    # Extra plotting options
    if(origin_lines){
        chrom_base <- chrom_base + geom_hline(yintercept = 0, color = "gray70")+
            geom_vline(xintercept = 0, color = "gray70")
    }
    if(smoothing){
        chrom_base <- chrom_base +
            geom_smooth(aes(x = dt[[1]], y = dt[[2]], colour = groups))
    }
    return(chrom_base)
}

findOrder <- function(object, pheno_var){
    samp <- sub(pattern = "\\.[mM][zZ]?[xX][mM][lL]$", replacement = "",
                x =  basename(fileNames(object)))
    ordering <- phenoData(object)@data[[pheno_var]]
    return(order(sapply(samp, function(x){ which(ordering == x)})))
}
