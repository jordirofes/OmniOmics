#'@import ggplot2
#'@importFrom plotly ggplotly
#'@import oligo
#'@import beadarray
#'@importFrom genefilter nsFilter
#'@import xcms
#'@import biosigner
#'@import xlsx
#'@import CAMERA
#'@import cliqueMS
#'@import dplyr
#'
#'@export
# Creates a ggplot2 chromatogram or EIC from an XCMS object
ggChromPlot <- function(object, filenum = NA, mz = NA, ppm = NA, pheno_var = 1,
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
    } else{
        logscale <- NA
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
        cbind(chrom_dt[[x]]@rtime, chrom_dt[[x]]@intensity)})
    p <- ggStandardPlot(dt = dt, groups = groups,
                        plottype = "line", ytrans = logscale,
                        ptitle = paste(mz[1], "-", mz[2], "mz",
                                        collapse = ""), xlab = "rt",
                        ylab = "Intensity", angle = 0)
    # if(class(object)[1] == "XCMSnExp"){
    #     if(hasChromPeaks(object)){
    #
    #     }
    # }
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@export
# Creates a violin/boxplot TIC plot from an XCMS object
ggTicQuality <- function(object, filenum = NA, pheno_var = 2, pheno_filter,
                        logscale = TRUE, violin = FALSE, interactive = TRUE,
                        order = FALSE){
    # Data preparation
    if(any(is.na(filenum))){
        filenum <- seq_along(fileNames(object))
    }
    if(logscale){
        logscale <- "log10"
    } else{
        logscale <- NA
    }
    groups <- phenoData(object)@data[[pheno_var]][filenum]
    if(!missing(pheno_filter) & (pheno_filter %in% groups)){
        filenum <- filenum[groups == pheno_filter]
        groups <- phenoData(object)@data[[pheno_var]][filenum]
    }
    if(order){
        plabs <- seq_along(fileNames(object))
    } else{
        plabs <- sub(pattern = "\\.[mM][zZ][xX]?[mM][lL]$", replacement = "",
                    x =  basename(fileNames(object)))[filenum]
    }
    # Creates de data
    bpdt <- split(tic(filterFile(object, filenum)),
                  f = fromFile(filterFile(object, filenum)))
    plottype <- ifelse(violin == TRUE, "violin", "boxplot")
    p <- ggStandardPlot(dt = bpdt, plabs = plabs, groups = groups,
                        plottype = plottype, ytrans = logscale,
                        ptitle = "Total Ion Count plot", xlab = "Sample",
                        ylab = "Intensity")
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@export
ggStandardPlot <- function(dt, plabs, groups, plottype, ptitle, ytrans = NA,
                            xtrans = NA, xlab, ylab, angle = 90,
                            origin_lines = FALSE, smoothing = FALSE){
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
        chrom_plot <- lapply(seq_along(dt), function(x){
            geom_violin(aes(x = plabs[x], y = dt[[x]], fill = groups[x]))
        })
    } else if(plottype == "boxplot"){
        chrom_plot <- lapply(seq_along(dt), function(x){
            geom_boxplot(aes(x = plabs[x], y = dt[[x]],
                                fill = groups[x]))
        })
    } else if(plottype == "scatter"){
        chrom_plot <- list(geom_point(aes(x = dt[[1]], y = dt[[2]],
                                    color = groups)))

    } else if(plottype == "density"){
        chrom_plot <- lapply(seq_along(dt), function(x){
            geom_density(aes(x = dt[[x]], color = groups[x]))
        })
    } else if(plottype == "line"){
        chrom_plot <- lapply(seq_along(dt), function(x){
            geom_line(aes(x = dt[[x]][,1], y = dt[[x]][,2], color = groups[x]))
        })
    } else if(plottype == "column"){
        chrom_plot <- list(geom_col(aes(x = plabs, y = dt, fill = plabs)))
    }
    # Unites the plots into a single plot
    for(i in seq_along(chrom_plot)){
        chrom_base <- chrom_base + chrom_plot[[i]]
    }
    if(!is.na(ytrans)){
        chrom_base <- chrom_base + ylab(paste0("Intensity ", "(", ytrans, ")")) +
            scale_y_continuous(trans = ytrans)
    }
    if(!is.na(xtrans)){
        chrom_base <- chrom_base + xlab(paste0("Intensity ", "(", xtrans, ")")) +
            scale_x_continuous(trans = xtrans)
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
#'@export
findOrder <- function(files, phenodata, pheno_var){
    samp <- sub(pattern = "\\.[mM][zZ][xX]?[mM][lL]$", replacement = "",
                x =  basename(files))
    if(length(pheno_var) == 1){
        ordering <- phenodata[[pheno_var]]
    } else{
        ordering <- pheno_var
    }
    if(!all(samp %in% ordering)){
        stop("Filenames do not agree with the phenodata variable. Change the names accordingly")
    }
    return(order(sapply(samp, function(x){ which(ordering == x)})))
}
