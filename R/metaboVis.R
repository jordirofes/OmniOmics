#'@import ggplot2
#'@importFrom plotly ggplotly
#'@import oligo
#'@import beadarray
#'@import genefilter
#'@import xcms
#'@import biosigner
#'@import xlsx
#'@import CAMERA
#'@import cliqueMS
#'@import dplyr
#'@import pmp
#'@import limma
#'@import caret
#'@import shiny
#'@import shinyWidgets
#'@import shinyFiles
#'@import SummarizedExperiment
#'
#'@export
# Creates a ggplot2 chromatogram or EIC from an XCMS object
#'@title Chromatogram plot
#'@author Jordi Rofes Herrera
#'@description Creates a ggplot2 chromatogram or EIC from an XCMS object
#'@param object A OnDiskXCMSnExp or XCMSnExp object
#'@param filenum An optional numeric vector indicating the files to plot
#'@param mz A numeric indicating the mz value to plot (requires the ppm param) or a numeric two with the mz limits to plot
#'@param ppm An optional numeric indicating the ppm error for the mz param
#'@param rtint An optional length two numeric vector indicating the retention time interval to plot
#'@param pheno_var A numeric or string indicating the variable from the phenodata
#'to use as file groups.
#'@param chromtype The type of chromatogram to create. "max" for the maximum value and "sum" for the sum of intensities.
#'@param logscale A boolean to apply a log10 transformation of the intensities
#'@param interactive A boolean indicating if the plot will be converted to a plotly
#'@return A ggplot2 or plotly line plot of the chromatogram.
#'@export
ggChromPlot <- function(object, filenum = NA, mz = NA, ppm = NA, rtint = NA, pheno_var = 1,
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
    if(any(is.na(rtint))){
        rtint <- c(min(xcms::rtime(object)), max(xcms::rtime(object)))
    }
    # Colors by sample names
    groups <- sub(object@phenoData@data[,pheno_var][filenum], pattern =
                                "\\.mzX?ML",replacement =  "")
    if(!any(is.na(filenum))){
        object <- filterFile(object, filenum)
    }
    # Chromatogram creation
    chrom_dt <- chromatogram(object, aggregationFun = chromtype, mz = mz, rt = rtint)
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
#'@title Total ion count plot
#'@author Jordi Rofes Herrera
#'@description Creates a violin/boxplot TIC plot from an XCMS object
#'@param object A OnDiskXCMSnExp or XCMSnExp object
#'@param filenum An optional numeric vector indicating the files to plot
#'@param pheno_var A numeric or string indicating the variable from the phenodata
#'to use as file groups.
#'@param pheno_filter An optiona string of a variable name in pheno_var to filter from plotting
#'@param violing A boolean indicating if a violin plot should be plot instead of a boxplot
#'@param logscale A boolean to apply a log10 transformation of the intensities
#'@param interactive A boolean indicating if the plot will be converted to a plotly
#'@param order A boolean indicating if the file order should be used for the order of the samples
#'@return A ggplot2 or plotly violin/boxplot TIC plot
#'@export
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
    if(!missing(pheno_filter)){
        if(pheno_filter %in% groups){
            filenum <- filenum[groups == pheno_filter]
            groups <- phenoData(object)@data[[pheno_var]][filenum]
        }
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
        if(pheno_var == "rownames"){
            ordering <- rownames(phenodata)
        } else{
            ordering <- phenodata[[pheno_var]]
        }
    } else{
        ordering <- pheno_var
    }
    if(!all(ordering %in% samp)){
        stop("Filenames do not agree with the phenodata variable. Change the names accordingly")
    }
    return(order(sapply(samp, function(x){ which(ordering == x)})))
}
