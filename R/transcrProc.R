#'@title Microarray intensity distribution boxplot/violin plot
#'@author Jordi Rofes Herrera
#'@description Creates a boxplot/violin plot from a microarray object with `exprs()` such as an ExpressionSet
#'@param object A microarray data in a ExpressionSet or GeneFeatureSet
#'@param violin A boolean indicating if the plot will be a violin plot/boxplot
#'@param groupvar A numeric indicating the column to use as grouping factor or a character indicating it's name
#'@param interactive A boolean indicating if the plot will be converted to an interactive `ggplotly()`
#'@param nsamp The number of samples used in the plotting to reduce computation time
#'@param ytrans A transformation applied to the represented intensities, usually log2 for microarray analysis
#'@return Returns a ggplot object or a ggplotly
#'@examples
#'plot_crayons()
#'@export
ggExprDistrPlot <- function(object, violin, groupvar, interactive,
                            nsamp = 100000, ytrans = "log2"){
    dt <- exprs(object)[sample(x = seq_len(nrow(exprs(object))),size = nsamp),]
    dt <- as.list.data.frame(as.data.frame(dt))
    plabs <- names(dt)
    groups <- object@phenoData@data[,groupvar]
    plottype <- ifelse(violin == TRUE, "violin", "boxplot")
    p <- ggStandardPlot(dt = dt,plabs = plabs, groups = groups,
                    plottype = plottype, ytrans = ytrans,
                    ptitle = "Array expression distribution", xlab = "Sample",
                    ylab = "Intensity")
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@title Microarray intensity density plot
#'@author Jordi Rofes Herrera
#'@description Creates a density distribution plot from a microarray object with `exprs()` such as an ExpressionSet
#'@param object A microarray data in a ExpressionSet or GeneFeatureSet
#'@param groupvar A numeric indicating the column to use as grouping factor or a character indicating it's name
#'@param interactive A boolean indicating if the plot will be converted to an interactive `ggplotly()`
#'@param nsamp The number of samples used in the plotting to reduce computation time
#'@param xtrans A transformation applied to the represented intensities, usually log2 for microarray analysis
#'@return Returns a ggplot object or a ggplotly
#'@examples
#'plot_crayons()
#'@export
ggDensityPlot <- function(object, groupvar, interactive, nsamp = 10000,
                          ytrans = NA, xtrans = "log2"){
    dt <- exprs(object)[sample(x = seq_len(nrow(exprs(object))),size = nsamp),]
    dt <- as.list.data.frame(as.data.frame(dt))
    plabs <- names(dt)
    groups <- phenoData(object)[[groupvar]]
    browser()
    p <- ggStandardPlot(dt = dt, plabs = plabs, groups = groups,
                   plottype = "density", ytrans = ytrans, xtrans = xtrans,
                   ptitle = "Expression density distribution", ylab = "Density",
                   xlab = "Values")
    if(interactive){
        return(ggplotly(p))
    }
    return(p)
}
#'@title Microarray sample PCA plot
#'@author Jordi Rofes Herrera
#'@description Creates a PCA plot from a microarray object with `exprs()` such as an ExpressionSet
#'@param object A microarray data in a ExpressionSet or GeneFeatureSet
#'@param pc A vector of length two indicanting the principal components to plot on each axis
#'@param groupvar A numeric indicating the column to use as grouping factor or a character indicating it's name
#'@param interactive A boolean indicating if the plot will be converted to an interactive `ggplotly()`
#'@param scale A boolan indicating if the PCA will be scaled, (the default is TRUE).
#'@param nsamp The number of samples used in the plotting to reduce computation time
#'@return Returns a ggplot object or a ggplotly
#'@examples
#'plot_crayons()
#'@export
ggpcaPlot <- function(object, pc = c(1,2), groupvar,
                        scale = TRUE, interactive = TRUE){
    if(length(pc) != 2){
        pc <- c(1,2)
    }
    dt <- exprs(object)
    pcdt <- prcomp(t(dt), scale = scale)
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

#'@title Wrapper function for .idat, .cel and geoDatasets importing
#'@author Jordi Rofes Herrera
#'@description Imports .idat, .cel files from microarray expression data into their respective formats and the download of geoDatasets.
#'To import illumina datasets it is required the appropiate annotation package
#'@param datapath A string with the directory where the files are located or a vector with all the files to import
#'@param sampledata A string indicating the phenodata location as a comma separated values or excel file.
#'@param geoid A string indicating a geoDataset accesion entry
#'@param header A boolean indicating if the sampledata has a header with the variable names
#'@param sep A string indicating the separator of the sampledata file if it's a .csv file, usually a space or a comma.
#'@param groupvar A numeric indicating the column to use as grouping factor or a character indicating it's name
#'@return Returns an object with the expression data.
#'@examples
#'plot_crayons()
#'@export
importRawTranscript <- function(datapath, sampledata, geoid, header = TRUE, sep = ","){
    if(!is.na(geoid)){
        dt <- GEOquery::getGEO(geoid)
        if(dim(featureData(dt[[1]]))[2] == 0){
            stop("The selected accession entry data has no features")
        } else{
            return(dt)
        }
    }
    if(!missing(sampledata)){
        if(is.character(sampledata)){
            if(grepl("\\.[cC][sS][vV]$", sampledata)){
                sampledata <- read.csv(sampledata, header = header, sep = sep)
            } else if(grepl("\\.[xX][lL][sS]?[xX]$", sampledata)){
                sampledata <- read_excel(sampledata)
            }
        }
    }
    if(dir.exists(pathfile)){
        if(length(list.celfiles(pathfile)) != 0){
            if(missing(sampledata)){
                files <- list.files(datapath)
                if(grepl("\\.[cC][sS][vV]$", files)){
                    sampledata <- read.csv(files[grepl("\\.[cC][sS][vV]$",
                                                        files)],
                                            header = header, sep = sep)
                } else if(grepl("\\.[xX][lL][sS]?[xX]$", files)){
                    sampledata <- read_excel(
                        files[grepl("\\.[xX][lL][sS]?[xX]$")])
                }
            }
            return(read.celfiles(list.celfiles(pathfile, full.names = TRUE),
                                phenoData = sampledata))
        }
        if(length(list.idatfiles(pathfile)) != 0){
            return(readIdatFiles(list.idatfiles(pathfile, full.names = TRUE)))
        }
    } else{
        if(any(grepl("\\.[cC][eE][lL]$", pathfile))){
            return(read.celfiles(pathfile, phenoData = sampledata))
        }
        if(any(grepl("\\.[iI][dD][aA][tT]$", pathfile))){
            return(readIdatFiles(pathfile))
        }
    }
}
# List of .IDAT files from a directory
list.idatfiles <- function(...){
    dt_files <- list.files(...)
    return(dt_files[grepl("\\.[iI][dD][aA][tT]$", dt_files)])
}

#'@title Processing function for illumina and affymetrix microarrays illumina expression normalization and RMA for affymetrix.
#'@author Jordi Rofes Herrera
#'@description Normalizes a GeneFeatureSet using the RMA method or an ExpressionSetIllumina using the selected methodology
#'@param object A microarray data in a ExpressionSet or GeneFeatureSet
#'@param method String specifying the illumina normalization method methods are:
#'"quantile", "qspline", "vsn", "rankInvariant", "median" and "none".
#'@return Returns a processed ExpressionSetIllumina for illumina data or an ExpressionSet
#'@examples
#'plot_crayons()
#'@export
setGeneric("procTranscript", function(object, method){
    standardGeneric("procTranscript")
})
#'@export
setMethod("procTranscript", "ExpressionSetIllumina",
            function(object, method = "quantile"){
    return(normaliseIllumina(object, method = method))
})
#'@export
setMethod("procTranscript", "GeneFeatureSet", function(object){
    return(oligo::rma(object))
})
