
# Creates a ggplot2 chromatogram or EIC from an XCMS object
ggChromPlot <- function(object, filenum, mzfilter, chromtype){
    if(any(is.null(mzfilter))){
        mzfilter <- c(-Inf, Inf)
    } else{
        if(any(is.na(mzfilter))){
            mzfilter <- c(-Inf, Inf)
        }
    }

    # Colors by sample names
    group_colors <- seq_along(object@phenoData@data$sample_name)
    names(group_colors) <- sub(object@phenoData@data$sample_name, pattern =
                                ".mzXML",replacement =  "", fixed = TRUE)
    names(group_colors) <- sub(x = names(group_colors), pattern = ".mzML",
                                replacement =  "", fixed = TRUE)
    if(!is.na(filenum)){
        object <- filterFile(object, filenum)
    }
    # Chromatogram creation
    chrom_dt <- chromatogram(object,
                             aggregationFun = chromtype, mz = mzfilter)
    # Plot creation
    chrom_plot <- lapply(seq_along(filenum), function(x){
        geom_line(aes(x = rtime(chrom_dt[1,x]), y = intensity(chrom_dt[1,x]),
        color = names(group_colors[filenum[x]])))
    })
    chrom_base <- ggplot() + theme_minimal() + xlab("Retention Time (RT)") +
        ylab("Intensity") + theme(axis.line = element_line(colour = "black",
                                                            size = 0.5,
                                                            linetype = "solid"),
                                    legend.title = element_blank())
    for(i in seq_along(chrom_plot)){
        chrom_base <- chrom_base + chrom_plot[[i]]
    }

    if(logscale){chrom_base <- chrom_base + scale_y_log10() +
        ylab("Intensity (log10)")}

    if(any(is.infinite(mzfilter))){
        mz_dt <- round(mz(chrom_dt[,1]), 4)
    } else{mz_dt <- mzfilter}

    p <- (ggplotly(chrom_base) %>%
        layout(title = paste(mz_dt[1], "-", mz_dt[2], "mz",collapse = ""),
                legend = list(x = 100, y = 1,
                            title = list(text = "Experimental\nCondition",
                                        x = 100, y = 1))))
    return(p)
}

# Creates a violing TIC plot from an XCMS object
ggTicViolin <- function(object, filenum){


bpdt <- split(tic(filterFile(object, filenum)),
              f = fromFile(filterFile(object, filenum)))

chrom_plot <- lapply(seq_along(bpdt), function(x){
    geom_violin(aes(x = x, y = bpdt[[x]]))
})

mz_dt <- round(mz(chrom_dt[,1]), 4)

chrom_base <- ggplot() + theme_minimal() + xlab("Sample") +
    ylab("Intensity") + theme(axis.line = element_line(colour = "black",
                                                        size = 0.5,
                                                        linetype = "solid"),
                                legend.title = element_blank())


for(i in seq_along(chrom_plot)){
    chrom_base <- chrom_base + chrom_plot[[i]]
}
if(logscale){chrom_base <- chrom_base + scale_y_log10() +
    ylab("Intensity (log10)")}

p <- ggplotly(chrom_base) %>%
            layout(title = "Total Ion Current plot per sample",
                    legend = list(x = 100, y = 1,
                                title = list(text = "Experimental\nCondition",
                                            x = 100, y = 1)))
return(p)
}



