# Main function that processes and annotates an MSnExp with the XCMS workflow
# camera and CliqueMS
metaboProc <- function(object,polarity = "positive", peakwidth, noise,
                        snthresh, ppm, expandrt, binsize, minFraction, bw,
                        annotation = "camera", cliqsamp = NA,mergepeaks = TRUE,
                        rtadjust = TRUE, group = TRUE, fill = TRUE){
    if(dim(object)[1] == 1){
        rtadjust <- FALSE
        group <- FALSE
        fill <- FALSE
    }
    # CentWave peak picking
    centParam <- CentWaveParam(peakwidth = peakwidth, noise = noise,
                                snthresh = snthresh, ppm = ppm)
    peakDt <- findChromPeaks(object, param = centParam)

    if(mergepeaks){
        # Merging close peaks
        mergeParam <- MergeNeighboringPeaksParam(expandRt = expandrt)
        peakDt <- refineChromPeaks(peakDt, mergeParam)
    }
    if(rtadjust){
        # Adjusting RT
        peakDt <- adjustRtime(peakDt, param = ObiwarpParam(binSize = binsize))
    }
    if(group){
        # Grouping features
        groupParam <- PeakDensityParam(sampleGroups = peakDt$sample_group,
                            minFraction = minFraction, bw = bw)
    }
    if(annotation == "cliqueMS"){
        clic_an <- cliqueAnot(object = peakDt, cliqsamp = cliqsamp,
                              polarity = polarity, ppm = ppm)
    }
    if(fill){
        peakDt <- groupChromPeaks(peakDt, param = groupParam)
        # Peak filling
        peakDt <- fillChromPeaks(peakDt)
    }
    if(annotation == "none"){
        return(peakDt)
    } else if(annotation == "camera"){
        xcms_set <- as(peakDt, "xcmsSet")
        sampclass(xcms_set) <- peakDt@phenoData@data$sample_group
        cam_an <- CAMERA::annotate(xcms_set, polarity = polarity)
        return(list(peakDt, cam_an))
    } else if(annotation == "cliqueMS"){
        return(list(peakDt, clic_an))
    }
    return(peakDt)
}

# wrapper function for cliqueMS full annotation
cliqueAnot <- function(object, cliqsamp = NA, seed = 1234, ppm, polarity){
    set.seed(seed)
    if(polarity == "positive"){
        data("positive.adinfo")
        adduct_list <- positive.adinfo
    } else if(polarity == "negative"){
        data("negative.adinfo")
        adduct_list <- negative.adinfo
    }
    if(!is.na(cliqsamp)){
        object <- filterFile(object, cliqsamp)
    }
    clic_an <- getCliques(object, filter = TRUE)
    clic_an <- getIsotopes(clic_an, ppm = ppm)
    clic_an <- getAnnotation(clic_an, ppm = ppm, adinfo = adduct_list,
                            polarity = polarity)
    # if(!is.na(cliqsamp)){
    #     featureDefinitions(object)$clicIndex <- objectlapply(
    #         1:nrow(featureDefinitions(object)), function(x){
    #         which(between(clic_an@peaklist$mz,
    #                     featureDefinitions(object)[x, "mzmin"],
    #                     featureDefinitions(object)[x, "mzmax"]) &
    #                 between(ex.Adducts2@peaklist$rt,
    #                     featureDefinitions(object)[x, "rtmin"],
    #                     featureDefinitions(object)[x, "rtmax"]))
    #     })
    return(clic_an)
}
