# Main function that processes and annotates an MSnExp with the XCMS workflow,
# camera and cliqueMS
#'@export
metaboImport <- function(filedir, filelabels){
    if(dir.exists(filedir)){
        met_files <- list.files(filedir, pattern = "\\.[mM][zZ]?[xX][mM][lL]$",
                                full.names = TRUE)
        pd <- data.frame(sample_name = basename(met_files),
                        sample_group = filelabels,
                        stringsAsFactors = FALSE)
        mz_dt <- readMSData(files = met_files, mode = "onDisk", pdata = new("NAnnotatedDataFrame",pd))
        return(mz_dt)
    } else{
        pd <- data.frame(sample_name = basename(filedir),
                         sample_group = grouplabels,
                         stringsAsFactors = FALSE)
        mz_dt <- readMSData(files = filedir, mode = "onDisk", pdata = new("NAnnotatedDataFrame",pd))
        return(mz_dt)
    }
}
#'@export
metaboProc <- function(object, polarity = "positive", groupvar, peakwidth,
                        noise, snthresh, ppm, expandrt, binsize, minFraction,
                        bw, annotation = "camera", cliqsamp = NA,
                        mergepeaks = TRUE, rtadjust = TRUE, group = TRUE,
                        fill = TRUE,summ = TRUE){
    # Changes some parameters if there is only one sample
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
        groupParam <- PeakDensityParam(
            sampleGroups = phenoData(peakDt)@data[,groupvar],
                            minFraction = minFraction, bw = bw)
        peakDt <- groupChromPeaks(peakDt, param = groupParam)

    }
    if(annotation == "cliqueMS"){
        # cliqueMS annotation
        clic_an <- cliqueAnot(object = peakDt, cliqsamp = cliqsamp,
                              polarity = polarity, ppm = ppm)
    }
    if(fill){
        # Peak filling
        peakDt <- fillChromPeaks(peakDt)
    }
    if(annotation == "none"){
        return(peakDt)
    } else if(annotation == "camera"){
        # camera annotation
        xcms_set <- as(peakDt, "xcmsSet")
        sampclass(xcms_set) <- phenoData(peakDt)@data[,groupvar]
        cam_an <- CAMERA::annotate(xcms_set, polarity = polarity)
        featureDefinitions(peakDt) <- cbind(
            featureDefinitions(peakDt), getPeaklist(cam_an)[,c("isotopes",
                                                                "adduct")])
        res <- quantify(peakDt)
        return(list(peakDt, cam_an, res))
    } else if(annotation == "cliqueMS"){

        if(!is.na(cliqsamp)){
            joint_cp <- left_join(as.data.frame(chromPeaks(peakDt)),
                                  getPeaklistanClique(clic_an))
            annot <- lapply(featureDefinitions(peakDt)$peakidx, function(x){
                l <- lapply(x, function(y){
                    if(!is.na(joint_cp[y, 13])){
                        return(joint_cp[y, 13:28])
                    }
                })
                if(all(is.na(l))){
                    return(NA)
                } else{
                    return(l)
                }
            })
            featureDefinitions(peakDt)$annotation <- annot
        }
        feat_list <- quantify(peakDt)
        return(list(peakDt, clic_an, feat_list))
    }
    feat_list <- quantify(peakDt)
    return(list(peakDt, feat_list))
}
#'@export
# wrapper function for cliqueMS full annotation
cliqueAnot <- function(object, cliqsamp = NA, seed = 1234, ppm, polarity){
    set.seed(seed)
    # cliqueMS adduct data
    if(polarity == "positive"){
        data("positive.adinfo", package = "cliqueMS")
        adduct_list <- positive.adinfo
    } else if(polarity == "negative"){
        data("negative.adinfo", package = "cliqueMS")
        adduct_list <- negative.adinfo
    }
    if(!is.na(cliqsamp)){
        object <- filterFile(object, cliqsamp)
    }
    # cliqueMS workflow
    clic_an <- cliqueMS::getCliques(object, filter = TRUE)
    clic_an <- cliqueMS::getIsotopes(clic_an, ppm = ppm)
    clic_an <- cliqueMS::getAnnotation(clic_an, ppm = ppm, adinfo = adduct_list,
                            polarity = polarity)
    return(clic_an)
}
