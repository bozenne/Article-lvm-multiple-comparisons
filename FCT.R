## * sinkDirectory
sinkDirectory <- function (path, string.keep = NULL, string.exclude = NULL, addMissingCol = FALSE, 
                           fixed = FALSE, trace = TRUE){
    if (!dir.exists(path)) {
        possible.dirs <- setdiff(list.dirs(file.path(path, ".."), 
                                           full.names = FALSE), "")
        stop("Directory ", path, " does not exists \n", "Existing dir. in parent dir.: \"", 
             paste0(possible.dirs, collapse = "\" \""), "\"\n")
    }
    allFiles <- list.files(path)
    if (length(allFiles) == 0) {
        warning("Empty directory \n")
        return(NULL)
    }
    index.file <- 1:length(allFiles)
    if (!is.null(string.keep)) {
        index.file <- intersect(index.file, grep(string.keep, 
                                                 allFiles, fixed = fixed))
    }
    if (!is.null(string.exclude)) {
        index.file <- setdiff(index.file, grep(string.exclude, 
                                               allFiles, fixed = fixed))
    }
    n.files <- length(index.file)
    if (n.files == 0) {
        warning("No file matching the query \n")
        return(NULL)
    }
    dt.merge <- NULL
    if (trace) {
        cat("read ", n.files, " files \n", sep = "")
        pb <- txtProgressBar(max = n.files)
    }
    for (iFile in 1:n.files) {
        iFilename <- allFiles[index.file[iFile]]
        iRes <- data.table::as.data.table(readRDS(file = file.path(path, 
                                                                   iFilename)))
        if (NROW(iRes) == 0) {
            next
        }
        iRes[, `:=`(iFile, iFilename)]
        if (!is.null(dt.merge) && addMissingCol == TRUE && NCOL(dt.merge) != 
            NCOL(iRes)) {
            missing <- setdiff(names(dt.merge), names(iRes))
            if (!is.null(missing)) {
                for (iMiss in missing) {
                    vec.tempo <- dt.merge[1, .SD, .SDcols = iMiss][[1]]
                    vec.tempo[1] <- NA
                    iRes[, `:=`(c(iMiss), vec.tempo[1])]
                }
            }
            missing <- setdiff(names(iRes), names(dt.merge))
            if (!is.null(missing)) {
                for (iMiss in missing) {
                    vec.tempo <- iRes[1, .SD, .SDcols = iMiss][[1]]
                    vec.tempo[1] <- NA
                    dt.merge[, `:=`(c(iMiss), vec.tempo[1])]
                }
            }
        }
        dt.merge <- rbind(dt.merge, iRes)
        if (trace) {
            setTxtProgressBar(pb, iFile)
        }
    }
    if (trace) {
        close(pb)
    }
    return(dt.merge)
}
