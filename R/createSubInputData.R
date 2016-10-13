#' createSubInputData
#'
#'
#' @param ioSubEnv a global environment variable

createSubInputData <- function( ioSubEnv ){

  # Create all the sub input datasets
  tmp.dirPath <- file.path(ioSubEnv$output.dirPath, "subGraphs")
  if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}
  if(!dir.exists(file.path(tmp.dirPath, "subInput"))){dir.create(file.path(tmp.dirPath, "subInput"))}

  # --!-- Time
  time.multiReconstruction.input.start <- proc.time()

  if(ioSubEnv$subset.select == "spectral"){

    cat("# ... using subset from the ordered eigen vectors...\n")
    createSubInputData.spectral(ioSubEnv = ioSubEnv, outDirPath = tmp.dirPath)

  } else if(ioSubEnv$subset.select == "random"){

    cat("# ... using random subsets...\n")
    createSubInputData.random(ioSubEnv = ioSubEnv, outDirPath = tmp.dirPath)

  } else if(ioSubEnv$subset.select == "spectralKmeans"){

    cat("# ... using subset from the spectral clustering...\n")
    createSubInputData.spectralKmeans(ioSubEnv = ioSubEnv, outDirPath = tmp.dirPath)

  } else if(ioSubEnv$subset.select == "spectralFuzzyCmeansOrder"){

    cat("# ... using subset from the spectral fuzzy cmeans clustering (order)...\n")
    createSubInputData.spectralFuzzyCmeansOrder(ioSubEnv = ioSubEnv, outDirPath = tmp.dirPath)

  } else if(ioSubEnv$subset.select == "spectralFuzzyCmeansSample"){

    cat("# ... using subset from the spectral fuzzy cmeans clustering (sample)...\n")
    createSubInputData.spectralFuzzyCmeansSample(ioSubEnv = ioSubEnv, outDirPath = tmp.dirPath)

  } else if(ioSubEnv$subset.select == "spectralBipartition"){

    cat("# ... using subset from the spectral bi partitionning...\n")
    createSubInputData.spectralBipartition(ioSubEnv = ioSubEnv, outDirPath = tmp.dirPath)

  } else { stop("# Not implemented yet...") }

  # --!-- Time
  ioSubEnv$inputs.time <- (proc.time() - time.multiReconstruction.input.start)
  ioSubEnv$total.time <- ioSubEnv$total.time + ioSubEnv$inputs.time
}
