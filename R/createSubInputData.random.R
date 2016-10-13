#' createSubInputData.random
#'
#' Create input sub datasets : random
#'
#' @param outDirPath output directory path where all ScaleNet outputs are saved
#' @param ioSubEnv a global environment variable

createSubInputData.random <- function( ioSubEnv, outDirPath ){

  doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
  tmpOut <- foreach(iEigen=c(2:ncol(ioSubEnv$eigen.vectors))) %dopar%{

    # Get the G+,m sub dataset, ie. with the m variables corresponding to the m highest eigen vector elements
    # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
    subInData = file.path(outDirPath, "subInput", paste("eigenVector", iEigen, "_", ioSubEnv$subset.m,"pos.tsv", sep=''))

    subset.randIdx <- sort(sample(seq_len(ncol(ioSubEnv$allData)), ioSubEnv$subset.m), decreasing = FALSE, method = "quick")
    write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[subset.randIdx]], file = subInData,
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

    # Get the G-,m sub dataset, ie. with the m variables corresponding to the m lowest eigen vector elements
    # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
    subInData = file.path(outDirPath, "subInput", paste("eigenVector", iEigen, "_", ioSubEnv$subset.m,"neg.tsv", sep=''))

    subset.randIdx <- sort(sample(seq_len(ncol(ioSubEnv$allData)), ioSubEnv$subset.m), decreasing = FALSE, method = "quick")
    write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[subset.randIdx]],
                file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  }
}
