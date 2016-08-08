#' createSubInputData.spectral
#'
#' Create input sub datasets : spectral
#'
#' @param outDirPath output directory path where all ScaleNet outputs are saved
#' @param ioSubEnv a global environment variable

createSubInputData.spectral <- function( ioSubEnv, outDirPath ){

  doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
  tmpOut <- foreach(iEigen=c(2:ncol(ioSubEnv$eigen.vectors))) %dopar% {

    # Sort the eigen vector elements in decreasing order
    tmp.elts <- sort(ioSubEnv$eigen.vectors[,iEigen], decreasing = TRUE, method = "quick")

    # Get the G+,m sub dataset, ie. with the m variables corresponding to the m highest eigen vector elements
    # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
    subInData = file.path(outDirPath, "subInput", paste("eigenVector", iEigen, "_", ioSubEnv$subset.m,"pos.tsv", sep=''))
    write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[which(colnames(ioSubEnv$allData) %in% names(tmp.elts)[1:ioSubEnv$subset.m])]],
                file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

    # Get the G-,m sub dataset, ie. with the m variables corresponding to the m lowest eigen vector elements
    # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
    subInData = file.path(outDirPath, "subInput", paste("eigenVector", iEigen, "_", ioSubEnv$subset.m,"neg.tsv", sep=''))
    write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[which(colnames(ioSubEnv$allData) %in% names(tmp.elts)[(length(tmp.elts)-ioSubEnv$subset.m+1):length(tmp.elts)])]],
                file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  }
}
