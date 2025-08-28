#' createSubInputData.spectralBipartition
#'
#' Create input sub datasets : spectral bi partitionning
#'
#' @param outDirPath output directory path where all scalenet outputs are saved
#' @param ioSubEnv a global environment variable

createSubInputData.spectralBipartition <- function( ioSubEnv, outDirPath ){

  # Consider the variables with deg>0
  # Set a cluster vector with all the variables in the cluster 1
  ioSubEnv$specBi.allClusters <- rep(1, ncol(ioSubEnv$affinity.mat))
  names(ioSubEnv$specBi.allClusters) <- colnames(ioSubEnv$affinity.mat)
  ioSubEnv$specBi.nextCluster <- 1

  # Make a recursive bi partitionning by splitting each sorted 2nd eigenvector in 2
  scaleNet.spectralBipartition(ioSubEnv = ioSubEnv, sub.varNames = names(ioSubEnv$specBi.allClusters))
  tmp.specBi.allClusterLabels <- unique(ioSubEnv$specBi.allClusters)
  # tmp.specBi.allClusterLabels

  # iClust = 7
  doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
  tmpOut <- foreach(iClust=tmp.specBi.allClusterLabels) %dopar%{
    print(iClust)
    # Keep the partition variables
    tmp.elts <- ioSubEnv$specBi.allClusters[which(ioSubEnv$specBi.allClusters==iClust)]

    # Get the Gclust,m sub dataset, ie. with the ~m variables corresponding to the m clustered variables
    # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
    subInData <- file.path(outDirPath, "subInput", paste("spectralBipartition", iClust, "_", ioSubEnv$subset.m,".tsv", sep=''))

    # tmp.elts <- res.pamk$pamobject$clustering[res.pamk$pamobject$clustering==iClust]
    write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[which(colnames(ioSubEnv$allData) %in% names(tmp.elts))], drop = FALSE],
                file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  }

}
