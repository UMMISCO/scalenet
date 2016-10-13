#' createSubInputData.spectralFuzzyCmeansOrder
#'
#' Create input sub datasets : spectral fuzzy cmeans (order)
#' - make 2*subset.k fuzzy clusters
#' - for each fuzzy cluster, keep only the m variables with the strongest membership
#'
#' @param outDirPath output directory path where all ScaleNet outputs are saved
#' @param ioSubEnv a global environment variable

createSubInputData.spectralFuzzyCmeansOrder <- function( ioSubEnv, outDirPath ){

  # Perform a cmeans on the k first eigen vectors
  # Ask for 2*subset.k clusters, and
  # take the m variables with the best membership from each of these clusters
  res.cmeans <- try(e1071::cmeans(ioSubEnv$eigen.vectors, centers = ioSubEnv$cmeans.k, dist = "manhattan",
                                  method = "cmeans", m = 2, iter.max = 100), silent = TRUE)

  if(!inherits(res.cmeans, "try-error")){

    # iClust = 3
    doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
    tmpOut <- foreach(iClust=seq_len(ioSubEnv$cmeans.k)) %dopar%{

      # Keep the m first variables
      tmp.elts <- (sort(res.cmeans$membership[,iClust], decreasing = T, method = "quick"))[c(1:ioSubEnv$subset.m)]

      # Get the Gclust,m sub dataset, ie. with the ~m variables corresponding to the m clustered variables
      # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
      subInData <- file.path(outDirPath, "subInput", paste("spectralFuzzyCmeansOrder", iClust, "_", ioSubEnv$subset.m,".tsv", sep=''))

      # tmp.elts <- res.pamk$pamobject$clustering[res.pamk$pamobject$clustering==iClust]
      write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[which(colnames(ioSubEnv$allData) %in% names(tmp.elts))], drop = FALSE],
                  file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
    }
  } else { warning("# --W-- No possible fuzzy c-means") }
}
