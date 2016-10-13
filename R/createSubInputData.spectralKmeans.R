#' createSubInputData.spectralKmeans
#'
#' Create input sub datasets : spectral K-means
#'
#' @param outDirPath output directory path where all ScaleNet outputs are saved
#' @param ioSubEnv a global environment variable

createSubInputData.spectralKmeans <- function( ioSubEnv, outDirPath ){

  # Perform a kmeans on the k first eigen vectors
  # set.seed(numSeed)
  res.pamk <- try(fpc::pamk(ioSubEnv$eigen.vectors, krange = ioSubEnv$kmeans.k, criterion="asw", critout = TRUE))

  if(!inherits(res.pamk, "try-error")){

    # Update the mean, median and sd of subset sizes
    tmp.table <- table(res.pamk$pamobject$clustering)
    ioSubEnv$subset.m.mean <- mean(tmp.table)
    ioSubEnv$subset.m.med <- median(tmp.table)
    ioSubEnv$subset.m.sd <- sd(tmp.table)
    ioSubEnv$cluster.asw <- res.pamk$crit[ioSubEnv$kmeans.k]

    # iClust = 2
    doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
    tmpOut <- foreach(iClust=seq_len(ioSubEnv$kmeans.k)) %dopar%{

      # Get the Gclust,m sub dataset, ie. with the ~m variables corresponding to the m clustered variables
      # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
      subInData <- file.path(outDirPath, "subInput", paste("spectralKmeans", iClust, "_", ioSubEnv$subset.m,".tsv", sep=''))

      tmp.elts <- res.pamk$pamobject$clustering[res.pamk$pamobject$clustering==iClust]
      write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[which(colnames(ioSubEnv$allData) %in% names(tmp.elts))], drop = FALSE],
                  file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
    }
  } else { warning("# --W-- k-means not possible") }
}
