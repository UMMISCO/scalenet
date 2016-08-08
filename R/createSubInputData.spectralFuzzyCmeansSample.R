#' createSubInputData.spectralFuzzyCmeansSample
#'
#' Create input sub datasets : spectral fuzzy cmeans (sample)
#' - make 1/m clusters
#' - assign each variable to a cluster based on its membership
#' (use 'sample' based on the membership as probabilities and perform 20 assignments)

#'
#' @param outDirPath output directory path where all ScaleNet outputs are saved
#' @param ioSubEnv a global environment variable

createSubInputData.spectralFuzzyCmeansSample <- function( ioSubEnv, outDirPath ){

  # Perform a cmeans on the k first eigen vectors
  # Ask for 1/m clusters, and
  # assign each variable to a cluster based on its membership
  res.cmeans <- try(e1071::cmeans(ioSubEnv$eigen.vectors, centers = ioSubEnv$cmeans.k, dist = "manhattan",
                                  method = "cmeans", m = 2, iter.max = 100), silent = TRUE)

  if(!inherits(res.cmeans, "try-error")){

    # iNbVarSpl = 1
    for(iNbVarSpl in seq_len(ioSubEnv$cmeans.nbVarSample)){

      # Make ioSubEnv$cmeans.nbVarSample assignments
      tmp.cluster.members.list <- vector("list", ioSubEnv$cmeans.k)

      # iVar = 1
      for(iVar in seq_len(nrow(res.cmeans$membership))){
        var.membership <- as.vector(as.matrix(res.cmeans$membership[iVar,]))
        var.assignedCluster <- sample(seq_len(length(var.membership)), 1, prob = var.membership)
        tmp.cluster.members.list[[as.numeric(var.assignedCluster)]] <- c(tmp.cluster.members.list[[as.numeric(var.assignedCluster)]], rownames(res.cmeans$membership)[iVar])
      }

      # iClust = 3
      doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
      tmpOut <- foreach(iClust=seq_len(length(tmp.cluster.members.list))) %dopar%{

        # Get the member variable name
        tmp.elts  <- tmp.cluster.members.list[[iClust]]

        # Get the Gclust,m sub dataset, ie. with the ~m variables corresponding to the m clustered variables
        # --!!--> make sure the pairs of variables are in the same order as in the original dataset (splData)
        subInData <- file.path(outDirPath, "subInput", paste("spectralFuzzyCmeansSample", iClust, "_", ioSubEnv$subset.m, "_", iNbVarSpl, ".tsv", sep=''))

        # tmp.elts <- res.pamk$pamobject$clustering[res.pamk$pamobject$clustering==iClust]
        write.table(ioSubEnv$allData[, colnames(ioSubEnv$allData)[which(colnames(ioSubEnv$allData) %in% tmp.elts)], drop = FALSE],
                    file = subInData, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
      }
    }

  } else { warning("# --W-- No possible fuzzy c-means") }

}
