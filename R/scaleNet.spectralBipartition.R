#' scaleNet.spectralBipartition
#'
#' Make clustering using a recursive bi-partitionning approach based on the 2nd eigenvector
#'
#' @param outDirPath output directory path where all scalenet outputs are saved
#' @param ioSubEnv a global environment variable

scaleNet.spectralBipartition <- function(ioSubEnv, sub.varNames){

  # Display the upstream cluster name
  if(ioSubEnv$verbose == TRUE){cat("# Prev nbr", ioSubEnv$specBi.nextCluster, "\n")}

  # Get an extract of the affinity matrix corresponding to the variables sub.varNames
  sub.affinity.mat <- ioSubEnv$affinity.mat[as.character(sub.varNames), as.character(sub.varNames)]

  # Compute the Laplacian matrix for these variables
  if(ioSubEnv$similarityType == "Lrw"){

    # Compute the 'normalized Laplacian matrix'
    tmp.degree = diag(apply(sub.affinity.mat, 2, sum))
    colnames(tmp.degree) = rownames(tmp.degree) = colnames(sub.affinity.mat)

    # Remove the variables not associated to any others (di == 0)
    tmp.diagNull.idx <- which(apply(tmp.degree, 1, sum) == 0)

    if(length(tmp.diagNull.idx) == ncol(tmp.degree)){

      if(ioSubEnv$verbose == TRUE){cat("# --> All variables have zero degree\n")}
      return()

    } else if(length(tmp.diagNull.idx)>0){

      if(ioSubEnv$verbose == TRUE){cat("# --> The variables [", paste(colnames(sub.affinity.mat)[tmp.diagNull.idx], collapse = ', '), "] have zero degree --> ignored!\n")}

      # Resize the degree and affinity matrix
      sub.affinity.mat <- sub.affinity.mat[-tmp.diagNull.idx,-tmp.diagNull.idx]
      tmp.degree <- tmp.degree[-tmp.diagNull.idx,-tmp.diagNull.idx]

      if(ncol(tmp.degree) < ioSubEnv$subset.m){

        if(ioSubEnv$verbose == TRUE){cat("# --> The number of nodes per subgraph is greater that the number of connected variables\n")}
        return()
      }

    } else { if(ioSubEnv$verbose == TRUE){cat("# --> [.spectralBipartition] All variables have a none zero degree\n")} }

    # Invert the degree matrix
    tmp.degree.inv <- solve(tmp.degree)

    # Compute the Lwr matrix
    sub.decompose.mat = diag(nrow(sub.affinity.mat)) - tmp.degree.inv %*% sub.affinity.mat

    # Compute eigen values and eigen vectors of Lwr and get the second eigenvector
    eigen.decompo <- eigen(sub.decompose.mat, symmetric = TRUE)
    sub.eigen.vectors <- eigen.decompo$vectors
    rownames(sub.eigen.vectors) <- rownames(sub.decompose.mat)
    tmp.second.eigenvect <- sub.eigen.vectors[, (ncol(eigen.decompo$vectors)-1)]
    # Order the coordinates of the second eigen vector
    tmp.second.eigenvect <- sort(tmp.second.eigenvect, decreasing = TRUE, method = "quick")

    # Get the left half of the variables
    left.ensVar <- names(tmp.second.eigenvect)[1:floor(length(tmp.second.eigenvect)/2)]

    if(length(left.ensVar)>ioSubEnv$subset.m){

      # Set the cluster number for these variables
      ioSubEnv$specBi.allClusters[as.character(left.ensVar)] <- (ioSubEnv$specBi.nextCluster+1)
      if(ioSubEnv$verbose == TRUE){
        cat("# --[", ioSubEnv$specBi.nextCluster, "]", paste(ioSubEnv$specBi.allClusters[as.character(left.ensVar)], collapse = ","),"\n")
      }

      # Set the next available cluster number and iterate the bipartition
      ioSubEnv$specBi.nextCluster <- (ioSubEnv$specBi.nextCluster + 1)
      scaleNet.spectralBipartition( ioSubEnv = ioSubEnv, sub.varNames = left.ensVar )

    } else {if(ioSubEnv$verbose == TRUE){cat("# --[", ioSubEnv$specBi.nextCluster, "] Left has length", length(left.ensVar), "\n")}}

    # Get the right half of the variables
    right.ensVar <- names(tmp.second.eigenvect)[(floor(length(tmp.second.eigenvect)/2)+1):length(tmp.second.eigenvect)]

    if(length(right.ensVar)>ioSubEnv$subset.m){

      # Set the cluster number for these variables
      ioSubEnv$specBi.allClusters[as.character(right.ensVar)] <- (ioSubEnv$specBi.nextCluster+1)
      if(ioSubEnv$verbose == TRUE){cat("# --[", ioSubEnv$specBi.nextCluster, "]", paste(ioSubEnv$specBi.allClusters[as.character(right.ensVar)], collapse = ","),"\n")}

      # Set the next available cluster number and iterate the bipartition
      ioSubEnv$specBi.nextCluster <- (ioSubEnv$specBi.nextCluster + 1)
      scaleNet.spectralBipartition( ioSubEnv = ioSubEnv, sub.varNames = right.ensVar )

    } else {if(ioSubEnv$verbose == TRUE){cat("# --[", ioSubEnv$specBi.nextCluster, "] Right has length", length(right.ensVar), "\n")}}


  } else {
    cat(stop("# Not implemented yet..."))
  }
}
