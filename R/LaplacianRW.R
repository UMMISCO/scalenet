#' LaplacianRW
#'
#' LaplacianRW computes the LAPLACIAN MATRIX, Lrw=D^-1.L=I-D^-1.W, where W is the
#' mutual information matrix, D=diag(di) and di=sum(Coli)
#'
#' @param ioSubEnv a global environment variable

LaplacianRW <- function(ioSubEnv){

  tmp.inputFileName <- NULL; tmp.path <- NULL
  tmp.path.affinity.mat <- NULL; tmp.path.affinity.mat.time <- NULL

  if(!is.null(ioSubEnv$inputData.filePath)){

    tmp.inputFileName <- basename(ioSubEnv$inputData.filePath)
    tmp.path <- gsub(tmp.inputFileName,'', ioSubEnv$inputData.filePath)
    tmp.path <- file.path(tmp.path, "scalenet")
    if(!dir.exists(tmp.path)){dir.create(tmp.path)}

    # MAKE AFFINITY
    tmp.path.affinity.mat <- file.path(tmp.path, gsub(".txt$", "_affinity.mat.rds", tmp.inputFileName))
    tmp.path.affinity.mat.time <- file.path(tmp.path, gsub(".txt$", "_affinity.time.rds", tmp.inputFileName))
  }

  if(!file.exists(ifelse(is.null(tmp.path.affinity.mat), "", tmp.path.affinity.mat))){

    # Compute the mutual information matrix if it does no exists yet
    cat("# --> ... the mutual information matrix (", ioSubEnv$mi.estimator, "estimator)\n# ---\n")

    # --!-- Time
    time.affinity.start <- proc.time()

    # Initialize the affinity matrix
    ioSubEnv$affinity.mat <- matrix(0, ncol = ncol(ioSubEnv$allData), nrow = ncol(ioSubEnv$allData))
    colnames(ioSubEnv$affinity.mat) = rownames(ioSubEnv$affinity.mat) = colnames(ioSubEnv$allData)

    # --> compute the mutual information matrix
    ioSubEnv$affinity.mat <- minet::build.mim(as.data.frame(ioSubEnv$allData),
                                              estimator = ioSubEnv$mi.estimator)

    # --> If some values are missing (NA/Nan), set to 0 the affinity
    tmp.naArrIdx <- which(is.na(ioSubEnv$affinity.mat), arr.ind = TRUE)
    if(nrow(tmp.naArrIdx) > 0){
      ioSubEnv$affinity.mat[tmp.naArrIdx] <- 0
    }

    # --!-- Time
    ioSubEnv$affinity.mat.time <- (proc.time() - time.affinity.start)

    if(!is.null(ioSubEnv$inputData.filePath)){
      # Save the affinity matrix and computation time for futur use
      saveRDS(ioSubEnv$affinity.mat, tmp.path.affinity.mat)
      saveRDS(ioSubEnv$affinity.mat.time, tmp.path.affinity.mat.time)
    }

  } else{

    # Load the computed affinity matrix and computation time
    ioSubEnv$affinity.mat <- readRDS(tmp.path.affinity.mat)
    ioSubEnv$affinity.mat.time <- readRDS(tmp.path.affinity.mat.time)
  }
  # Add to total time (in fact set here...)
  ioSubEnv$total.time <- ioSubEnv$affinity.mat.time

  if(ioSubEnv$verbose){
    cat("\n# --| Head of MIM data:\n")
    cat("# -----------------------------\n")
    print(ioSubEnv$affinity.mat[1:5,1:5])
    print(dim(ioSubEnv$affinity.mat))
    print(ioSubEnv$affinity.mat.time)
    cat("# -----------------------------\n")
  }

  tmp.path.decompose.mat <- NULL; tmp.path.decompose.mat.time <- NULL

  if(!is.null(ioSubEnv$inputData.filePath)){
    # MAKE DECOMPOSE
    # --!-- The cached matrix depends on WHICH matrix we decompose (Lsym vs the
    # legacy asymmetric Lrw), so the variant has to be part of the file name.
    # Without this a cache written before the fix would be reloaded after it and
    # silently reinstate the column-order dependence.
    tmp.variant <- if(isTRUE(ioSubEnv$legacy.lrw)) "_decompLrwLegacy" else "_decompLsym"
    tmp.path.decompose.mat <- file.path(tmp.path, gsub(".txt$", paste0(tmp.variant, ".mat.rds"), tmp.inputFileName))
    tmp.path.decompose.mat.time <- file.path(tmp.path, gsub(".txt$", paste0(tmp.variant, ".time.rds"), tmp.inputFileName))
  }

  if(!file.exists(ifelse(is.null(tmp.path.affinity.mat), "", tmp.path.decompose.mat))){

    # --> init the matrice to be decomposed with the similarity matrix
    ioSubEnv$decompose.mat <- ioSubEnv$affinity.mat

    if(ioSubEnv$similarityType == "Lrw"){

      # --!-- Time
      time.laplacian.start <- proc.time()

      # Compute the 'normalized Laplacian matrix'
      tmp.degree = diag(apply(ioSubEnv$affinity.mat, 2, sum, na.rm = TRUE))
      colnames(tmp.degree) = rownames(tmp.degree) = colnames(ioSubEnv$allData)

      # Remove the variables not associated to any others (di == 0)
      tmp.diagNull.idx <- which(apply(tmp.degree, 1, sum) == 0)

      if(length(tmp.diagNull.idx) == ncol(tmp.degree)){

        cat("# --> All variables have zero degree\n")
        quit(save = "no", status = 0)

      } else if(length(tmp.diagNull.idx)>0){

        cat("# --> The variables [", paste(colnames(ioSubEnv$affinity.mat)[tmp.diagNull.idx],
                                           collapse = ', '), "] have zero degree --> ignored!\n")

        # Resize the degree and affinity matrix
        ioSubEnv$affinity.mat <- ioSubEnv$affinity.mat[-tmp.diagNull.idx,-tmp.diagNull.idx]
        tmp.degree <- tmp.degree[-tmp.diagNull.idx,-tmp.diagNull.idx]

        if(ncol(tmp.degree) < ioSubEnv$subset.m){

          cat("# --> The number of nodes per subgraph is greater that the number of connected variables\n")
          quit(save = "no", status = 0)
        }

      } else { cat("# --> [.LaplacianRW] All variables have a none zero degree\n") }

      if(ioSubEnv$verbose){
        cat("\n# --| Head of MIM degree:\n")
        cat("# -----------------------------\n")
        print(tmp.degree[1:5,1:5])
        print(dim(tmp.degree))
        cat("# -----------------------------\n")
      }

      # Keep the degree vector: computeEigenVectVal needs it to map the
      # eigenvectors of Lsym back to those of Lrw.
      ioSubEnv$decompose.degree <- diag(tmp.degree)
      tmp.n <- nrow(ioSubEnv$affinity.mat)
      tmp.names <- colnames(ioSubEnv$affinity.mat)

      if(isTRUE(ioSubEnv$legacy.lrw)){

        # --!-- LEGACY, column-order dependent. Lrw = I - D^-1.W is NOT
        # symmetric, yet computeEigenVectVal decomposes it with
        # eigen(symmetric = TRUE), which reads only the LOWER triangle: entry
        # (i,j) is taken as -W[i,j]/D[i] with i the larger index. Permuting the
        # input columns therefore changes the matrix being decomposed, not just
        # its row/column order, so the eigenvectors, the variable subsets and
        # every reconstructed edge change with it. Kept only to regenerate
        # previously published results.
        tmp.degree.inv <- solve(tmp.degree)
        ioSubEnv$decompose.mat = (diag(tmp.n) - tmp.degree.inv %*% ioSubEnv$affinity.mat)

      } else {

        # Lsym = I - D^-1/2.W.D^-1/2 is symmetric by construction, so permuting
        # the input permutes the matrix and nothing else. It has the same
        # eigenvalues as Lrw, and computeEigenVectVal recovers Lrw's
        # eigenvectors as v = D^-1/2.u (since Lrw = D^-1/2.Lsym.D^1/2).
        tmp.dh <- diag(1/sqrt(ioSubEnv$decompose.degree), nrow = tmp.n)
        ioSubEnv$decompose.mat = (diag(tmp.n) - tmp.dh %*% ioSubEnv$affinity.mat %*% tmp.dh)
      }
      dimnames(ioSubEnv$decompose.mat) <- list(tmp.names, tmp.names)

      # --!-- Time
      ioSubEnv$decompose.mat.time <- (proc.time() - time.laplacian.start)

      if(!is.null(ioSubEnv$inputData.filePath)){
        # Save the decompose matrix and computation time for futur use
        saveRDS(ioSubEnv$decompose.mat, tmp.path.decompose.mat)
        saveRDS(ioSubEnv$decompose.mat.time, tmp.path.decompose.mat.time)
      }

    }else{cat(stop("# Not implemented yet..."))}

  } else {

    # Load the computed decompose matrix and computation time
    ioSubEnv$decompose.mat <- readRDS(tmp.path.decompose.mat)
    ioSubEnv$decompose.mat.time <- readRDS(tmp.path.decompose.mat.time)
  }
  # Add to total time
  ioSubEnv$total.time <- ioSubEnv$total.time + ioSubEnv$decompose.mat.time

  if(ioSubEnv$verbose){
    cat("\n# --| Head of Lwr:\n")
    cat("# -----------------------------\n")
    print(ioSubEnv$decompose.mat[1:5,1:5])
    print(dim(ioSubEnv$decompose.mat))
    print(ioSubEnv$decompose.mat.time)
    cat("# -----------------------------\n")
  }
}
