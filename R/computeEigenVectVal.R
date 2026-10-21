#' computeEigenVectVal
#'
#' computeEigenVectVal computes the eigen vectors and eigen values from the laplacian matrix.
#' Only a percentage of eigen vectors is kepts (either following input arguments or automatically
#' with an elbow-like heuristic)
#'
#' @param ioSubEnv a global environment variable

computeEigenVectVal <- function(ioSubEnv){

  tmp.inputFileName <- NULL; tmp.path <- NULL
  tmp.path.eigen.vect <- NULL; tmp.path.eigen.values <- NULL; tmp.path.eigen.time <- NULL

  if(!is.null(ioSubEnv$inputData.filePath)){
    tmp.inputFileName <- basename(ioSubEnv$inputData.filePath)
    tmp.path <- gsub(tmp.inputFileName,'', ioSubEnv$inputData.filePath)
    tmp.path <- file.path(tmp.path, "scalenet")
    if(!dir.exists(tmp.path)){dir.create(tmp.path)}

    # MAKE EIGEN
    # --!-- As in LaplacianRW, the cache key must carry the variant: eigenvectors
    # of Lsym mapped back to Lrw are not the eigenvectors the legacy path
    # produces, so a cache from one must never be reloaded for the other.
    tmp.variant <- if(isTRUE(ioSubEnv$legacy.lrw)) "_lrwLegacy" else "_lsym"
    tmp.path.eigen.vect <- file.path(tmp.path, gsub(".txt$", paste0("_eigenvect", tmp.variant, ".rds"), tmp.inputFileName))
    tmp.path.eigen.values <- file.path(tmp.path, gsub(".txt$", paste0("_eigenval", tmp.variant, ".rds"), tmp.inputFileName))
    tmp.path.eigen.time <- file.path(tmp.path, gsub(".txt$", paste0("_eigen", tmp.variant, ".time.rds"), tmp.inputFileName))
  }

  if(!file.exists(ifelse(is.null(tmp.path.eigen.vect), "", tmp.path.eigen.vect))){

    # --!-- Time
    time.eigen.start <- proc.time()

    # Compute eigen values and eigen vectors of Lwr
    eigen.decompo <- eigen(ioSubEnv$decompose.mat, symmetric = TRUE)

    ioSubEnv$eigen.values <- eigen.decompo$values
    ioSubEnv$eigen.vectors <- eigen.decompo$vectors
    rownames(ioSubEnv$eigen.vectors) <- rownames(ioSubEnv$decompose.mat)
    colnames(ioSubEnv$eigen.vectors) <- paste("e", seq_len(ncol(ioSubEnv$eigen.vectors)), sep='')

    # If laplacian, sort the values and vectors in increasing order of values
    if(ioSubEnv$similarityType == "Lrw"){

      ioSubEnv$eigen.values <- sort(ioSubEnv$eigen.values, decreasing = F, method = "quick")
      ioSubEnv$eigen.vectors <- ioSubEnv$eigen.vectors[,rev(seq_len(ncol(ioSubEnv$eigen.vectors)))]

      if(!isTRUE(ioSubEnv$legacy.lrw)){

        # LaplacianRW decomposed Lsym, so these are Lsym's eigenvectors. Map them
        # back to Lrw's: Lrw = D^-1/2.Lsym.D^1/2, hence v_Lrw = D^-1/2.u_Lsym for
        # the same eigenvalue.
        tmp.degree <- ioSubEnv$decompose.degree
        if(is.null(tmp.degree)){
          # decompose.mat came from the cache; the degrees are recoverable from
          # the (already zero-degree-filtered) affinity matrix.
          tmp.degree <- colSums(ioSubEnv$affinity.mat, na.rm = TRUE)
        }
        ioSubEnv$eigen.vectors <- ioSubEnv$eigen.vectors / sqrt(tmp.degree)

        # An eigenvector is only defined up to sign and the solver's choice is
        # not permutation-stable; a flip swaps an eigenvector's G+ and G- subsets.
        # Pin it with a rule that does not depend on the row order: the sign of
        # the column sum, falling back to the largest-magnitude component (ties
        # broken by variable name) when that sum is ~0.
        tmp.names <- rownames(ioSubEnv$eigen.vectors)
        for(iCol in seq_len(ncol(ioSubEnv$eigen.vectors))){
          tmp.v <- ioSubEnv$eigen.vectors[, iCol]
          tmp.s <- sum(tmp.v)
          if(abs(tmp.s) > 1e-12){
            tmp.sign <- sign(tmp.s)
          } else {
            tmp.max <- max(abs(tmp.v))
            tmp.cand <- which(abs(abs(tmp.v) - tmp.max) <= 1e-15)
            tmp.pick <- tmp.cand[order(tmp.names[tmp.cand])][1]
            tmp.sign <- sign(tmp.v[tmp.pick])
          }
          if(tmp.sign < 0){ ioSubEnv$eigen.vectors[, iCol] <- -tmp.v }
        }
      }

    } else { cat(stop("# Not implemented yet...")) }

    # --!-- Time
    ioSubEnv$eigen.time <- (proc.time() - time.eigen.start)

    if(!is.null(ioSubEnv$inputData.filePath)){

      # Save the eigen vectors, values and eigen decomposition time
      saveRDS(ioSubEnv$eigen.vectors, tmp.path.eigen.vect)
      saveRDS(ioSubEnv$eigen.values, tmp.path.eigen.values)
      saveRDS(ioSubEnv$eigen.time, tmp.path.eigen.time)

      # Save also the eigenvalues as a table
      # tmp.dirPath <- gsub(basename(tmp.path), "eigenPlots", tmp.path)
      tmp.dirPath <- file.path(tmp.path, "eigenPlots")
      if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}
      tmp.out <- file.path( tmp.dirPath, paste(gsub(".txt", '', tmp.inputFileName),
                                             "eigenvalues.txt", sep = '-'))
      write.table(ioSubEnv$eigen.values, file=tmp.out, col.names = F, row.names = F, sep='\t')
    }

  } else {

    # Load the computed decompose matrix and computation time
    ioSubEnv$eigen.vectors <- readRDS(tmp.path.eigen.vect)
    ioSubEnv$eigen.values <- readRDS(tmp.path.eigen.values)
    ioSubEnv$eigen.time <- readRDS(tmp.path.eigen.time)
  }
  # Add to total time
  ioSubEnv$total.time <- ioSubEnv$total.time + ioSubEnv$eigen.time

  # If no eigen vector percent is given, compute it using the elbow-like heuristic
  if(ioSubEnv$subset.k.perc == -1){
    ioSubEnv$subset.k.perc <- scalenet:::computeEigenVectBestK(ioSubEnv = ioSubEnv)
    ioSubEnv$subset.k <- ceiling(ioSubEnv$subset.k.perc*dim(ioSubEnv$allData)[2])

    # In case we need FuzzyCmeansOrder, we set the number of clusters
    ioSubEnv$cmeans.k <- 2*ioSubEnv$subset.k
    if(ioSubEnv$cmeans.k < 2){stop("# --Eslib10-- Two cluster at least are expected!")}
  }

  # Take only subset.k eigenvectors, but make sure that this subset.k number is less than the number of vectors...
  if(ioSubEnv$subset.k <= ncol(ioSubEnv$eigen.vectors)){
    ioSubEnv$eigen.vectors <- ioSubEnv$eigen.vectors[,1:ioSubEnv$subset.k]
  }

  if(ioSubEnv$verbose){
    cat("# ----| Row norm:\n")
    print(head(apply(ioSubEnv$eigen.vectors, 1, function(x){sqrt(sum(x^2))})))
    cat("# ----| Col norm:\n")
    print(head(apply(ioSubEnv$eigen.vectors, 2, function(x){sqrt(sum(x^2))})))
  }
}
