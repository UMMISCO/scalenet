#' plotSaveEigenVect
#'
#'
#' @param ioSubEnv a global environment variable

plotSaveEigenVect <- function( ioSubEnv ){

  if(!is.null(ioSubEnv$inputData.filePath)){
    # Plot the sorted elements of each eigen vector
    # --> prepare a split of the eigen vectors for plot conveniency...
    # --> here, consider the Lwr matrix colnames, as some variables might have been ignored earlier due to 0 degree
    myVector.idx.list <- split(seq_len(ncol(ioSubEnv$eigen.vectors)),
                               ceiling(seq_along(seq_len(ncol(ioSubEnv$eigen.vectors)))/8))

    tmp.dirPath <- file.path(gsub(basename(ioSubEnv$inputData.filePath), "scalenet",
                                  ioSubEnv$inputData.filePath), "eigenPlots", "eigenVectors")
    if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}

    # --> plot/save the eigen vector elements
    for(iIdx in seq_len(length(myVector.idx.list))){

      tmp.plot.path <- file.path(tmp.dirPath, paste(paste("eigenVectors", min(myVector.idx.list[[iIdx]]),
                                                          "to", max(myVector.idx.list[[iIdx]]), sep = "_"),
                                                    "pdf", sep = "."))

      if(!file.exists(tmp.plot.path)){

        pdf(tmp.plot.path)

        par(mfrow = c(4,2))
        for(iVect in myVector.idx.list[[iIdx]]){

          # Get a vector (ie., a column) and sort its elements
          tmp.elts <- sort(ioSubEnv$eigen.vectors[,iVect], decreasing = TRUE, method = "quick")

          # Plot
          plot(tmp.elts, main = paste("Eigen vector", colnames(ioSubEnv$eigen.vectors)[iVect],
                                      "(signed)"),
               ylab = "Elt. values", xaxt = 'n', ylim = c(-1.0, 1.1), xlab = '')

          abline(h=0, col = "blue", lwd = 1.5)
          abline(h=10^-1, col = "blue", lwd = 1, lty=2)
          abline(h=-10^-1, col = "blue", lwd = 1, lty=2)
          abline(h=seq(-1,1,by=0.1), col="gray", lty=3)
          abline(v=seq_len(ncol(ioSubEnv$affinity.mat)), col="gray", lty=3)
          axis(side = 1, at = seq_len(length(tmp.elts)), labels = names(tmp.elts),
               las = 2, cex.axis = 0.65)
        }
        par(mfrow = c(1,1)); dev.off()
      }
    }
  }
}
