#' plotSaveEigenVal
#'
#'
#' @param ioSubEnv a global environment variable

plotSaveEigenVal <- function( ioSubEnv ){

  if(!is.null(ioSubEnv$inputData.filePath)){
    tmp.dirPath <- gsub(basename(ioSubEnv$inputData.filePath), "scalenet",
                        ioSubEnv$inputData.filePath)
    if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}

    tmp.plot.path <- file.path( tmp.dirPath, "eigenPlots",
                                paste(gsub(".txt", '', basename(ioSubEnv$inputData.filePath)),
                                      "eigenvalues.pdf", sep = '-'))
    if(!file.exists(tmp.plot.path)){
      pdf(file = tmp.plot.path, paper = 'a4r')
      par(oma=c(0,0,2,0))
      plot(round(ioSubEnv$eigen.values, digits=5), xlab = "#eigen vector", ylab = "eigen value", xaxt = "n")
      title(paste(basename(ioSubEnv$inputData.filePath), " - Eigen values of ",
                  paste(ioSubEnv$mi.estimator, ioSubEnv$similarityType, sep = ' - '), sep=''), outer = TRUE)
      axis(1, at = seq_len(ncol(ioSubEnv$allData)), labels = seq_len(ncol(ioSubEnv$allData)), cex.axis = 0.75)
      abline(v=seq_len(ncol(ioSubEnv$allData)), col="gray", lty=3)
      dev.off()
    }
  }
}
