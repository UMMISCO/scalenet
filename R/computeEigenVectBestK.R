#' computeEigenVectBestK
#'
#' Compute the best number of eigen vectors using a two regression fit
#'
#' @param ioSubEnv a global environment variable
#'
#'

computeEigenVectBestK <- function(ioSubEnv){

  # Set a step to fit the data
  myStep <- 10

  # Get a copy of the eigen values
  eigenVal.vect <- sort(ioSubEnv$eigen.values, decreasing = FALSE)
  p <- length(eigenVal.vect)

  # Perform the first regression --> [p/2-n <---- p/2 ----> p/2+n]
  # ------
  myRSquared <- c()

  for(n in seq(from = myStep, to = p, by=myStep)){

    n.start <- (ceiling(p/2) - ceiling(n/2))
    n.end <- (ceiling(p/2) + ceiling(n/2))
    if(n.start <= 0 | n.end > p){break;}

    fit <- lm(eigenVal.vect[n.start:n.end] ~ c(n.start:n.end))
    myRSquared <- c(myRSquared, summary(fit)$r.squared)
  }
  n1 <- max(which(round(max(myRSquared)-myRSquared, digits = 4)<10^-1))*myStep
  n1.start <- (ceiling(p/2) - ceiling(n1/2))
  n1.end <- (ceiling(p/2) + ceiling(n1/2))
  fit.n1 <- lm(eigenVal.vect[n1.start:n1.end] ~ c(n1.start:n1.end))

  # Perform the second regression --> [1:n1.start]
  # ------
  fit.n2 <- lm(eigenVal.vect[1:n1.start] ~ seq_len(length(eigenVal.vect))[1:n1.start])
  best.k <- (fit.n2$coefficients[[1]]-fit.n1$coefficients[[1]])/((fit.n1$coefficients[[2]]-fit.n2$coefficients[[2]]))

  # Keep a plot of the eigen values, the two regression and the best k
  # ------
  tmp.dirPath <- file.path(ioSubEnv$output.dirPath, "eigenVectBestK")
  if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}
  tmp.filePath <- file.path(tmp.dirPath, "eigenVectBestK.pdf")
  pdf(file = tmp.filePath, width=10, height = 5)
  plot(y=eigenVal.vect, x= seq_len(p))
  abline(fit.n1, col="darkblue")
  abline(v=n1.start, col="darkblue", lty=2)
  abline(v=n1.end, col="darkblue", lty=2)
  abline(fit.n2, col="magenta")
  suppressWarnings(mtext(paste("best k:", as.character(format(best.k/p, digits = 4)), sep="\t"), col = "red"))

  # Two first fit intersection
  abline(v=best.k,col="red")
  dev.off()

  return(as.numeric(format(best.k/p, digits = 4)))
}
