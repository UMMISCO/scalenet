#' discretize
#'
#' discretize discretizes numerical column of a data.frame. The data.frame can contain NA values.
#' Outlier values are identified using the Tukey's test. A pamk clustering is done on the non NA/no outlier values.
#' Upper outlier values are inserted in the cluster with the largest medoid, and lower outlier values with
#' the cluster of smallest medoid. If non outlier values are all equal, the values are put in one cluster, and the outlier
#' values, if any, are put in different clusters.
#'
#' @param argInData input data file path (a tab separated file) or a data.frame, with features in columns and observations in rows
#' @param argDisc, the discretizing method, "pamk", "equalfreq" or "equalwidth"; default is "pamk"
#' @param argMaxClusters, the expected maximum number of clusters (only used with the "pamk" method)
#' @param argBins, the required number of bins (used with "equalfreq" or "equalwidth"; default is 5)
#'
#' @export

# discretize(argInData = "~/Projects/Projects_largeScale/data/microbaria/input/rawData/clin/test_discretize/clin.txt")
# discretize.plot(argInData = "~/Projects/Projects_largeScale/data/microbaria/input/rawData/clin/test_discretize/clin.txt")

discretize <- function(argInData, argDisc = "pamk", argBins = 5,
                           argMaxClusters = 5, argVerbose = FALSE){
  # Read/load the data
  myData <- NULL
  if(class(argInData) == "character"){

    if(!file.exists(argInData)){stop("# --Err-- Input file does not exist")}
    cat("# Loading data:", basename(argInData), "...\n")
    myData <- read.table(argInData, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # myData <- suppressWarnings(data.table::fread(input = argInData, sep = "\t", stringsAsFactors = F,
    #                                              data.table = F, showProgress = F,
    #                                              header = TRUE))
  } else if(class(argInData) == "data.frame") {

    cat("# Loading data...\n")
    myData <- argInData

  } else { stop("# --Err-- Input should be a file path or a data.frame") }

  # Set a copy where the discrete values will be stored
  myData.disc <- myData

  # Explore data and find numerical features
  prettyData <- prettyR::describe(myData, num.desc=c("min","mean","median","max","var","sd","valid.n"))
  numeric.colNames <- names(prettyData$Numeric)

  if(length(numeric.colNames)>0){

    cat("# Input data has ", length(numeric.colNames),
        " numeric variables to be discretized\n")

    if(argDisc == "pamk"){

      # Set the max number of cluster
      myMaxNbrCluster = argMaxClusters

      # Set the outlier criteria
      myOutlier.k = 1.5 # outlier
      # myOutlier.k = 3   # extreme-outlier

      # Graphical parameters
      # iGraphCount <- 1; myCexAxis = 0.5; nbRowGraph = 10; nbColGraph = 9

      # if(class(argInData) == "character"){
      #   pdf(file = paste(argInData, "pdf", sep = "."), paper = "a4", width = 7, height = 25)
      # }

      # For each numerical feature
      counter = 0
      for(strNumVar in numeric.colNames){

        counter = counter + 1

        # Load the continuous data
        myVar <- strNumVar
        myCont.values <- myData[, myVar]
        if(argVerbose){
          cat("# ------------------------------------\n")
          cat("# --> counter:", counter, "--> var:", myVar,"\n")
          cat("# ------------------------------------\n")
          print(str(myCont.values))
          cat("# ------------------------------------\n")
        }

        # Get the index of the values that are not NA
        # ----
        myCont.values.notNa.idx <- which(!is.na(myCont.values))

        # Check if there are enough not NA values,
        # otherwise keep all NA, or set all to 1 and go next
        if(length(myCont.values.notNa.idx) < (myMaxNbrCluster+1)){

          # If all NAs, go next
          if(length(myCont.values.notNa.idx) == 0){next}

          # If only one or two non NA, set to 1 and go next
          if(length(myCont.values.notNa.idx) < 3){
            myData.disc[myCont.values.notNa.idx, strNumVar] <- 1
            next
          }
        }

        # Check if the variable is not a constant
        if(length(unique(myCont.values[myCont.values.notNa.idx])) ==  1){

          # If yes, set to 1 and go next
          myData.disc[myCont.values.notNa.idx, strNumVar] <- 1
          next
        }

        # Identify the outliers
        # ----
        # Compute IQR
        lowerq = quantile(myCont.values, na.rm = TRUE)[2]
        upperq = quantile(myCont.values, na.rm = TRUE)[4]
        iqr = upperq - lowerq #Or use IQR(data)

        # Identify 'positive' outliers
        extreme.threshold.upper = (iqr * myOutlier.k) + upperq
        myCont.values.posOutlier.idx <- which(myCont.values > extreme.threshold.upper)
        # Identify 'negative' outliers
        extreme.threshold.lower = lowerq - (iqr * myOutlier.k)
        myCont.values.negOutlier.idx <- which(myCont.values < extreme.threshold.lower)
        # Gather all outliers index
        myCont.values.outlier.idx <- unique(c(myCont.values.posOutlier.idx, myCont.values.negOutlier.idx))

        # Get the non NA/outliers values
        myCont.values.noNAOutliers <- myCont.values
        tmp.idx <- c(myCont.values.outlier.idx, which(!seq_len(length(myCont.values)) %in% myCont.values.notNa.idx))
        if(length(tmp.idx)>0){ myCont.values.noNAOutliers <- myCont.values[-tmp.idx] }

        # perform clustering
        # ----
        # Init the clustering and medoid vectors
        currClustering.vect <- c()
        currMedoid.vect <- c()

        # If all values are identical, set 1 cluster and one medoid
        if(length(unique(myCont.values.noNAOutliers)) == 1){

          currClustering.vect <- rep(1, length(myCont.values.noNAOutliers))
          currMedoid.vect <- mean(myCont.values.noNAOutliers)

        } else{ # else do pamk

          # Clustering the non outliers
          # ----
          tmp.maxNbrCluster <- myMaxNbrCluster
          # Set the range of cluster number to explore
          myKrange <- c(2:tmp.maxNbrCluster)
          # If number of non NA <= 5, set max n-1 clusters
          if(length(myCont.values.noNAOutliers) <= max(myKrange)){
            myKrange <- c(2:(length(myCont.values.noNAOutliers)-1))
            tmp.maxNbrCluster <- (length(myCont.values.noNAOutliers)-1)
          }

          # Perform a medoid clustering on the not NA/outlier values
          if(argVerbose){
            cat("# --> pamk, range :", paste(myKrange, collapse = ','), "\n# ------------------------------------\n")
          }

          # --> the maximum number of clusters is set to max myKrange
          # --!!--| if the number of none NAs is lower than max myKrange, there is an error...
          # --!!--| thus, set either max myKrange to the number of none NAs values...
          pam.out <- pamk(myCont.values.noNAOutliers, krange = myKrange,
                          criterion="asw", metric = "manhattan", diss = FALSE, critout=argVerbose)

          # Get the clustering and medoid vectors
          currClustering.vect <- as.vector(pam.out$pamobject$clustering)
          currMedoid.vect <- as.vector(pam.out$pamobject$medoids)

        }
        # Keep the number of clusters for the non NA/outlier values
        currMedoid.noOutliers.length <- length(currMedoid.vect)

        # Reshape the clustering vector to insert the upper and lower outliers
        currClustering.noNaShape.vect <- rep(NA, length(myCont.values))

        # -> insert not NA and non outliers
        tmp.idx <- sort(c(myCont.values.outlier.idx,
                          which(!seq_len(length(myCont.values)) %in% myCont.values.notNa.idx)),
                        decreasing = FALSE)
        if(length(tmp.idx)>0){ currClustering.noNaShape.vect[-tmp.idx] <- currClustering.vect
        } else { currClustering.noNaShape.vect <- currClustering.vect }

        # -> insert upper outliers
        if(length(myCont.values.posOutlier.idx)>0){

          if(currMedoid.noOutliers.length > 1){

            maxMedoid.idx <- which(currMedoid.vect == max(currMedoid.vect))
            currClustering.noNaShape.vect[myCont.values.posOutlier.idx] <- maxMedoid.idx

          } else {

            currClustering.noNaShape.vect[myCont.values.posOutlier.idx] <- 2
            currMedoid.vect <- c(currMedoid.vect, mean(myCont.values[myCont.values.posOutlier.idx]))
          }
        }

        # -> insert lower outliers
        if(length(myCont.values.negOutlier.idx)>0){

          if(currMedoid.noOutliers.length > 1){

            minMedoid.idx <- which(currMedoid.vect == min(currMedoid.vect))
            currClustering.noNaShape.vect[myCont.values.negOutlier.idx] <- minMedoid.idx

          } else {

            currClustering.noNaShape.vect[myCont.values.negOutlier.idx] <- 3
            currMedoid.vect <- c(currMedoid.vect, mean(myCont.values[myCont.values.negOutlier.idx]))
          }
        }

        # Finally, rm NA
        currClustering.noNaShape.vect <- as.vector(na.omit(currClustering.noNaShape.vect))

        # Merge small cluster with most similar cluster
        # ----

        # Define smallest cluster size
        small.size <- ifelse( currMedoid.noOutliers.length == 1,
                              ceiling(0.75*length(currClustering.noNaShape.vect)/3),
                              ceiling(0.75*length(currClustering.noNaShape.vect)/max(currClustering.noNaShape.vect)))

        # Check if some clusters are smaller than small size
        clust.small <- which(table(currClustering.noNaShape.vect)<small.size)

        # While there are too small clusters, do merge (until 2 clusters left)
        while(length(clust.small) > 0 & max(currClustering.noNaShape.vect) > 2){

          # For all clust pairs, compute medoid difference
          allPairs <- combn(seq_len(max(currClustering.noNaShape.vect)), 2)
          allPairs.medDiff <- apply(allPairs, 2, function(currPair){ abs(currMedoid.vect[currPair[1]]-currMedoid.vect[currPair[2]]) })
          medDiff.mat <- matrix(NA, ncol=3, nrow=ncol(allPairs))
          medDiff.mat[,1] <- allPairs[1,]
          medDiff.mat[,2] <- allPairs[2,]
          medDiff.mat[,c(3)] <- allPairs.medDiff
          medDiff.mat <- medDiff.mat[order(medDiff.mat[,3], decreasing = FALSE),]

          # Find the first pair with smallest difference and one cluster with too small size
          tmp.bestMerge.idx <- sort(which(medDiff.mat[,1] %in% clust.small | medDiff.mat[,2] %in% clust.small), decreasing = FALSE)[1]
          tmp.bestMerge <- as.vector(medDiff.mat[tmp.bestMerge.idx, c(1,2)])

          # As two clusters are merged, cluster numbers should be reorganized
          currClustering.conversion.vect <- numeric(length = max(currClustering.noNaShape.vect))
          currClustering.conversion.vect[tmp.bestMerge] <- (max(currClustering.noNaShape.vect)-1)
          currClustering.conversion.vect[currClustering.conversion.vect!=(max(currClustering.noNaShape.vect)-1)] <- seq_len(max(currClustering.noNaShape.vect)-2)
          for(iIdx in seq_len(length(currClustering.noNaShape.vect))){
            currClustering.noNaShape.vect[iIdx] <- currClustering.conversion.vect[currClustering.noNaShape.vect[iIdx]]
          }

          # As two clusters are merged, define a new medoid vector
          newMedoid.vect <- numeric(length = (max(currClustering.noNaShape.vect)))
          newMedoid.vect[max(currClustering.noNaShape.vect)] <- mean(currMedoid.vect[tmp.bestMerge])
          for(iIdx in seq_len(max(currClustering.noNaShape.vect)-1)){
            newMedoid.vect[iIdx] <- currMedoid.vect[which(currClustering.conversion.vect == iIdx)]
          }
          currMedoid.vect <- newMedoid.vect

          # Check if some clusters are still smaller than small size
          clust.small <- which(table(currClustering.noNaShape.vect)<small.size)
        }

        # TEST TEST
        # Make two plots:
        # --> left: unsorted myCont.values colored with the cluster number
        # --> middle: hist on the number of clusters
        # --> right: Sturges hist

        # Set a vector with the cluster number or NA of the continuous value was NA
        myCont.clust <- myCont.values
        # print(myCont.values)
        # myCont.clust[myCont.values.notNa.idx] <- pam.out$pamobject$clustering
        myCont.clust[myCont.values.notNa.idx] <- currClustering.noNaShape.vect

        # --!!-- order the cluster number by medoid values
        tmp.clustering <- myCont.clust
        # tmp.medoids.order <- order(pam.out$pamobject$medoids, decreasing = FALSE)
        tmp.medoids.order <- order(currMedoid.vect, decreasing = FALSE)
        # print(tmp.medoids.order)
        for(iMed in seq_len(length(tmp.medoids.order))){
          myCont.clust[tmp.clustering == tmp.medoids.order[iMed]] <- iMed
        }
        # print(myCont.clust)
        # plot(x=myCont.clust, y=myCont.values)

        # Copy into discretize dataset
        myData.disc[myCont.values.notNa.idx, strNumVar] <- myCont.clust[myCont.values.notNa.idx]

        # # Set a vector with a color for each group, and black for NA
        # myColor.clust <- rep("black", length(myCont.values))
        # myColors.ref <- c("blue", "red", "magenta", "orange", "green")
        # myColor.clust <- myColors.ref[myCont.clust]
        # myColor.clust[is.na(myColor.clust)] <- "black"

        # # Bind continuous, cluster and color
        # myCont.idx <- seq_len(length(myCont.values))
        # myCont.clust.color <- cbind.data.frame(myCont.values, myCont.idx, myCont.clust, myColor.clust)
        # myCont.clust.color; colnames(myCont.clust.color) <- c("values", "index", "cluster", "color")
        # myCont.notNA.idx <- which(!is.na(myCont.clust.color$values))

        # ----
        #   if(iGraphCount == 1){
        #     par(mar=c(2,2,2,2), mfrow=c(nbRowGraph, nbColGraph))
        #     iGraphCount <- (iGraphCount+1)
        #   }

        # if(iGraphCount == 1){
        #   layout(matrix(1:(nbRowGraph*nbColGraph), nbRowGraph, nbColGraph, byrow = TRUE), respect = TRUE)
        #   par(omi=c(0,0,0,0), mar=c(1.5, 1.5, 1.5, 1.5))
        #   iGraphCount <- (iGraphCount+1)
        # }
        # # ----

        # nbrClusters <- max(myCont.clust.color$cluster, na.rm = T)

        # # --> left: Sturges hist
        # hist(myCont.clust.color$values[myCont.notNA.idx],
        #      # main = paste("Sturges classes"), cex.main = 0.7,
        #      xlab = "value",
        #      xlim = c(min(myCont.values, na.rm = T), max(myCont.values, na.rm = T)),
        #      ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
        # mtext(side = 3, text = paste("Sturges"), line = 0.3, cex = myCexAxis)

        # # --> middle: unsorted myCont.values colored with the cluster number
        # myCont.clust.color.orderByValues <- myCont.clust.color[order(myCont.clust.color$values, decreasing = FALSE),]
        # # plot(y = myCont.clust.color$values[myCont.notNA.idx],
        # #      x = myCont.clust.color$index[myCont.notNA.idx],
        # #      col = myCont.clust.color$color[myCont.notNA.idx],
        # #      main = paste(nbrClusters, "(pamk) clusters"),
        # #      xlab = "sample", ylab = "value")
        # plot(y = myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)],
        #      x = seq_len(length(which(!is.na(myCont.clust.color.orderByValues$values)))),
        #      col = as.vector(myCont.clust.color.orderByValues$color[!is.na(myCont.clust.color.orderByValues$values)]),
        #      xlab = "sample", ylab = "",
        #      ylim = c(min(myCont.values, na.rm = T), max(myCont.values, na.rm = T)),
        #      ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
        # mtext(side = 3, text = myVar, line = 0.2, cex = myCexAxis)
        #
        # Get the original values 25%, 50% and 75 % values
        # myQuantile <- quantile(myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)])
        # if(argVerbose){
        #   cat("# --> quantile :", paste(myQuantile[c(2:4)], collapse = '  |  '), "\n# ------------------------------------\n")
        # }
        # abline(h = myQuantile[2], lty = 1, lwd = 2, col = "gray")
        # abline(h = myQuantile[3], lty = 1, lwd = 2, col = "gray")
        # abline(h = myQuantile[4], lty = 1, lwd = 2, col = "gray")
        # abline(h = mean(myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)]),
        #        lty = 2, lwd = 2, col = "black")
        #
        # # --> right: hist on the pamk clusters
        # myCont.maxByClust.agg <- aggregate(values ~ cluster, data = myCont.clust.color, max)
        # myCont.minByClust.agg <- aggregate(values ~ cluster, data = myCont.clust.color, min)
        #
        # if(argVerbose){
        #   cat("# --> min value per cluster :\n")
        #   print(myCont.minByClust.agg)
        #   cat("# --> max value per cluster :\n")
        #   print(myCont.maxByClust.agg)
        #   cat("\n# ------------------------------------\n")
        # }

        # hist(myCont.clust.color$values[myCont.notNA.idx], breaks = unique(c( min(myCont.minByClust.agg$values),
        #                                                                      sort(myCont.maxByClust.agg$values,
        #                                                                           decreasing = F))),
        #      xlab = "value",
        #      xlim = c(min(myCont.values, na.rm = T), max(myCont.values, na.rm = T)),
        #      ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
        # mtext(side = 3, text = paste(nbrClusters, "ext-pamk"), line = 0.2, cex = myCexAxis)
        #
        # if(iGraphCount == (nbRowGraph*nbColGraph)){
        #   # par(mfrow=c(1,1))
        #
        #   iGraphCount = 1
        # }
        # # ----
      }

      # if(class(argInData) == "character"){dev.off()}

      # -- Make sure all columns are set to factor
      myData.disc[, colnames(myData.disc)] <- lapply(myData.disc[, colnames(myData.disc)], as.factor)

      if(class(argInData) == "character"){
        # - with rownames
        write.table(myData.disc, file=paste(argInData, "disc.txt", sep = "_"),
                    col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
        # - without rownames
        write.table(myData.disc, file=paste(argInData, "disc_noRowNames.txt", sep = "_"),
                    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
      }
    } else if( argDisc %in% c("equalwidth", "equalfreq")){

      # Order the cluster number by the mean of the corresponding continuous values
      for(strCol in numeric.colNames){

        notNA.idx <- which(!is.na(myData[, strCol]))
        continuousValues <- as.vector(as.matrix(myData[notNA.idx, strCol]))
        discreteValues <- as.vector(as.matrix(infotheo::discretize(continuousValues,
                                                                   disc = argDisc, nbins = argBins)))
        # For each cluster, compute the mean over the corresponding continuous values
        possibleClusterNbr <- unique(discreteValues)
        clusterMean <- rep(NA, length(possibleClusterNbr))
        names(clusterMean) <- as.character(possibleClusterNbr)
        for(iClust in possibleClusterNbr){
          clusterMean[as.character(iClust)] <- mean(continuousValues[which(discreteValues == iClust)])
        }

        # Reorder the cluster number based on the order of the cluster mean values
        clusterMean.ordered <- clusterMean[order(clusterMean, decreasing = F)]
        tmp.discreteValues <- discreteValues

        for(iClust in seq_len(length(clusterMean.ordered))){
          discreteValues[which(tmp.discreteValues == as.numeric(names(clusterMean.ordered)[iClust]))] <- iClust
        }

        # Set the disc.data data.frame with the cluster number
        myData.disc[notNA.idx, strCol] <- discreteValues

      }
    } else {stop("Unknown discretize method")}

  } else {
    cat("# Nothing to be discretized...\n")
  }

  return(list(discData = myData.disc, numColName = numeric.colNames))
}

# discretize.bck <- function(argInData, argDisc = "pamk", argBins = 5,
#                        argMaxClusters = 5, argVerbose = FALSE){
#
#   # Read/load the data
#   myData <- NULL
#   if(class(argInData) == "character"){
#
#     if(!file.exists(argInData)){stop("# --Err-- Input file does not exist")}
#     cat("# Loading data:", basename(argInData), "...\n")
#     myData <- read.table(argInData, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#     # myData <- suppressWarnings(data.table::fread(input = argInData, sep = "\t", stringsAsFactors = F,
#     #                                              data.table = F, showProgress = F,
#     #                                              header = TRUE))
#   } else if(class(argInData) == "data.frame") {
#
#     cat("# Loading data...\n")
#     myData <- argInData
#
#   } else { stop("# --Err-- Input should be a file path or a data.frame") }
#
#   # Set a copy where the discrete values will be stored
#   myData.disc <- myData
#
#   # Explore data and find numerical features
#   prettyData <- prettyR::describe(myData, num.desc=c("min","mean","median","max","var","sd","valid.n"))
#   numeric.colNames <- names(prettyData$Numeric)
#
#   if(length(numeric.colNames)>0){
#
#     cat("# Input data has ", length(numeric.colNames), " numeric variables to be discretized\n")
#
#     if(argDisc == "pamk"){
#
#       # Set the max number of cluster
#       myMaxNbrCluster = argMaxClusters
#
#       # Set the outlier criteria
#       myOutlier.k = 1.5 # outlier
#       # myOutlier.k = 3   # extreme-outlier
#
#       # Graphical parameters
#       iGraphCount <- 1; myCexAxis = 0.5; nbRowGraph = 10; nbColGraph = 9
#
#       if(class(argInData) == "character"){
#         pdf(file = paste(argInData, "pdf", sep = "."), paper = "a4", width = 7, height = 25)
#       }
#
#       # For each numerical feature
#       for(strNumVar in numeric.colNames){
#
#         # Load the continuous data
#         myVar <- strNumVar
#         myCont.values <- myData[, myVar]
#         if(argVerbose){
#           cat("# ------------------------------------\n# --> var:", myVar, "\n# ------------------------------------\n")
#           print(myCont.values)
#           cat("# ------------------------------------\n")
#         }
#
#         # Get the index of the values that are not NA
#         # ----
#         myCont.values.notNa.idx <- which(!is.na(myCont.values))
#
#         # Check if there are enough not NA values,
#         # otherwise keep all NA, or set all to 1 and go next
#         if(length(myCont.values.notNa.idx) < (myMaxNbrCluster+1)){
#
#           # If all NAs, go next
#           if(length(myCont.values.notNa.idx) == 0){next}
#
#           # If only one or two non NA, set to 1 and go next
#           if(length(myCont.values.notNa.idx) < 3){
#             myData.disc[myCont.values.notNa.idx, strNumVar] <- 1
#             next
#           }
#         }
#
#         # Check if the variable is not a constant
#         if(length(unique(myCont.values[myCont.values.notNa.idx])) ==  1){
#
#           # If yes, set to 1 and go next
#           myData.disc[myCont.values.notNa.idx, strNumVar] <- 1
#           next
#         }
#
#
#         # Identify the outliers
#         # ----
#         # Compute IQR
#         lowerq = quantile(myCont.values, na.rm = TRUE)[2]
#         upperq = quantile(myCont.values, na.rm = TRUE)[4]
#         iqr = upperq - lowerq #Or use IQR(data)
#
#         # Identify 'positive' outliers
#         extreme.threshold.upper = (iqr * myOutlier.k) + upperq
#         myCont.values.posOutlier.idx <- which(myCont.values > extreme.threshold.upper)
#         # Identify 'negative' outliers
#         extreme.threshold.lower = lowerq - (iqr * myOutlier.k)
#         myCont.values.negOutlier.idx <- which(myCont.values < extreme.threshold.lower)
#         # Gather all outliers index
#         myCont.values.outlier.idx <- unique(c(myCont.values.posOutlier.idx, myCont.values.negOutlier.idx))
#
#         # Get the non NA/outliers values
#         myCont.values.noNAOutliers <- myCont.values
#         tmp.idx <- c(myCont.values.outlier.idx, which(!seq_len(length(myCont.values)) %in% myCont.values.notNa.idx))
#         if(length(tmp.idx)>0){ myCont.values.noNAOutliers <- myCont.values[-tmp.idx] }
#
#         # perform clustering
#         # ----
#         # Init the clustering and medoid vectors
#         currClustering.vect <- c()
#         currMedoid.vect <- c()
#
#         # If all values are identical, set 1 cluster and one medoid
#         if(length(unique(myCont.values.noNAOutliers)) == 1){
#
#           currClustering.vect <- rep(1, length(myCont.values.noNAOutliers))
#           currMedoid.vect <- mean(myCont.values.noNAOutliers)
#
#         } else{ # else do pamk
#
#           # Clustering the non outliers
#           # ----
#           tmp.maxNbrCluster <- myMaxNbrCluster
#           # Set the range of cluster number to explore
#           myKrange <- c(2:tmp.maxNbrCluster)
#           # If number of non NA <= 5, set max n-1 clusters
#           if(length(myCont.values.noNAOutliers) <= max(myKrange)){
#             myKrange <- c(2:(length(myCont.values.noNAOutliers)-1))
#             tmp.maxNbrCluster <- (length(myCont.values.noNAOutliers)-1)
#           }
#
#           # Perform a medoid clustering on the not NA/outlier values
#           if(argVerbose){
#             cat("# --> pamk, range :", paste(myKrange, collapse = ','), "\n# ------------------------------------\n")
#           }
#
#           # --> the maximum number of clusters is set to max myKrange
#           # --!!--| if the number of none NAs is lower than max myKrange, there is an error...
#           # --!!--| thus, set either max myKrange to the number of none NAs values...
#           pam.out <- pamk(myCont.values.noNAOutliers, krange = myKrange,
#                           criterion="asw", metric = "manhattan", diss = FALSE, critout=argVerbose)
#
#           # Get the clustering and medoid vectors
#           currClustering.vect <- as.vector(pam.out$pamobject$clustering)
#           currMedoid.vect <- as.vector(pam.out$pamobject$medoids)
#
#         }
#         # Keep the number of clusters for the non NA/outlier values
#         currMedoid.noOutliers.length <- length(currMedoid.vect)
#
#         # Reshape the clustering vector to insert the upper and lower outliers
#         currClustering.noNaShape.vect <- rep(NA, length(myCont.values))
#
#         # -> insert not NA and non outliers
#         tmp.idx <- sort(c(myCont.values.outlier.idx,
#                           which(!seq_len(length(myCont.values)) %in% myCont.values.notNa.idx)),
#                         decreasing = FALSE)
#         if(length(tmp.idx)>0){ currClustering.noNaShape.vect[-tmp.idx] <- currClustering.vect
#         } else { currClustering.noNaShape.vect <- currClustering.vect }
#
#         # -> insert upper outliers
#         if(length(myCont.values.posOutlier.idx)>0){
#
#           if(currMedoid.noOutliers.length > 1){
#
#             maxMedoid.idx <- which(currMedoid.vect == max(currMedoid.vect))
#             currClustering.noNaShape.vect[myCont.values.posOutlier.idx] <- maxMedoid.idx
#
#           } else {
#
#             currClustering.noNaShape.vect[myCont.values.posOutlier.idx] <- 2
#             currMedoid.vect <- c(currMedoid.vect, mean(myCont.values[myCont.values.posOutlier.idx]))
#           }
#         }
#
#         # -> insert lower outliers
#         if(length(myCont.values.negOutlier.idx)>0){
#
#           if(currMedoid.noOutliers.length > 1){
#
#             minMedoid.idx <- which(currMedoid.vect == min(currMedoid.vect))
#             currClustering.noNaShape.vect[myCont.values.negOutlier.idx] <- minMedoid.idx
#
#           } else {
#
#             currClustering.noNaShape.vect[myCont.values.negOutlier.idx] <- 3
#             currMedoid.vect <- c(currMedoid.vect, mean(myCont.values[myCont.values.negOutlier.idx]))
#           }
#         }
#
#         # Finally, rm NA
#         currClustering.noNaShape.vect <- as.vector(na.omit(currClustering.noNaShape.vect))
#
#         # Merge small cluster with most similar cluster
#         # ----
#
#         # Define smallest cluster size
#         small.size <- ifelse( currMedoid.noOutliers.length == 1,
#                               ceiling(0.75*length(currClustering.noNaShape.vect)/3),
#                               ceiling(0.75*length(currClustering.noNaShape.vect)/max(currClustering.noNaShape.vect)))
#
#         # Check if some clusters are smaller than small size
#         clust.small <- which(table(currClustering.noNaShape.vect)<small.size)
#
#         # While there are too small clusters, do merge (until 2 clusters left)
#         while(length(clust.small) > 0 & max(currClustering.noNaShape.vect) > 2){
#
#           # For all clust pairs, compute medoid difference
#           allPairs <- combn(seq_len(max(currClustering.noNaShape.vect)), 2)
#           allPairs.medDiff <- apply(allPairs, 2, function(currPair){ abs(currMedoid.vect[currPair[1]]-currMedoid.vect[currPair[2]]) })
#           medDiff.mat <- matrix(NA, ncol=3, nrow=ncol(allPairs))
#           medDiff.mat[,1] <- allPairs[1,]
#           medDiff.mat[,2] <- allPairs[2,]
#           medDiff.mat[,c(3)] <- allPairs.medDiff
#           medDiff.mat <- medDiff.mat[order(medDiff.mat[,3], decreasing = FALSE),]
#
#           # Find the first pair with smallest difference and one cluster with too small size
#           tmp.bestMerge.idx <- sort(which(medDiff.mat[,1] %in% clust.small | medDiff.mat[,2] %in% clust.small), decreasing = FALSE)[1]
#           tmp.bestMerge <- as.vector(medDiff.mat[tmp.bestMerge.idx, c(1,2)])
#
#           # As two clusters are merged, cluster numbers should be reorganized
#           currClustering.conversion.vect <- numeric(length = max(currClustering.noNaShape.vect))
#           currClustering.conversion.vect[tmp.bestMerge] <- (max(currClustering.noNaShape.vect)-1)
#           currClustering.conversion.vect[currClustering.conversion.vect!=(max(currClustering.noNaShape.vect)-1)] <- seq_len(max(currClustering.noNaShape.vect)-2)
#           for(iIdx in seq_len(length(currClustering.noNaShape.vect))){
#             currClustering.noNaShape.vect[iIdx] <- currClustering.conversion.vect[currClustering.noNaShape.vect[iIdx]]
#           }
#
#           # As two clusters are merged, define a new medoid vector
#           newMedoid.vect <- numeric(length = (max(currClustering.noNaShape.vect)))
#           newMedoid.vect[max(currClustering.noNaShape.vect)] <- mean(currMedoid.vect[tmp.bestMerge])
#           for(iIdx in seq_len(max(currClustering.noNaShape.vect)-1)){
#             newMedoid.vect[iIdx] <- currMedoid.vect[which(currClustering.conversion.vect == iIdx)]
#           }
#           currMedoid.vect <- newMedoid.vect
#
#           # Check if some clusters are still smaller than small size
#           clust.small <- which(table(currClustering.noNaShape.vect)<small.size)
#         }
#
#         # TEST TEST
#         # Make two plots:
#         # --> left: unsorted myCont.values colored with the cluster number
#         # --> middle: hist on the number of clusters
#         # --> right: Sturges hist
#
#         # Set a vector with the cluster number or NA of the continuous value was NA
#         myCont.clust <- myCont.values
#         # print(myCont.values)
#         # myCont.clust[myCont.values.notNa.idx] <- pam.out$pamobject$clustering
#         myCont.clust[myCont.values.notNa.idx] <- currClustering.noNaShape.vect
#
#         # --!!-- order the cluster number by medoid values
#         tmp.clustering <- myCont.clust
#         # tmp.medoids.order <- order(pam.out$pamobject$medoids, decreasing = FALSE)
#         tmp.medoids.order <- order(currMedoid.vect, decreasing = FALSE)
#         # print(tmp.medoids.order)
#         for(iMed in seq_len(length(tmp.medoids.order))){
#           myCont.clust[tmp.clustering == tmp.medoids.order[iMed]] <- iMed
#         }
#         # print(myCont.clust)
#         # plot(x=myCont.clust, y=myCont.values)
#
#         # Copy into discretize dataset
#         myData.disc[myCont.values.notNa.idx, strNumVar] <- myCont.clust[myCont.values.notNa.idx]
#
#         # Set a vector with a color for each group, and black for NA
#         myColor.clust <- rep("black", length(myCont.values))
#         myColors.ref <- c("blue", "red", "magenta", "orange", "green")
#         myColor.clust <- myColors.ref[myCont.clust]
#         myColor.clust[is.na(myColor.clust)] <- "black"
#
#         # Bind continuous, cluster and color
#         myCont.idx <- seq_len(length(myCont.values))
#         myCont.clust.color <- cbind.data.frame(myCont.values, myCont.idx, myCont.clust, myColor.clust)
#         myCont.clust.color; colnames(myCont.clust.color) <- c("values", "index", "cluster", "color")
#         myCont.notNA.idx <- which(!is.na(myCont.clust.color$values))
#
#         # ----
#         #   if(iGraphCount == 1){
#         #     par(mar=c(2,2,2,2), mfrow=c(nbRowGraph, nbColGraph))
#         #     iGraphCount <- (iGraphCount+1)
#         #   }
#         #
#         if(iGraphCount == 1){
#           layout(matrix(1:(nbRowGraph*nbColGraph), nbRowGraph, nbColGraph, byrow = TRUE), respect = TRUE)
#           par(omi=c(0,0,0,0), mar=c(1.5, 1.5, 1.5, 1.5))
#           iGraphCount <- (iGraphCount+1)
#         }
#         # ----
#
#         nbrClusters <- max(myCont.clust.color$cluster, na.rm = T)
#
#         # --> left: Sturges hist
#         hist(myCont.clust.color$values[myCont.notNA.idx],
#              # main = paste("Sturges classes"), cex.main = 0.7,
#              xlab = "value",
#              xlim = c(min(myCont.values, na.rm = T), max(myCont.values, na.rm = T)),
#              ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
#         mtext(side = 3, text = paste("Sturges"), line = 0.3, cex = myCexAxis)
#
#         # --> middle: unsorted myCont.values colored with the cluster number
#         myCont.clust.color.orderByValues <- myCont.clust.color[order(myCont.clust.color$values, decreasing = FALSE),]
#         # plot(y = myCont.clust.color$values[myCont.notNA.idx],
#         #      x = myCont.clust.color$index[myCont.notNA.idx],
#         #      col = myCont.clust.color$color[myCont.notNA.idx],
#         #      main = paste(nbrClusters, "(pamk) clusters"),
#         #      xlab = "sample", ylab = "value")
#         plot(y = myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)],
#              x = seq_len(length(which(!is.na(myCont.clust.color.orderByValues$values)))),
#              col = as.vector(myCont.clust.color.orderByValues$color[!is.na(myCont.clust.color.orderByValues$values)]),
#              xlab = "sample", ylab = "",
#              ylim = c(min(myCont.values, na.rm = T), max(myCont.values, na.rm = T)),
#              ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
#         mtext(side = 3, text = myVar, line = 0.2, cex = myCexAxis)
#
#         # Get the original values 25%, 50% and 75 % values
#         myQuantile <- quantile(myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)])
#         if(argVerbose){
#           cat("# --> quantile :", paste(myQuantile[c(2:4)], collapse = '  |  '), "\n# ------------------------------------\n")
#         }
#         abline(h = myQuantile[2], lty = 1, lwd = 2, col = "gray")
#         abline(h = myQuantile[3], lty = 1, lwd = 2, col = "gray")
#         abline(h = myQuantile[4], lty = 1, lwd = 2, col = "gray")
#         abline(h = mean(myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)]),
#                lty = 2, lwd = 2, col = "black")
#
#
#         # --> right: hist on the pamk clusters
#         myCont.maxByClust.agg <- aggregate(values ~ cluster, data = myCont.clust.color, max)
#         myCont.minByClust.agg <- aggregate(values ~ cluster, data = myCont.clust.color, min)
#
#         if(argVerbose){
#           cat("# --> min value per cluster :\n")
#           print(myCont.minByClust.agg)
#           cat("# --> max value per cluster :\n")
#           print(myCont.maxByClust.agg)
#           cat("\n# ------------------------------------\n")
#         }
#
#         hist(myCont.clust.color$values[myCont.notNA.idx], breaks = unique(c( min(myCont.minByClust.agg$values),
#                                                                              sort(myCont.maxByClust.agg$values, decreasing = F))),
#              xlab = "value",
#              xlim = c(min(myCont.values, na.rm = T), max(myCont.values, na.rm = T)),
#              ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
#         mtext(side = 3, text = paste(nbrClusters, "ext-pamk"), line = 0.2, cex = myCexAxis)
#
#         if(iGraphCount == (nbRowGraph*nbColGraph)){
#           # par(mfrow=c(1,1))
#
#           iGraphCount = 1
#         }
#         # ----
#       }
#
#       if(class(argInData) == "character"){dev.off()}
#
#       # -- Make sure all columns are set to factor
#       myData.disc[, colnames(myData.disc)] <- lapply(myData.disc[, colnames(myData.disc)], as.factor)
#
#       if(class(argInData) == "character"){
#         # - with rownames
#         write.table(myData.disc, file=paste(argInData, "disc.txt", sep = "_"),
#                     col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
#         # - without rownames
#         write.table(myData.disc, file=paste(argInData, "disc_noRowNames.txt", sep = "_"),
#                     col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
#       }
#     } else if( argDisc %in% c("equalwidth", "equalfreq")){
#
#       myData.disc[, numeric.colNames] <- infotheo::discretize(myData.disc[, numeric.colNames], disc = argDisc,
#                                                               nbins = argBins)
#     } else {
#
#
#       stop("Unknown discretize method")
#     }
#
#   } else {
#     cat("# Nothing to be discretized...\n")
#   }
#
#   return(myData.disc)
# }





