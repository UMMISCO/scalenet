#' discretize.plot

#' @export
#'
# Examples:
# retList = discretize(argInData = "~/Projects/Projects_largeScale/data/microbaria/input/rawData/clin/test_discretize/clin.txt")
# discretize.plot(argInData = "~/Projects/Projects_largeScale/data/microbaria/input/rawData/clin/test_discretize/clin.txt", retList$discData, retList$numCol, F)

discretize.plot <- function(argInData, disc.data, numeric.colNames, argVerbose = FALSE){

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

  # Graphical parameters
  iGraphCount <- 1; myCexAxis = 0.5; nbRowGraph = 10; nbColGraph = 9

  if(class(argInData) == "character"){
    pdf(file = paste(argInData, "pdf", sep = "."), paper = "a4", width = 7, height = 25)
  }

  # For each numerical feature
  for(strNumVar in numeric.colNames){

    # Make three plots:
    # --> left: unsorted continuous values colored with the cluster number
    # --> middle: hist on the number of clusters
    # --> right: Sturges hist

    # Set a vector with a color for each group, and black for NA
    myColor.clust <- rep("black", nrow(disc.data))
    myColors.ref <- c("blue", "red", "magenta", "orange", "green")
    myColor.clust <- myColors.ref[disc.data[, strNumVar]]
    myColor.clust[is.na(myColor.clust)] <- "black"

    # print("------------------------")
    # cat(length(which(is.na(myData[, strNumVar]))), ":\t", myData[, strNumVar], "\n")
    # cat(length(which(is.na(disc.data[, strNumVar]))), ":\t", disc.data[, strNumVar], "\n")
    # # cat(length(which(is.na(myColor.clust))), ":\t", myColor.clust, "\n")
    # print("------------------------")

    # Bind continuous, cluster and color
    myCont.idx <- seq_len(nrow(disc.data))
    myCont.clust.color <- cbind.data.frame(myData[, strNumVar], myCont.idx,
                                             disc.data[, strNumVar], myColor.clust)
    colnames(myCont.clust.color) <- c("values", "index", "cluster", "color")
    myCont.notNA.idx <- which(!is.na(myCont.clust.color$values))

    if(iGraphCount == 1){
      layout(matrix(1:(nbRowGraph*nbColGraph), nbRowGraph, nbColGraph, byrow = TRUE), respect = TRUE)
      par(omi=c(0,0,0,0), mar=c(1.5, 1.5, 1.5, 1.5))
      iGraphCount <- (iGraphCount+1)
    }

    nbrClusters <- max(as.numeric(myCont.clust.color$cluster), na.rm = T)

    # --> left: Sturges hist
    hist(myCont.clust.color[myCont.notNA.idx, "values"],
         # main = paste("Sturges classes"), cex.main = 0.7,
         xlab = "value",
         xlim = c(min(as.numeric(myData[, strNumVar]), na.rm = T),
                  max(as.numeric(myData[, strNumVar]), na.rm = T)),
         ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
    mtext(side = 3, text = paste("Sturges"), line = 0.3, cex = myCexAxis)

    # --> middle: unsorted continuous values colored with the cluster number
    myCont.clust.color.orderByValues <- myCont.clust.color[order(myCont.clust.color$values,
                                                                 decreasing = FALSE),]
    plot(y = myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)],
         x = seq_len(length(which(!is.na(myCont.clust.color.orderByValues$values)))),
         col = as.vector(myCont.clust.color.orderByValues$color[!is.na(myCont.clust.color.orderByValues$values)]),
         xlab = "sample", ylab = "",
         ylim = c(min(as.numeric(myCont.clust.color[myCont.notNA.idx, "values"])),
                      max(as.numeric(myCont.clust.color[myCont.notNA.idx, "values"]))),
         ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
    mtext(side = 3, text = strNumVar, line = 0.2, cex = myCexAxis)

    # Get the original values 25%, 50% and 75 % values
    myQuantile <- quantile(myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)])
    if(argVerbose){
      cat("# --> quantile :", paste(myQuantile[c(2:4)], collapse = '  |  '), "\n# ------------------------------------\n")
    }
    abline(h = myQuantile[2], lty = 1, lwd = 2, col = "gray")
    abline(h = myQuantile[3], lty = 1, lwd = 2, col = "gray")
    abline(h = myQuantile[4], lty = 1, lwd = 2, col = "gray")
    abline(h = mean(myCont.clust.color.orderByValues$values[!is.na(myCont.clust.color.orderByValues$values)]),
           lty = 2, lwd = 2, col = "black")

    # --> right: hist on the pamk clusters
    myCont.maxByClust.agg <- aggregate(values ~ cluster, data = myCont.clust.color, max)
    myCont.minByClust.agg <- aggregate(values ~ cluster, data = myCont.clust.color, min)

    if(argVerbose){
      cat("# --> min value per cluster :\n")
      print(myCont.minByClust.agg)
      cat("# --> max value per cluster :\n")
      print(myCont.maxByClust.agg)
      cat("\n# ------------------------------------\n")
    }

    hist(myCont.clust.color$values[myCont.notNA.idx], breaks = unique(c( min(myCont.minByClust.agg$values),
                                                                         sort(myCont.maxByClust.agg$values,
                                                                                decreasing = F))),
         xlab = "value",
         xlim = c(min(myData[, strNumVar], na.rm = T), max(myData[, strNumVar], na.rm = T)),
         ann = FALSE, mgp = c(3, 0.6, 0), cex.axis = myCexAxis, las = 2)
    mtext(side = 3, text = paste(nbrClusters, "self-disc"), line = 0.2, cex = myCexAxis)

    if(iGraphCount == (nbRowGraph*nbColGraph)){
        # par(mfrow=c(1,1))

      iGraphCount = 1
    }

  }

  if(class(argInData) == "character"){dev.off()}

}
