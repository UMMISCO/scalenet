#' scaleNet.convertToScaleNetFormat
#'
#' Convert into a scaleNet edgeList
#'
#' @param ioSubEnv a global environment variable

scaleNet.convertToScaleNetFormat <- function(ioSubEnv, iVarSpl){

  # For all pairs x-y, create a scaleNet key
  # --> get the properties
  inputData.header <- colnames(ioSubEnv$allData)

  # Load the globalNet
  tmp.dirPath <- file.path(ioSubEnv$output.dirPath, "globalGraph")
  if(!dir.exists(tmp.dirPath)){stop("# Not global network found!")}

  # iVarSpl = 1
#   tmp.prefix <- ifelse((ioSubEnv$subset.select == "spectralFuzzyCmeansSample"),
#                        paste("_globalNet_", iVarSpl, ".txt", sep=''),
#                        paste("_globalNet", ioSubEnv$recons.method,"txt", sep = '.'))
#   globalNet.filePath <- file.path(tmp.dirPath, gsub(".txt", tmp.prefix, basename(ioSubEnv$inputData.filePath)))
  tmp.prefix <- ifelse((ioSubEnv$subset.select == "spectralFuzzyCmeansSample"),
                       paste("globalNet_", iVarSpl, ".txt", sep=''),
                       paste("globalNet", ioSubEnv$recons.method,"txt", sep = '.'))
  globalNet.filePath <- file.path(tmp.dirPath, tmp.prefix)
  globalNet <- suppressWarnings(data.table::fread(input = globalNet.filePath, header = T, sep = "\t",
                                                  stringsAsFactors = F, data.table = F, showProgress = F))
  rownames(globalNet) <- paste(globalNet[,1], globalNet[,2], sep = "_<<_")

  # Create a scaleNet-like edgeList output (template)
  # ------
  edges.scale.col <- c("x", "y", "epresenceScore", "epresence", "eorientScore", "eorient", "ecorr")
  edges.scale <- as.data.frame(matrix(NA, ncol=length(edges.scale.col), nrow=nrow(globalNet)))
  colnames(edges.scale) <- edges.scale.col; rownames(edges.scale) <- rownames(globalNet)
  edges.scale[, "x"] <- globalNet[,"x"]
  edges.scale[, "y"] <- globalNet[,"y"]
  # Define the edgeList filename
  ioSubEnv$edgesList.fileName <- "edgesList.txt"

  # Build different global graph depending on the growing presence frequency threshold
  tmp.presFreq.thresh <- c()
  if(ioSubEnv$subset.select %in% c("spectral", "random")){

    tmp.presFreq.thresh <- ioSubEnv$presFreqThresh

  } else if(ioSubEnv$subset.select %in% c("spectralKmeans", "spectralFuzzyCmeansOrder", "spectralBipartition", "spectralFuzzyCmeansSample")){

    tmp.presFreq.thresh <- c(1)
  }

  # iPresFreq = 1
  for(iPresFreq in tmp.presFreq.thresh){

    # Set the output directory path
    # tmp.dirPath <- gsub(".txt", paste("_presFreq", iPresFreq, sep = ''), globalNet.filePath)
    tmp.dirPath <- gsub( paste(".", ioSubEnv$recons.method, ".txt", sep = ''), paste("_presFreq", iPresFreq, sep = ''), globalNet.filePath)
    if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}

    # Get the edges with more than iPresFreq
    globalNet.edges.idx <- which(globalNet[, 'presence.freq'] >= iPresFreq)
    if(length(globalNet.edges.idx)==0){warning("# --Wlib2-- No edges in global net!"); next;}

    # Extract the corresponding rows from the globalNet
    tmp.globalNet.subFreq <- globalNet[globalNet.edges.idx,]
    if(nrow(tmp.globalNet.subFreq)==0){break;}

    # Set a scaleNet-like edgeList output
    # ------
    # Make a copy of the 'template' discoNet edgeList
    tmp.edges.scale <- edges.scale

    # Set all edges as phantom
    # Then, for each learned edge, set the presence to 1
    tmp.edges.scale[, "epresence"] <- rep(0, nrow(tmp.edges.scale))
    tmp.edges.scale[rownames(tmp.globalNet.subFreq), "epresence"] <- 1

    # Set for all inferred edges the orientation
    # ----
    # --> set nonoriented pairs
    tmp.pairs <- tmp.globalNet.subFreq[which(tmp.globalNet.subFreq[, "presence.ort"] == 1), c("x", "y")]
    if(nrow(tmp.pairs)>0){tmp.edges.scale[rownames(tmp.pairs), "eorient"] <- 1}
    # --> set forward oriented pairs
    tmp.pairs <- tmp.globalNet.subFreq[which(tmp.globalNet.subFreq[, "presence.ort"] == 2), c("x", "y")]
    if(nrow(tmp.pairs)>0){tmp.edges.scale[rownames(tmp.pairs), "eorient"] <- 2}
    # --> set backward oriented pairs
    tmp.pairs <- tmp.globalNet.subFreq[which(tmp.globalNet.subFreq[, "presence.ort"] == -2), c("x", "y")]
    if(nrow(tmp.pairs)>0){tmp.edges.scale[rownames(tmp.pairs), "eorient"] <- -2}

    # Set for all inferred edges the score if given (ie: not all NA, as for bayes_hc)
    # NB1: average over all scores
    # NB2: should be identical if not influence by neighbours
    # ----
    # --> get pairs with a score
    tmp.pairs <- tmp.globalNet.subFreq[which(!is.na(tmp.globalNet.subFreq[, "epresenceScore"])), c("x", "y")]
    if(!all(is.na(tmp.globalNet.subFreq$epresenceScore))){

      tmp.epresenceScore.avg <- tmp.globalNet.subFreq[, "epresenceScore"]
      if(class(tmp.globalNet.subFreq$epresenceScore) == "character")
      {
        tmp.epresenceScore.list <- sapply(tmp.globalNet.subFreq[, "epresenceScore"], strsplit, ",")
        tmp.epresenceScore.list <- sapply(tmp.epresenceScore.list, unique)
        tmp.epresenceScore.list <- suppressWarnings(sapply(tmp.epresenceScore.list, as.numeric))
        tmp.epresenceScore.avg <- unlist(lapply(tmp.epresenceScore.list, mean, na.rm = TRUE))
      }
      tmp.edges.scale[rownames(tmp.globalNet.subFreq), "epresenceScore"] <- tmp.epresenceScore.avg
    }

    # Compute pairwise correlation (to get a sign)
    tmp.pair.idx <- which(tmp.edges.scale[, "epresence"] != 0)
    tmp.corr <- apply(tmp.edges.scale[tmp.pair.idx,], MARGIN = 1, function(myRow){
      cor(ioSubEnv$allData[,myRow["x"]], ioSubEnv$allData[,myRow["y"]], method = "spearman", use = "pairwise.complete.obs")
    })
    tmp.edges.scale[tmp.pair.idx, "ecorr"] <- tmp.corr

    # Save the edge list
    write.table(file=file.path(tmp.dirPath, gsub(".txt", paste(".", ioSubEnv$recons.method, ".txt", sep = ''), ioSubEnv$edgesList.fileName)), tmp.edges.scale,
                col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
  }
}
