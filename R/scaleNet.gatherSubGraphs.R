#' scaleNet.gatherSubGraphs
#'
#' Gather the sub graphs
#' - if the input datasets are based on spectralFuzzyCmeansSample
#' - save a globalNet for each assignment (up to ioSubEnv$nbVarSample)
#'
#'
#' @param ioSubEnv a global environment variable

scaleNet.gatherSubGraphs <- function(ioSubEnv, iVarSpl){

  # --!-- Time
  time.globalReconstruction.start <- proc.time()

  # Set all possible pairs for the global network
  allPairs <- combn(colnames(ioSubEnv$allData), 2)

  # Define an output table for the gathered information
  # --> here, the counts are made on the pair presence, orientation, or possible (ie. pair in the dataset)
  cat("# Set globalNet....\n")
  globalNet.colnames <- c("presence", "forward", "backward", "possible")
  globalNet <- matrix(0, ncol = length(globalNet.colnames), nrow = ncol(allPairs))
  colnames(globalNet) <- globalNet.colnames
  rownames(globalNet) <- paste(allPairs[1,], allPairs[2,], sep = "_<<_")

  if(ioSubEnv$verbose){
    cat("# --|", paste(head(colnames(globalNet)), collapse = ','), "...\n")
    cat("# --| nbPairs:", dim(globalNet)[1],"\n")
  }

  # Keep complementary info
  # --> For all reconstruction method, keep the eigenvector or cluster number that gives the presence
  globalNet.cmpl.colnames <- c("eigen", "cluster", "epresenceScore")

  # --> Create a data frame for complementary information with the rownames corresponding to globalNet
  globalNet.cmpl <- as.data.frame(matrix(ncol = length(globalNet.cmpl.colnames), nrow = nrow(globalNet)))
  rownames(globalNet.cmpl) <- rownames(globalNet)
  colnames(globalNet.cmpl) <- globalNet.cmpl.colnames

  # Set the directory path where all subnetworks can be found
  tmp.dirPath <- file.path(ioSubEnv$output.dirPath, "subGraphs")
  if(dir.exists(tmp.dirPath) == FALSE){stop("# Not local network found!")}

  if(ioSubEnv$subset.select %in% c("spectral", "random")){

    # iEigen <- 2
    for(iEigen in c(2:ncol(ioSubEnv$eigen.vectors))){

      # strSign <- "pos"
      for(strSign in c("pos", "neg")){

        # Define the network summary file path
        summary.fileName <- paste("edgesList", ioSubEnv$recons.method, "txt", sep = '.')

        subOutData.summaryPath.template = file.path(tmp.dirPath, "subOutput",
                                                    paste("eigenVector", "XXX_IEIGEN_XXX", "_",
                                                          ioSubEnv$subset.m, strSign, sep=''),
                                                    summary.fileName)
        # Load the sub graph summary
        tmp.subNet.filePath <- gsub( "XXX_IEIGEN_XXX", as.character(iEigen), subOutData.summaryPath.template)
        if(!file.exists(tmp.subNet.filePath)){next}
        tmp.subNet <- read.table(file = tmp.subNet.filePath, header = TRUE, row.names = 1, sep = '\t', as.is = TRUE)

        # Create a new format for the sub graph summary that fits the globalNet data file format
        # --> make the addition easier!!
        # --> this is why the colnames order of sub datasets is important!!
        tmp.subNet.newFormat <- matrix(0, ncol = ncol(globalNet), nrow = nrow(tmp.subNet) )
        # --> same columns as globalNet
        colnames(tmp.subNet.newFormat) <- colnames(globalNet)
        # --> same row as subNet
        rownames(tmp.subNet.newFormat) <- paste(tmp.subNet[,"x"], tmp.subNet[,"y"], sep = "_<<_")

        # Set all pairs of the subgraph as 'possible' (ie. visited)
        tmp.subNet.newFormat[,"possible"] <- rep(1,nrow(tmp.subNet.newFormat))

        # Set the presence of edges, and their possible orientation
        if(length(which(tmp.subNet[,"epresence"] == 1)) > 0){

          # If an edge has been inferred, set the 'presence' to 1
          tmp.subNet.newFormat[(tmp.subNet[,"epresence"] == 1),"presence"] <- 1

          # Set also the orientation
          # --> as the variable names order is respected, the fwd/bck sign is also the same
          if(length(which(tmp.subNet[,"eorient"] == 2)) > 0){
            tmp.subNet.newFormat[(tmp.subNet[,"eorient"] == 2),"forward"] <- 1
          }

          if(length(which(tmp.subNet[,"eorient"] == -2)) > 0){
            tmp.subNet.newFormat[(tmp.subNet[,"eorient"] == -2),"backward"] <- 1
          }
        }

        # Find the rownames of the new format subgraph in globalNet data file
        tmp.match <- match(rownames(tmp.subNet.newFormat), rownames(globalNet))

        # Add the counts from the subgraph to global net
        globalNet[tmp.match,] <-(globalNet[tmp.match,] + tmp.subNet.newFormat[c(1:length(tmp.match)),])

        # Keep also the eigen vector number and the edge score into a separated data frame when an edge has been inferred
        tmp.presence.key <- rownames(tmp.subNet.newFormat)[(tmp.subNet[,"epresence"] == 1)]
        if(length(tmp.presence.key) > 0){
          globalNet.cmpl[tmp.presence.key, "eigen"] <- paste(globalNet.cmpl[tmp.presence.key, "eigen"],
                                                             rep(iEigen, length(tmp.presence.key)), sep = ',')

          globalNet.cmpl[tmp.presence.key, "epresenceScore"] <- paste(globalNet.cmpl[tmp.presence.key, "epresenceScore"],
                                                              tmp.subNet[tmp.presence.key, "epresenceScore"], sep = ',')
        }
      }
    }

  } else if(ioSubEnv$subset.select %in% c("spectralKmeans", "spectralFuzzyCmeansOrder",
                                          "spectralBipartition", "spectralFuzzyCmeansSample")){
    # Set the cluster numbers
    nbClusters <- c()
    if(ioSubEnv$subset.select == "spectralKmeans"){
      nbClusters <- seq_len(ioSubEnv$kmeans.k)

    } else if(ioSubEnv$subset.select %in% c("spectralFuzzyCmeansOrder", "spectralFuzzyCmeansSample")){
      nbClusters <- seq_len(ioSubEnv$cmeans.k)

    } else if(ioSubEnv$subset.select == "spectralBipartition"){
      nbClusters <- unique(ioSubEnv$specBi.allClusters)
    }

    # iClust = nbClusters[1]
    for(iClust in nbClusters){

      # Define the network summary file path
      summary.fileName <- paste("edgesList", ioSubEnv$recons.method, "txt", sep = '.')

      subOutData.summaryPath.template <- ""
      if(ioSubEnv$subset.select == "spectralFuzzyCmeansSample"){

        subOutData.summaryPath.template <- file.path(tmp.dirPath, "subOutput",
                                                     paste(ioSubEnv$subset.select, "XXX_ICLUST_XXX", "_",
                                                           ioSubEnv$subset.m, "_", iVarSpl, sep=''),
                                                     summary.fileName)

      } else {

        subOutData.summaryPath.template <- file.path(tmp.dirPath, "subOutput",
                                                     paste(ioSubEnv$subset.select, "XXX_ICLUST_XXX", "_",
                                                           ioSubEnv$subset.m, sep=''),
                                                     summary.fileName)
      }

      # Load the sub graph summary
      tmp.subNet.filePath <- gsub( "XXX_ICLUST_XXX", as.character(iClust), subOutData.summaryPath.template)
      if(!file.exists(tmp.subNet.filePath)){next}
      tmp.subNet <- read.table(file = tmp.subNet.filePath, header = TRUE, row.names = 1, sep = '\t', as.is = TRUE)

      # Create a new format for the sub graph summary that fits the globalNet data file format
      # --> make the addition easier!!
      # --> this is why the colnames order of sub datasets are important!!
      tmp.subNet.newFormat <- matrix(0, ncol = ncol(globalNet), nrow = nrow(tmp.subNet) )
      # --> same columns as globalNet
      colnames(tmp.subNet.newFormat) <- colnames(globalNet)

      # --> same row as subNet
      rownames(tmp.subNet.newFormat) <- paste(tmp.subNet[,"x"], tmp.subNet[,"y"], sep = "_<<_")

      # Set all pairs of the subgraph as 'possible' (ie. visited)
      tmp.subNet.newFormat[,"possible"] <- rep(1,nrow(tmp.subNet.newFormat))

      # Set the presence of edges, and their possible orientation
      if(length(which(tmp.subNet[,"epresence"] == 1)) > 0){

        # If an edge has been inferred, set the 'presence' to 1
        tmp.subNet.newFormat[(tmp.subNet[,"epresence"] == 1),"presence"] <- 1

        # Set also the orientation
        # --> as the variable names order is respected, the fwd/bck sign is also the same
        if(length(which(tmp.subNet[,"eorient"] == 2)) > 0){
          tmp.subNet.newFormat[(tmp.subNet[,"eorient"] == 2),"forward"] <- 1
        }

        if(length(which(tmp.subNet[,"eorient"] == -2)) > 0){
          tmp.subNet.newFormat[(tmp.subNet[,"eorient"] == -2),"backward"] <- 1
        }
      }

      # Find the rownames of the new format subgraph in globalNet data file
      tmp.match <- match(rownames(tmp.subNet.newFormat), rownames(globalNet))

      # Add the counts from the subgraph to global net
      globalNet[tmp.match,] <-(globalNet[tmp.match,] + tmp.subNet.newFormat[c(1:length(tmp.match)),])

      # Keep also the cluster number into a separated data frame when an edge has been inferred
      tmp.presence.key <- rownames(tmp.subNet.newFormat)[(tmp.subNet[,"epresence"] == 1)]
      if(length(tmp.presence.key) > 0){
        globalNet.cmpl[tmp.presence.key, "cluster"] <- paste(globalNet.cmpl[tmp.presence.key, "cluster"],
                                                             rep(iClust, length(tmp.presence.key)), sep = ',')

        globalNet.cmpl[tmp.presence.key, "epresenceScore"] <- paste(globalNet.cmpl[tmp.presence.key, "epresenceScore"],
                                                            tmp.subNet[tmp.presence.key, "epresenceScore"], sep = ',')
      }
    }
  }

  # Compute supplementary stats from the gathered data
  # --> presence.freq: frequency of edge 'presence' over the number of edge 'possible'
  # --> presence.ort: the consensus orientation (using a simple max rule...no orientation if 2 max!)
  globalNet.supp.col <- c("presence.freq", "presence.ort")
  globalNet.supp <- as.data.frame(matrix(NA, ncol = length(globalNet.supp.col), nrow = nrow(globalNet)))
  rownames(globalNet.supp) <- rownames(globalNet)
  colnames(globalNet.supp) <- globalNet.supp.col

  # --> presence.freq:
  globalNet.supp[, "presence.freq"] <- (globalNet[, "presence"]/globalNet[, "possible"])

  # --> presence.ort:
  max.ort.list <- plyr::alply(globalNet[,c("forward", "backward")], 1, function(x){which(x==max(x))})
  names(max.ort.list) <- c()
  max.ort.list.length <- unlist(lapply(max.ort.list, length))

  # single max ==>
  # --> set the name of the corresponding orientation type (forward, backward)
  # --> then replace with integer 2, -2
  if(length(which(max.ort.list.length == 1))>0){

    globalNet.supp[which(max.ort.list.length == 1), "presence.ort"] <- names(unlist(max.ort.list[which(max.ort.list.length == 1)]))

    tmp.idx <- which(globalNet.supp[, "presence.ort"] == "forward")
    if(length(tmp.idx) > 0){ globalNet.supp[tmp.idx, "presence.ort"] <- rep("2", length(tmp.idx)) }

    tmp.idx <- which(globalNet.supp[, "presence.ort"] == "backward")
    if(length(tmp.idx) > 0){ globalNet.supp[tmp.idx, "presence.ort"] <- rep("-2", length(tmp.idx)) }

    globalNet.supp[, "presence.ort"] <- as.numeric(globalNet.supp[, "presence.ort"])
  }

  # Set 1 for the orientation otherwise
  tmp.idx <- which(is.na(globalNet.supp[, "presence.ort"]) & (globalNet.supp[, "presence.freq"] > 0) )
  if(length(tmp.idx) > 0){ globalNet.supp[tmp.idx, "presence.ort"] <- rep(1, length(tmp.idx)) }

  # --> presence eigenvect pos/neg
  # ----> remove first "NA," for !NA data
  tmp.col <- ifelse((ioSubEnv$subset.select %in% c("spectral", "random")), "eigen", "cluster")
  tmp.key <- rownames(globalNet.cmpl)[which(!is.na(globalNet.cmpl[, tmp.col]))]
  tmp.vect <- sapply(globalNet.cmpl[tmp.key, tmp.col], gsub, pattern = "^NA,", replacement = "")
  names(tmp.vect) <- tmp.key
  globalNet.cmpl[tmp.key, tmp.col] <- tmp.vect[tmp.key]

  # --> epresenceScore
  # ----> remove first "NA," for !NA data
  tmp.col <- c("epresenceScore")
  tmp.key <- rownames(globalNet.cmpl)[which(!is.na(globalNet.cmpl[, tmp.col]))]
  tmp.vect <- sapply(globalNet.cmpl[tmp.key, tmp.col], gsub, pattern = "^NA,", replacement = "")
  names(tmp.vect) <- tmp.key
  globalNet.cmpl[tmp.key, tmp.col] <- tmp.vect[tmp.key]

  # Use the rownames to write the x and y columns, and cbind all the information
  globalNet <- cbind.data.frame(data.frame(do.call('rbind', strsplit(rownames(globalNet),"_<<_",
                                                                     fixed=TRUE))),
                                globalNet, globalNet.supp, globalNet.cmpl)
  colnames(globalNet) <- c("x", "y",
                           globalNet.colnames,
                           colnames(globalNet.supp),
                           colnames(globalNet.cmpl))
  # --!-- Time
  ioSubEnv$gather.time <- (proc.time() - time.globalReconstruction.start)
  ioSubEnv$total.time <- ioSubEnv$total.time + ioSubEnv$gather.time

  # Create output directory for the global graph
  tmp.dirPath <- file.path(ioSubEnv$output.dirPath, "globalGraph")
  if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}

#   tmp.prefix <- ifelse((ioSubEnv$subset.select == "spectralFuzzyCmeansSample"),
#                        paste("_globalNet_", iVarSpl, ".txt", sep=''),
#                        paste("_globalNet", ioSubEnv$recons.method,"txt", sep = '.'))
#
#   globalNet.filePath <- file.path(tmp.dirPath, gsub(".txt", tmp.prefix, basename(ioSubEnv$inputData.filePath)))
  tmp.prefix <- ifelse((ioSubEnv$subset.select == "spectralFuzzyCmeansSample"),
                       paste("globalNet_", iVarSpl, ".txt", sep=''),
                       paste("globalNet", ioSubEnv$recons.method,"txt", sep = '.'))

  globalNet.filePath <- file.path(tmp.dirPath, tmp.prefix)

  write.table(globalNet, file = globalNet.filePath, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}
