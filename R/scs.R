#' scs
#'
#' Create a consensus from scaleNet embedded method results
#'
#' @param workspaceDir output directory file path (to be created)
#' @param argInData input data file path (a tab separated file) or a data.frame, with features in columns and observations in rows
#' @param argReconsMeth string vector with reconstruction method name
#' @param argReconsMethInfo list of complementary information on each reconstruction method. The
#' provided information are: "ort" (the method provides orientation n/y) and edge weight "eweight"
#' (epresenceScore/ecorr)
#' @param argEmbReconsParam list specific parameters for each reconstruction method and complementary parameters
#' specific to scalenet, namely: "nbSamples", "numSeed", "nbCPU", "eigenPerc", "varPerc" and "subsetType"
#' @param argPresFreqThresh a numeric vector of threshold for the frequency of edges from all subgraphs
#' @param clean.workspace if FALSE, all output files are kept
#' @param argVerbose display details
#'
#' @export

# Error code from 10,000 to 10,099
# ----------

# scs( workspaceDir = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck_spectral",
#      argInData = "~/Projects/Projects_largeScale/data/benchmark/andes/input/rawData/andes_20000/andes_20000_0001.txt",
#      argReconsMeth = c("aracne", "bayes_hc"),
#      argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"),
#                               bayes_hc = list(ort = "y", eweight = "ecorr")),
#      argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001),
#                              bayes_hc = list(score="bde", restart=21), varPerc = 0.2),
#      argPresFreqThresh = c(0.3, 0.8), clean.workspace = TRUE, discretize = TRUE, argVerbose = TRUE)

# scs( workspaceDir = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/output/clinicalMGS_foldChange_raw_sigNameBigTargets",
#      argInData = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/input/clinicalMGS_foldChange_raw_sigNameBigTargets.txt",
#      argReconsMeth = c("aracne", "bayes_hc"),
#      argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"), bayes_hc = list(ort = "y", eweight = "ecorr")),
#      argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001), bayes_hc = list(score="bde", restart=21), varPerc = 0.2),
#      argPresFreqThresh = c(0.3, 0.8), clean.workspace = TRUE, argDiscretize = TRUE, argVerbose = TRUE)

# inlineData <- read.table(file="~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/input/clinicalMGS_foldChange_raw_sigNameBigTargets.txt",
#                          header = TRUE, sep = '\t', as.is = TRUE)
# scs( workspaceDir = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/output/clinicalMGS_foldChange_raw_sigNameBigTargets",
#      argInData = inlineData,
#      argReconsMeth = c("aracne", "bayes_hc"),
#      argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"), bayes_hc = list(ort = "y", eweight = "ecorr")),
#      argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001), bayes_hc = list(score="bde", restart=21), varPerc = 0.2),
#      argPresFreqThresh = c(0.3, 0.8), clean.workspace = FALSE, argDiscretize = TRUE, argVerbose = FALSE)
#
# scs( workspaceDir = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/output/clinicalMGS_foldChange_raw_sigNameBigTargets",
#      argInData = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/input/clinicalMGS_foldChange_raw_sigNameBigTargets.txt",
#      argReconsMeth = c("aracne", "bayes_hc"),
#      argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"), bayes_hc = list(ort = "y", eweight = "ecorr")),
#      argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001), bayes_hc = list(score="bde", restart=21), varPerc = 0.2),
#      argPresFreqThresh = c(0.3, 0.8), clean.workspace = FALSE, argDiscretize = TRUE, argVerbose = FALSE)

scs <- function( workspaceDir, argInData,
                 argReconsMeth = c("aracne", "bayes_hc"),
                 argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"),
                                          bayes_hc = list(ort = "y", eweight = "ecorr")),
                 argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001),
                                          bayes_hc = list(score="bde", restart=21),
                                          varPerc = 0.2),
                 argPresFreqThresh = c(0.5, 1), clean.workspace = TRUE,
                 argDiscretize = FALSE, argVerbose = FALSE ){

  # Check if the scalenet outputs already exists
  # if not, do scalenet for the method given to scs
  # ----------------------------------------------------
  #if(!file.exists(workspaceDir)){
  if(!file.exists(paste(workspaceDir,"globalGraph",sep=""))){

    # Check if argRconsParam is given
    # if not, quit...
    if(is.null(argEmbReconsParam)){stop("# --Err-- 1000")}

    # Loop on the reconstruction methods
    # strMeth = argReconsMeth[1]
    for(strMeth in argReconsMeth){

      scalenet( argInData = argInData,
                argOutDir = workspaceDir,
                argReconsMeth = strMeth,
                argReconsParam = argEmbReconsParam[[strMeth]],
                argPresFreqThresh = argPresFreqThresh,
                argNbSamples = argEmbReconsParam[["nbSamples"]],
                argNumSeed = argEmbReconsParam[["numSeed"]],
                argNbCPU = argEmbReconsParam[["nbCPU"]],
                argEigenPerc = argEmbReconsParam[["eigenPerc"]],
                argVarPerc = argEmbReconsParam[["varPerc"]],
                argSubsetType = argEmbReconsParam[["subsetType"]],
                argDiscretize = argDiscretize,
                argVerbose = argVerbose)
    }
  }

  # Get all variable names, pairs and make labels for the edgesList/adjMat template files
  # ----------------------------------------------------
  # --> get the properties
  argInData.header <- NULL
  if(class(argInData) == "character"){

    con <- file(argInData, open = "r")
    argInData.header <- readLines(con, n=1)
    close(con)
    argInData.header <- unlist(strsplit(argInData.header, split = "\t"))

  } else if(class(argInData) == "data.frame") {

    argInData.header <- colnames(argInData)

  } else { stop("# --Err-- 1002") }

  # --> generare all pairs
  allPairs <- combn(argInData.header, 2)
  allPairs.key <- paste(allPairs[1,], allPairs[2,], sep = "_<<_")

  # Create a scaleNet-like edgeList output (template)
  # ------
  edges.scale.col <- c("x", "y", "epresenceScore", "epresence", "eorientScore", "eorient", "ecorr")
  edges.scale <- as.data.frame(matrix(NA, ncol=length(edges.scale.col), nrow=ncol(allPairs)))
  colnames(edges.scale) <- edges.scale.col; rownames(edges.scale) <- allPairs.key
  edges.scale[, "x"] <- allPairs[1,]
  edges.scale[, "y"] <- allPairs[2,]
  # Define the edgeList filename
  edgesList.fileName <- "edgesList.txt"

  # Create a output consensus directory for this scalenet output
  # ------
  tmp.consensus.dirPath <- file.path(workspaceDir, "consensusGraph")
  if(!dir.exists(tmp.consensus.dirPath)){dir.create(tmp.consensus.dirPath)}

  # Create a directory path template to load the network learnt by scalenet
  #   workspaceDir.filtered.template <- file.path( workspaceDir, "globalGraph",
  #                                           paste(gsub(".txt", "", basename(argInData)),
  #                                                 "_globalNet_presFreq_XXX_PRESFREQ_XXX", sep = ""))
  workspaceDir.filtered.template <- file.path( workspaceDir, "globalGraph", "globalNet_presFreq_XXX_PRESFREQ_XXX")

  # Loop of filtered global subgraphes
  # iPresFreq <- 0.8
  for(iPresFreq in argPresFreqThresh){

    # Update the presFreq
    workspaceDir.filtered.freq <- gsub("_XXX_PRESFREQ_XXX", iPresFreq, workspaceDir.filtered.template)

    # Create a consensus directory for this edge frequency
    tmp.consensus.filtered.dirPath <- file.path(tmp.consensus.dirPath,
                                                gsub("_globalNet_", "_consensusNet", basename(workspaceDir.filtered.freq)))
    if(!dir.exists(tmp.consensus.filtered.dirPath)){dir.create(tmp.consensus.filtered.dirPath)}

    # Prepare a rank list, an orientation list and a totalNbrEdges list
    score.list <- vector("list", length(argReconsMeth)); names(score.list) <- argReconsMeth
    orient.list <- vector("list", length(argReconsMeth)); names(orient.list) <- argReconsMeth
    corr.list <- vector("list", length(argReconsMeth)); names(corr.list) <- argReconsMeth

    # Keep all the labels of the edges discovered by at least one recontruction method
    tmp.consensus.edgeLabels <- c()

    # Get the edge score from each reconstruction method
    # Order the scores following the given metric
    # strMeth <- argReconsMeth[1]
    for(strMeth in argReconsMeth){

      tmp.recons.summary.filePath <- file.path(workspaceDir.filtered.freq,
                                               gsub(".txt", paste(".", strMeth, ".txt", sep = ''), edgesList.fileName))

      if(file.exists(tmp.recons.summary.filePath)){

        # Load the data
        tmp.data <- read.table(tmp.recons.summary.filePath, header = TRUE, as.is = TRUE, sep = '\t')
        # head(tmp.data);dim(tmp.data)

        # Get the index of edges
        tmp.edges.idx <- which(tmp.data$epresence == 1)

        if(length(tmp.edges.idx)>0){

          # Keep only the inferred edges
          tmp.data <- tmp.data[tmp.edges.idx,]

          # Order by the given measure (always take abs measure)
          tmp.data <- tmp.data[order(x = abs(tmp.data[, argReconsMethInfo[[strMeth]][["eweight"]]]),
                                     decreasing = TRUE),]

          # Get the scores
          score.list[[strMeth]] <- tmp.data[, argReconsMethInfo[[strMeth]][["eweight"]]]
          names(score.list[[strMeth]]) <- rownames(tmp.data)

          # Get the orientations
          orient.list[[strMeth]] <- tmp.data[, "eorient"]
          names(orient.list[[strMeth]]) <- rownames(tmp.data)

          # Get the correlations
          corr.list[[strMeth]] <- tmp.data[, "ecorr"]
          names(corr.list[[strMeth]]) <- rownames(tmp.data)

          # Insert these labels in the consensus
          tmp.consensus.edgeLabels <- unique(c(tmp.consensus.edgeLabels, rownames(tmp.data)))
        }
      }
    }

    # Create an empty consensus table
    tmp.consensus.col <- c(paste(argReconsMeth,"score", sep = '.'),
                           paste(argReconsMeth,"ort",sep = '.'), "consensus.score.mean", "consensus.ort.wmean", "ecorr")
    tmp.consensus <- matrix(NA, ncol=length(tmp.consensus.col), nrow = length(tmp.consensus.edgeLabels))
    rownames(tmp.consensus) <- tmp.consensus.edgeLabels; colnames(tmp.consensus) <- tmp.consensus.col
    # head(tmp.consensus); dim(tmp.consensus)

    # strMeth <- argReconsMeth[1]
    for(strMeth in argReconsMeth){

      if(all(is.na(score.list[[strMeth]]))){next;}

      # Set the rank corresponding to ordered scores in consensus table
      tmp.consensus[names(score.list[[strMeth]]), paste(strMeth, "score", sep = '.')] <- seq_len(length(score.list[[strMeth]]))

      # Set the orientations corresponding to ordered scores in consensus table
      tmp.consensus[names(score.list[[strMeth]]),paste(strMeth, "ort", sep = '.')] <- orient.list[[strMeth]]

      # Set the correlations corresponding to ordered scores in consensus table
      tmp.consensus[names(score.list[[strMeth]]), "ecorr"] <- corr.list[[strMeth]]
      # head(tmp.consensus); dim(tmp.consensus)

      # Set a default rank for edges that have not been inferred by the method
      # (corresponding to max + 1)
      tmp.max <- (max(tmp.consensus[, paste(strMeth, "score", sep = '.')], na.rm = T) + 1)
      tmp.noInf.idx <- which(!rownames(tmp.consensus) %in% names(score.list[[strMeth]]))
      if(length(tmp.noInf.idx)>0){
        tmp.consensus[tmp.noInf.idx, paste(strMeth, "score", sep = '.')] <- tmp.max
      }
      # Rescale to the max
      tmp.consensus[, paste(strMeth, "score", sep = '.')] <- (tmp.consensus[, paste(strMeth, "score", sep = '.')]/tmp.max)
      # head(tmp.consensus); dim(tmp.consensus)
    }

    # Set the average score in the consensus.score.mean
    # ----
    tmp.score.col <- grep(pattern = '.score$', tmp.consensus.col)
    tmp.consensus[, "consensus.score.mean"] <- unlist(apply(tmp.consensus[, tmp.score.col], 1, mean))

    # Set the average orientation in the consensus.ort.wmean
    # ----
    reconsMeth.orient.vect <- unlist(lapply(argReconsMethInfo, `[[`, "ort"))
    tmp.ortMeth.idx <- which(reconsMeth.orient.vect == "y")

    if(length(tmp.ortMeth.idx) > 0){

      tmp.score.col <- paste(names(reconsMeth.orient.vect)[tmp.ortMeth.idx], "score", sep = '.')
      tmp.orient.col <- paste(names(reconsMeth.orient.vect)[tmp.ortMeth.idx], "ort", sep = '.')

      # Replace 1 with 0 in the tmp.orient.col columns
      tmp.consens.orient.col <- tmp.consensus[,tmp.orient.col, drop = FALSE]
      tmp.consens.orient.col[tmp.consens.orient.col == 1] <- 0
      tmp.consensus[,tmp.orient.col] <- tmp.consens.orient.col[,tmp.orient.col]
      rm(tmp.consens.orient.col)

      # Replace 2/-2 with 1/-1 in the tmp.orient.col columns
      tmp.consensus[,tmp.orient.col] <- (tmp.consensus[,tmp.orient.col]/2)

      # Multiply each (1-score) column by the corresponding orientation column,
      # Sum over the mutiplications
      # and divide by sum(1-scores)
      tmp.weigtedOrt <- unlist(apply((1-tmp.consensus[,tmp.score.col, drop = FALSE])*tmp.consensus[,tmp.orient.col, drop = FALSE], 1, sum, na.rm = T))
      tmp.sumScore <- unlist(apply((1-tmp.consensus[,tmp.score.col, drop = FALSE]), 1, sum, na.rm = T))

      # It's possible that all methods that can give orientation did not predict edges...
      # Thus, they have score at 1 (ie weight at (0))
      # This will produce orientation NA that can be replace by 0
      tmp.consensus[, "consensus.ort.wmean"] <- (tmp.weigtedOrt/tmp.sumScore)
      tmp.sumScore.null.idx <- which(tmp.sumScore == 0)
      # if(length(tmp.sumScore.null.idx) > 0){tmp.consensus[tmp.sumScore.null.idx, "consensus.ort.wmean"] <- 1}
      if(length(tmp.sumScore.null.idx) > 0){tmp.consensus[tmp.sumScore.null.idx, "consensus.ort.wmean"] <- 0}

    } else {
      tmp.consensus[, "consensus.ort.wmean"] <- 0
    }

    # Order consensus by consensus.score.mean
    tmp.consensus <- tmp.consensus[order(tmp.consensus[, "consensus.score.mean"], decreasing = FALSE ),]
    # head(tmp.consensus)

    # Bind ranks and orientations
    avg.rank.ort.globalNet.df <- cbind.data.frame(data.frame(do.call('rbind',
                                                                     strsplit(as.character(rownames(tmp.consensus)),"_<<_",fixed=TRUE))),
                                                  tmp.consensus[, "consensus.score.mean"],
                                                  tmp.consensus[, "consensus.ort.wmean"],
                                                  tmp.consensus[, "ecorr"])
    colnames(avg.rank.ort.globalNet.df) <- c("x", "y", "avg.rank", "avg.ort", "ecorr")
    rownames(avg.rank.ort.globalNet.df) <- rownames(tmp.consensus)
    # head(avg.rank.ort.globalNet.df)

    # Write the consensus interactions
    write.table(avg.rank.ort.globalNet.df,
                file = file.path(tmp.consensus.filtered.dirPath, "edgesList.consensus.rawAvg.txt"),
                col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")


    # ----------------------------------------
    # Convert into a scalenet edgeList
    # ----------------------------------------
    if(file.exists(file.path(tmp.consensus.filtered.dirPath, "edgesList.consensus.rawAvg.txt")) &
       (!file.exists(file.path(tmp.consensus.filtered.dirPath, edgesList.fileName)))){

      # Load the consensus net
      consensusNet <- read.table(file = file.path(tmp.consensus.filtered.dirPath, "edgesList.consensus.rawAvg.txt"),
                                 header = T, sep = "\t", as.is = T)
      # head(consensusNet); dim(consensusNet)

      # Set a scalenet-like edgeList output
      # ------
      # Make a copy of the 'template' scalenet edgeList
      tmp.edges.scale <- edges.scale
      # head(tmp.edges.scale); dim(tmp.edges.scale)

      # Set all edges as absent
      tmp.edges.scale[, "epresence"] <- rep(0, nrow(tmp.edges.scale))

      # Then, for each learned edge, set the epresence to 1, keep the avg.rank and the avg.ort
      tmp.edges.scale[rownames(consensusNet), "epresence"] <- 1
      tmp.edges.scale[rownames(consensusNet), "epresenceScore"] <- consensusNet[rownames(consensusNet), "avg.rank"]
      tmp.edges.scale[rownames(consensusNet), "eorientScore"] <- consensusNet[rownames(consensusNet), "avg.ort"]
      tmp.edges.scale[rownames(consensusNet), "ecorr"] <- consensusNet[rownames(consensusNet), "ecorr"]

      # Set nonoriented pairs
      tmp.idx <- which(consensusNet[, "avg.ort"] == 0  & !is.na(consensusNet[, "avg.rank"]))
      if(length(tmp.idx) > 0){
        tmp.pairs <- rownames(consensusNet[tmp.idx, c("x", "y")])
        tmp.edges.scale[tmp.pairs, "eorient"] <- 0
      }

      # Get forward oriented pairs
      # ----
      tmp.idx <- which(consensusNet[, "avg.ort"] > 0  & !is.na(consensusNet[, "avg.rank"]))
      if(length(tmp.idx) > 0){
        tmp.pairs <- rownames(consensusNet[tmp.idx, c("x", "y")])
        tmp.edges.scale[tmp.pairs, "eorient"] <- 1
      }

      # Get backward oriented pairs
      # ----
      tmp.idx <- which(consensusNet[, "avg.ort"] < 0  & !is.na(consensusNet[, "avg.rank"]))
      if(length(tmp.idx) > 0){
        tmp.pairs <- rownames(consensusNet[tmp.idx, c("x", "y")])
        tmp.edges.scale[tmp.pairs, "eorient"] <- (-1)
      }

      # Save the edge list
      write.table(file=file.path(tmp.consensus.filtered.dirPath, edgesList.fileName),
                  tmp.edges.scale, col.names = TRUE, row.names = TRUE, sep = "\t",
                  quote = FALSE)
    }

  }

  # Load each consensus graph in a list, with each element of the list
  # corresponding to an edge frequency
  consensusGraph.list <- list()

  # Get the list of directory in the "consensusGraph" subdir
  tmpDirPath <- list.dirs(path = file.path(workspaceDir, "consensusGraph"),
                          full.names = TRUE, recursive = FALSE)

  for(dirPath in tmpDirPath){

    # Split to get the frequency
    tmpSplit <- unlist(strsplit(basename(dirPath), split = "_consensusNetpresFreq"))

    if(file.exists(file.path(dirPath, "edgesList.txt"))){
      # Insert in the consensus graph list
      # consensusGraph.list[[as.character(tmpSplit[2])]] <- read.table(file = file.path(dirPath, "edgesList.txt"),
      #                                                                header = TRUE, row.names = 1, as.is = TRUE)
      consensusGraph.list[[as.character(tmpSplit)]] <- read.table(file = file.path(dirPath, "edgesList.txt"),
                                                                     header = TRUE, row.names = 1, as.is = TRUE)
    } else {stop("# --Err-- 1001")}

  }

  if(clean.workspace == TRUE){unlink(workspaceDir, recursive = TRUE, force = FALSE)}

  return(consensusGraph.list)
}

