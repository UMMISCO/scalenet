#' rMethod.hc
#'
#' rMethod.hc reconstruct network using the hhill climbing approach (from bnlearn package)
#'
#' @param argInData input data file, a tab separated file with features in columns and observations in rows
#' @param argOutDir output directory file, will be created
#' @param argScore scoring function
#' @param argRestart search and score restart
#' @param argOrders number of random variable orderings to fuse (ordering consensus).
#' Hill climbing keeps the first best-scoring move it finds, so the order in which
#' bnlearn visits candidates decides every tie and steers the greedy path. The base
#' order is first made canonical from the variables' own data, which already makes
#' the result invariant to renaming and to column order; argOrders > 1 then averages
#' over that many random orderings so the network does not rest on one arbitrary
#' scan order either, keeping arcs that recur in at least argOrderFreq of them.
#' @param argOrderFreq fraction of the orderings an arc must appear in to be kept
#' @param argSeed seed for the orderings; the caller derives it from the run seed and
#' the subgraph identity, neither of which depends on variable names or column
#' positions, so the consensus stays reproducible without reintroducing a dependence
#' on the input layout
#' @param argVerbose display details
#'
#' @export

# rMethod.hc(argInData = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/output/clinicalMGS_foldChange_raw_sigNameBigTargets/subGraphs/subInput/eigenVector5_34pos.tsv",
#          argOutDir = "~/Projects/Projects_largeScale/data/clinicalMGS_foldChange/output/clinicalMGS_foldChange_raw_sigNameBigTargets/subGraphs/subOutput/eigenVector5_34pos",
#          argScore = "bde", argRestart = 20, argVerbose = TRUE)

rMethod.hc <- function(argInData, argOutDir, argScore = "bde", argRestart = 20, argVerbose = FALSE,
                       maxIter = 10000, argOrders = 1, argOrderFreq = 0.5, argSeed = 6196) {

  if(argVerbose){cat( "# --------\n# -> START Bayes HC...\n" )}

  if(argVerbose){
    #### Recall the main parameters
    cat( "# --------\n# Inputs:\n# ----\n"
         , "# Input data file -->", argInData, "\n"
         , "# Output directory -->", argOutDir, "\n" )
    # ----
  }

  #### Load the raw data file
  inputData.df <- read.table(file = argInData, header = TRUE, stringsAsFactor = FALSE, sep = "\t", check.names = FALSE)
  if(ncol(inputData.df)<=2){quit(save = "no", status = 0)}

  head(inputData.df)

  #### Remove columns with more than 80% NA
  tmp.col.idx <- which(unlist(apply(inputData.df, 2, function(x){ length(which(is.na(x))) }))/nrow(inputData.df) > 0.75)
  if(length(tmp.col.idx)>0){inputData.df <- inputData.df[,-tmp.col.idx]}

  #### Keep only complet cases
  naRows <- unique(which(is.na(inputData.df), arr.ind = T)[,1])
  if( length(naRows) > 0 ){ inputData.df <- inputData.df[-naRows,] }
inputData.df

  # Make sure all columns are factors
  inputData.df[, colnames(inputData.df)] <- as.data.frame(lapply(inputData.df[, colnames(inputData.df)] , factor))

  #### Recall the main parameters
  if(argVerbose){
    cat( "# All properties -->", paste( colnames(inputData.df), collapse = ", " ), "\n"
         , "# Score -->", argScore, "\n"
         , "# Nbr. restarts -->", argRestart, "\n"
         , "# --------\n" )
  }
  # ----

  #### The hc method requires that all variables have more than 1 level
  allLevels = sapply( inputData.df, nlevels )
  varWithSingleLevel = which(allLevels < 2)
  tmp.allProperties = colnames(inputData.df)

  #### Start time
  startTime <- proc.time()

  #### The data hc actually sees (single-level variables cannot be scored)
  hcData <- inputData.df
  if( length( varWithSingleLevel ) > 0 )
  {
    hcData <- inputData.df[, -varWithSingleLevel]
    if(argVerbose){cat("# ------| Properties with less than 2 levels: ", paste( tmp.allProperties[varWithSingleLevel], collapse = "," ), "\n" )}
    tmp.allProperties = tmp.allProperties[-varWithSingleLevel]
  }
  rm(tmp.allProperties)

  #### Canonical base order, derived from the DATA rather than from any label.
  #### Sorting the variables by their own value vector is invariant to renaming
  #### (the data does not move) *and* to column permutation (the sorted sequence
  #### is the same), so the orderings sampled below are identical in both cases
  #### and the invariance is exact rather than statistical. Ordering by name is
  #### invariant to column order but not to renaming; ordering by column position
  #### is the reverse -- names and positions are equally arbitrary labels.
  #### Columns with identical data are interchangeable and fall back to the name.
  tmp.dataKey <- vapply(hcData, function(col) paste(as.integer(col), collapse = ","), character(1))
  hcData <- hcData[, order(tmp.dataKey, colnames(hcData)), drop = FALSE]

  #### Ordering consensus. One hc run per random column ordering; arcs are
  #### tallied by variable NAME, so the tally itself does not depend on any
  #### ordering. argOrders = 1 reproduces the single-ordering behaviour.
  nOrders <- max(1, as.integer(argOrders))
  arcTally <- new.env(hash = TRUE, parent = emptyenv())
  bump <- function(key){
    cur <- mget(key, envir = arcTally, ifnotfound = list(0L))[[1]]
    assign(key, cur + 1L, envir = arcTally)
  }
  for(iOrder in seq_len(nOrders)){
    tmp.data <- hcData
    if(nOrders > 1){
      # A permutation drawn from argSeed alone, reproducible across runs.
      set.seed(argSeed + iOrder)
      tmp.data <- hcData[, sample(ncol(hcData)), drop = FALSE]
    }
    tmp.net <- bnlearn::hc( tmp.data, score = argScore, restart = argRestart,
                            debug = FALSE, max.iter = maxIter )
    # bnlearn::arcs() returns the learned arcs by name, so single-level columns
    # dropped above cannot shift the mapping (the previous graphNEL readback
    # indexed node names positionally against the *unfiltered* column list).
    tmp.arcs <- bnlearn::arcs(tmp.net)
    for(iArc in seq_len(nrow(tmp.arcs))){
      bump(paste(tmp.arcs[iArc, 1], tmp.arcs[iArc, 2], sep = "\r"))
    }
  }
  tmp.keys <- ls(arcTally)
  tmp.keep <- character(0)
  if(length(tmp.keys) > 0){
    tmp.freq <- vapply(tmp.keys, function(k) get(k, envir = arcTally), integer(1)) / nOrders
    tmp.keep <- tmp.keys[tmp.freq + 1e-12 >= min(max(argOrderFreq, 0), 1)]
  }
  consensusArcs <- if(length(tmp.keep) == 0) {
    matrix(character(0), ncol = 2)
  } else {
    do.call(rbind, strsplit(tmp.keep, "\r", fixed = TRUE))
  }

  #### Save the edges as an edge table
  #### ----
  allProperties = colnames( inputData.df )

  #### Compute the list of edges and their key
  allPairs <- combn( colnames( inputData.df ), 2)
  allKeys <- paste( allPairs[1,], allPairs[2,], sep = "_<<_"  )

  #### Prepare a data frame to insert the information
  myEdgesN = ncol(allPairs)
  edgesInfo.df <- data.frame( x = character(myEdgesN), y = character(myEdgesN)
                              , epresenceScore = rep(NA, myEdgesN), epresence = rep(0,myEdgesN)
                              , eorientScore = rep(NA, myEdgesN), eorient = rep(0, myEdgesN)
                              , ecorr = rep(NA,myEdgesN)
                              , stringsAsFactors = FALSE )

  edgesInfo.df[, "x"] = allPairs[1,]
  edgesInfo.df[, "y"] = allPairs[2,]

  #### Write the consensus arcs into the edge table. The pair keys keep the
  #### original column order, which is what scaleNet.gatherSubGraphs matches on.
  for( iArc in seq_len(nrow(consensusArcs)) )
  {
    myX <- consensusArcs[iArc, 1]; myY <- consensusArcs[iArc, 2]; myOrt <- 2

    tmp.idx <- which(edgesInfo.df[, "x"] == myX & edgesInfo.df[, "y"] == myY)
    if(length(tmp.idx) == 0){
      myOrt <- (-2)
      tmp.idx <- which(edgesInfo.df[, "y"] == myX & edgesInfo.df[, "x"] == myY)
    }
    if(length(tmp.idx) == 0){stop("Edge not found ?!")}

    edgesInfo.df[tmp.idx, "epresence"] = 1

    if(edgesInfo.df[tmp.idx, "eorient"] == 0){
      edgesInfo.df[tmp.idx, "eorient"] <- myOrt
    } else if(sign(myOrt*edgesInfo.df[tmp.idx, "eorient"]) == -1){
      # Both directions survived the consensus: report it as unoriented.
      edgesInfo.df[tmp.idx, "eorient"] <- 1
    }
  }
  rownames(edgesInfo.df) = allKeys

  #### Write edges
  if(!dir.exists(argOutDir)){dir.create(argOutDir)}
  inferredEdgesFileName = "edgesList.bayes_hc.txt"
  write.table(edgesInfo.df, file = file.path(argOutDir, inferredEdgesFileName), col.names = TRUE, row.names = TRUE, quote = FALSE,sep = "\t")

  if(argVerbose){cat( "# --------\n# -> END Bayes HC...\n" )}

}
