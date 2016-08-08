#' rMethod.hc
#'
#' rMethod.hc reconstruct network using the hhill climbing approach (from bnlearn package)
#'
#' @param argInData input data file, a tab separated file with features in columns and observations in rows
#' @param argOutDir output directory file, will be created
#' @param argScore scoring function
#' @param argRestart search and score restart
#' @param argVerbose display details
#'
#' @export

# rMethod.hc(argInData = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck/subGraphs/subInput/eigenVector2_45neg.tsv",
#          argOutDir = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck/subGraphs/subOutput",
#          argScore = "bde", argRestart = 20, argVerbose = TRUE)

rMethod.hc <- function(argInData, argOutDir, argScore = "bde", argRestart = 20, argVerbose = FALSE) {

  if(argVerbose){cat( "# --------\n# -> START Bayes HC...\n" )}

  if(argVerbose){
    #### Recall the main parameters
    cat( "# --------\n# Inputs:\n# ----\n"
         , "# Input data file -->", argInData, "\n"
         , "# Output directory -->", argOutDir, "\n" )
    # ----
  }

  #### Load the raw data file
  inputData.df <- read.table(file = argInData, header = TRUE, stringsAsFactor = FALSE, sep = "\t")
  if(ncol(inputData.df)<=2){quit(save = "no", status = 0)}

  #### Remove NA values if any
  naRows <- unique(which(is.na(inputData.df), arr.ind = T)[,1])
  if( length(naRows) > 0 ){ inputData.df <- inputData.df[-naRows,] }

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

  myNet = NULL
  #### hc
  if( length( varWithSingleLevel ) > 0 )
  {
    myNet = bnlearn::hc( inputData.df[, -varWithSingleLevel], score = argScore,
                         restart =  argRestart, debug = FALSE )

    if(argVerbose){cat("# ------| Properties with less than 2 levels: ", paste( tmp.allProperties[varWithSingleLevel], collapse = "," ), "\n" )}
    tmp.allProperties = tmp.allProperties[-varWithSingleLevel]
  } else {
    myNet = bnlearn::hc( inputData.df, score = argScore, restart =  argRestart, debug = FALSE )
  }
  rm(tmp.allProperties)

  #### Save the edges as an edge table
  #### ----
  allProperties = colnames( inputData.df )

  #### Compute the list of edges and their key
  myNet.graphNEL = bnlearn::as.graphNEL(myNet)
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
  # iNode <- 1
  for( iNode in seq_len(length(graph::edges(myNet.graphNEL))) )
  {
    myNode = colnames( inputData.df )[iNode]
    myPartner.vect = graph::edges(myNet.graphNEL)[[myNode]]

    if( length( myPartner.vect ) > 0 )
    {
      # myNode ---> myPartner
      # iPartner <- 1
      for( iPartner in seq_len(length(myPartner.vect)) )
      {
        myPartner = myPartner.vect[iPartner]

        # (1) Suppose first that the edge can be found forward
        myX <- myNode; myY <- myPartner; myOrt <- 2
        tmp.idx <- which(edgesInfo.df[, "x"] == myX & edgesInfo.df[, "y"] == myY)

        # (2) If not, set the opposite direction
        if(length(tmp.idx) == 0){
          myOrt <- (-2)
          tmp.idx <- which(edgesInfo.df[, "y"] == myX & edgesInfo.df[, "x"] == myY)
        }

        # (3) If nothing...problem?!
        if(length(tmp.idx) == 0){stop("Edge not found ?!")}

        # Set the presence
        edgesInfo.df[tmp.idx, "epresence"] = 1

        if(edgesInfo.df[tmp.idx, "eorient"] == 0){

          # If no orientation already exists, set the orientation
          edgesInfo.df[tmp.idx, "eorient"] <- myOrt

        } else if(sign(myOrt*edgesInfo.df[tmp.idx, "eorient"]) == -1){

          # If an opposite orientation exists, set 1 (should never occur...)
          edgesInfo.df[tmp.idx, "eorient"] <- 1
        }
      }
    }
  }
  rownames(edgesInfo.df) = allKeys

  #### Write edges
  if(!dir.exists(argOutDir)){dir.create(argOutDir)}
  inferredEdgesFileName = "edgesList.bayes_hc.txt"
  write.table(edgesInfo.df, file = file.path(argOutDir, inferredEdgesFileName), col.names = TRUE, row.names = TRUE, quote = FALSE,sep = "\t")

  if(argVerbose){cat( "# --------\n# -> END Bayes HC...\n" )}

}
