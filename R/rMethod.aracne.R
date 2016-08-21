#' rMethod.aracne
#'
#' rMethod.aracne reconstruct network using the aracne approach (from minet package)
#'
#' @param argInData input data file, a tab separated file with features in columns and observations in rows
#' @param argOutDir output directory file, will be created
#' @param argEstimator mutual information estimator
#' @param argEpsilon constraint for third edge removal
#' @param argVerbose display details
#'
#' @export

# rMethod.aracne(argInData = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck_spectral/subGraphs/subInput/eigenVector2_45neg.tsv",
#                 argOutDir = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck_spectral/subGraphs/subOutput",
#                 argEstimator = "mi.mm", argEpsilon = 0, argVerbose = TRUE)

rMethod.aracne <- function(argInData, argOutDir, argEstimator = "mi.mm", argEpsilon = 0, argVerbose = FALSE) {

  if(argVerbose){cat( "# --------\n# -> START aracne ...\n" )}

  #### Recall the main parameters
  if(argVerbose){
    cat( "# --------\n# Inputs:\n# ----\n"
         , "# Input data file -->", argInData, "\n"
         , "# Output directory -->", argOutDir, "\n" )
  }
  # ----

  #### Load the raw data file & Make sure all columns are factors, then convert to numeric
  inputData.df <- read.table(file = argInData, header = TRUE, stringsAsFactor = FALSE, sep = "\t", check.names = FALSE)
  inputData.df[, colnames(inputData.df)] <- as.data.frame(lapply(inputData.df[, colnames(inputData.df)] , factor))

  # if the estimator is not "mi.", the input data should be numeric
  if( length( grep( "mi.", argEstimator ) ) == 0 )
  { inputData.df[, colnames(inputData.df)] <- as.data.frame(lapply(inputData.df[, colnames(inputData.df)] , as.numeric)) }

  #### Recall the main parameters
  if(argVerbose){
    cat( "# All properties -->", paste( colnames(inputData.df), collapse = ", " ), "\n"
         , "# MI estimator type -->", argEstimator, "\n"
         , "# Eps -->", argEpsilon, "\n"
         , "# --------\n" )
  }
  # ----

  #### Compute the MIM
  print("HERE")
  myMIM <- minet::build.mim( inputData.df, argEstimator )

  # Compute the weighted matrix
  myNet = minet::aracne( myMIM, eps = argEpsilon )

  #### Save the edges as an edges table and an adj mat
  #### ----
  gV <- new.env()
  gV$allProperties <- colnames(inputData.df)

  #### Compute the list of edges and their key
  allPairs <- combn( colnames( inputData.df ), 2)
  allKeys <- paste( allPairs[1,], allPairs[2,], sep = "_<<_"  )
  rm(gV)

  #### Insert the information into a data frame
  myEdgesN = ncol(allPairs)
  edgesInfo.df <- data.frame( x = character(myEdgesN), y = character(myEdgesN)
                              , epresenceScore = rep(NA, myEdgesN), epresence = rep(0,myEdgesN)
                              , eorientScore = rep(NA, myEdgesN), eorient = rep(0, myEdgesN)
                              , ecorr = rep(NA,myEdgesN)
                              , stringsAsFactors = FALSE )

  edgesInfo.df[, "x"] = allPairs[1,]
  edgesInfo.df[, "y"] = allPairs[2,]
  for( iEdge in seq_len(myEdgesN) )
  {
    #### Copy from the mim mat
    edgesInfo.df[iEdge, "epresenceScore"] = myNet[edgesInfo.df[iEdge, "x"],edgesInfo.df[iEdge, "y"]]
  }
  rownames(edgesInfo.df) = allKeys

  #### Order the data frame
  edgesInfo.df = edgesInfo.df[order(edgesInfo.df[, "epresenceScore"], decreasing = TRUE),]

  #### Init the edges as phantom (0) and if the weight is > 0, set them as non phantom (1)
  edgesInfo.df[, "epresence"] = 0
  myInfEdges = which( edgesInfo.df[, "epresenceScore"] > 0 )
  if( length( myInfEdges ) > 0 ){
    edgesInfo.df[myInfEdges, "epresence"] = 1
    edgesInfo.df[myInfEdges, "eorient"] = 1
  }

  #### Write the mutual information
  if(!dir.exists(argOutDir)){dir.create(argOutDir)}
  inferredEdgesFileName = "edgesList.aracne.txt"
  write.table( edgesInfo.df, file = file.path(argOutDir, inferredEdgesFileName), col.names = TRUE, row.names = TRUE, quote = FALSE,sep = "\t")

  if(argVerbose){cat( "# --------\n# -> END aracne...\n" )}

}
