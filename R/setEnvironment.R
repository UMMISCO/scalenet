#' setEnvironment
#'
#' setEnvironment sets all required information in en env() variables for reconstruction
#'
#' @param inData input data file path (a tab separated file) or a data.frame, with features in columns and observations in rows
#' @param outDir output directory file, will be created
#' @param eigenPerc percentage of eigen vectors to use,
#' if -1 this percentage is set with an internal elbow-like heuristic method
#' @param subsetType variable subset selection, subset of variables can be selected from:
#' (1) eigen vector elements ("spectral"), (2) spectral k-means clusters ("spectralKmeans"),
#' (3) spectral fuzzy c-means clusters ("spectralFuzzyCmeansOrder" or "spectralFuzzyCmeansSample")
#' (4) spectral bipartition clusters ("spectralBipartition" )
#' @param reconsMeth reconstruction method, (1) bayesian hill climbing ("bayes_hc"),
#' (2) aracne ("aracne")
#' @param reconsParam reconstruction method parameters, (1) [bayes_hc] "-s score -r restartNumber",
#' (2) [aracne] "-p epsilon -w estimator"
#' @param nbSamples number of observation to sample from the original dataset, if more than the
#' original number of observations, the maximum number of available samples is considered
#' @param numSeed random seed to sample observations from the original dataset
#' @param nbCPU number of CPU to be used for local recontructions
#' @param verbose display details
#' @param ioSubEnv a global environment variable

# Error code from 1,000 to 1,099

setEnvironment <- function(inData, outDir, eigenPerc, varPerc, subsetType,
                           reconsMeth, reconsParam, presFreqThresh, nbSamples,
                           numSeed, nbCPU, disc, verbose, ioSubEnv){

  # Check some input arguments
  if(is.null(inData)){stop("# --Err-- 1000")}
  if(is.null(outDir)){stop("# --Err-- 1001")}
  if(is.null(eigenPerc)){eigenPerc <- -1}
  if(is.null(varPerc)){stop("# --Err-- 1002")}
  if(is.null(subsetType)){subsetType <- "spectral"}
  if(is.null(reconsMeth)){stop("# --Err-- 1003")}
  if(is.null(reconsParam)){stop("# --Err-- 1004")}
  if(is.null(presFreqThresh)){presFreqThresh <- c(0.5, 1)}
  if(is.null(nbSamples)){nbSamples <- 50000}
  if(is.null(numSeed)){numSeed <- 6196}
  if(is.null(nbCPU)){nbCPU <- 2}
  if(is.null(verbose)){verbose <- FALSE}
  if(is.null(ioSubEnv)){stop("# --Err-- 1005")}

  ioSubEnv$verbose <- verbose
  ioSubEnv$total.time <- NULL

  # --> Mutual information estimator
  # ----| The affinity matrix is based on the mutual information values
  # ----| mi.emp, mi.shrink, mi.mm, spearman, kendall, pearson
  ioSubEnv$mi.estimator <- "mi.mm" # Miller-Madow asymptotic bias corrected empirical estimator
  # --> Type of matrix to decompose
  ioSubEnv$similarityType <- "Lrw"

  # --> input data file / output directory path / input data
  ioSubEnv$inputData.filePath <- NULL
  ioSubEnv$allData <- NULL
  if(class(inData) == "character"){

    ioSubEnv$inputData.filePath <- inData
    if(!file.exists(ioSubEnv$inputData.filePath)){stop("# --Err-- 1007")}
    cat("# Loading data:", basename(ioSubEnv$inputData.filePath), "...\n")
    ioSubEnv$allData <- suppressWarnings(data.table::fread(input = ioSubEnv$inputData.filePath, header = T,
                                                           sep = "\t", stringsAsFactors = F, data.table = F,
                                                           showProgress = F))
  } else if(class(inData) == "data.frame") {

    cat("# Loading data...\n")
    ioSubEnv$allData <- inData

  } else { stop("# --Err-- 1006") }

  # Discretize the data is required
  if(disc == TRUE){
    if(class(inData) == "character"){
      ioSubEnv$allData <- discretize(argInData = ioSubEnv$inputData.filePath, argMaxClusters = 5, argVerbose = verbose)
    } else if(class(inData) == "data.frame"){
      ioSubEnv$allData <- discretize(argInData = ioSubEnv$allData, argMaxClusters = 5, argVerbose = verbose)
    }
  }

  #### Make sure all columns are factors, then convert to numeric
  if(!ioSubEnv$mi.estimator %in% c("pearson", "spearman", "kendall")){
    # print(head(ioSubEnv$allData))
    ioSubEnv$allData[, colnames(ioSubEnv$allData)] <- as.data.frame(lapply(ioSubEnv$allData[, colnames(ioSubEnv$allData)] , factor))
  }
  ioSubEnv$allData[, colnames(ioSubEnv$allData)] <- as.data.frame(lapply(ioSubEnv$allData[, colnames(ioSubEnv$allData)] , as.numeric))
  ioSubEnv$allData <- as.matrix(ioSubEnv$allData)

  if(ioSubEnv$verbose){
    cat("# --|", paste(head(colnames(ioSubEnv$allData)), collapse = ','), "...\n")
    cat("# --| nbVar:", dim(ioSubEnv$allData)[2], ", nbSamples:", dim(ioSubEnv$allData)[1],"\n")

    cat("# --| Head of sampled data:\n")
    cat("# -----------------------------\n")
    print(ioSubEnv$allData[1:5,1:5])
    print(dim(ioSubEnv$allData))
    cat("# -----------------------------\n")
  }

  # --> Create the output directory
  ioSubEnv$output.dirPath <- outDir
  if(dir.exists(ioSubEnv$output.dirPath)){
    if(ioSubEnv$verbose){
      print(ioSubEnv$output.dirPath)
      print("# --Wlib3-- Output directory path already exists!")
    }
  }else{dir.create(ioSubEnv$output.dirPath)}

  # --> Percent of vertices per subgraph / of eigen vectors
  ioSubEnv$subset.k.perc <- eigenPerc
  if(ioSubEnv$subset.k.perc > 1){stop("# --Eslib4-- Precent of eigen vectors should be less than 1!")}
  ioSubEnv$subset.m.perc <- varPerc
  if(ioSubEnv$subset.m.perc > 1){stop("# --Eslib5-- Precent of variables per subgraph should be less than 1!")}

  # --> How the subset variables are choosen (spectral, random, spectralKmeans, spectralFuzzyCmeansOrder)
  ioSubEnv$subset.select <- subsetType
  if( !( ioSubEnv$subset.select %in% c( "random", "spectral", "spectralKmeans",
                                        "spectralFuzzyCmeansOrder", "spectralFuzzyCmeansSample",
                                        "spectralBipartition" ) ) )
  { stop( "# --Eslib6-- Unkown subset selection: ", ioSubEnv$subset.select ) }

  # Number of CPU to be used
  ioSubEnv$nbCPU <- nbCPU

  # --> Set the reconstruction method and associated parameters
  ioSubEnv$recons.method <- reconsMeth
  if( !( ioSubEnv$recons.method %in% c( "aracne", "bayes_hc") ) )
  { stop( "# --Eslib7-- Unkown method: ", ioSubEnv$recons.method ) }
  ioSubEnv$recons.param <- reconsParam
  ioSubEnv$recons.script <- paste(ioSubEnv$recons.method, "R", sep = ".")

#   # --> File path to save stats summary
#   ioSubEnv$statSummary.filePath <- ifelse(nchar(argList[["argStatSummary"]]) == 0,
#                                           file.path(ioSubEnv$output.dirPath, paste(basename(ioSubEnv$output.dirPath),
#                                                                                    "statSummary.tsv", sep='_')),
#                                           argList[["argStatSummary"]])
  ioSubEnv$numSeed <- numSeed
  ioSubEnv$nbSamples <- nbSamples
  ioSubEnv$presFreqThresh <- presFreqThresh

  # Complementary information to keep from the subnet depending on the reconstruction method
  # ioSubEnv$recons.method.vect <- c("aracne")
  # ioSubEnv$net.cmplInfo <- vector("list", length(ioSubEnv$recons.method.vect))
  # ioSubEnv$net.cmplInfo[["aracne"]] <- c("info")

  # --> subset.m: nbr variables per subgraph (by default, set mean, median and sd of subset sizes)
  ioSubEnv$subset.m <- ceiling(ioSubEnv$subset.m.perc*dim(ioSubEnv$allData)[2])
  ioSubEnv$subset.m.mean <- ioSubEnv$subset.m
  ioSubEnv$subset.m.med <- ioSubEnv$subset.m
  ioSubEnv$subset.m.sd <- 0
  # --> subset.k: nbr of eigen vector
  if(ioSubEnv$subset.k.perc != -1){
    ioSubEnv$subset.k <- ceiling(ioSubEnv$subset.k.perc*dim(ioSubEnv$allData)[2])}

  # --> kmeans.k: nbr of clusters (for 'spectralKmeans', 'spectralFuzzyCmeansOrder')
  # --> cluster.asw: cluster average silhouette width (for 'spectralKmeans', given by pamk)
  ioSubEnv$kmeans.k = -1
  ioSubEnv$cmeans.k = -1
  ioSubEnv$cluster.asw = -1
  ioSubEnv$cmeans.nbVarSample = -1

  # if( ioSubEnv$subset.select %in% c("spectralKmeans", "spectralFuzzyCmeansOrder") ){
  if( ioSubEnv$subset.select == "spectralKmeans" ){

    ioSubEnv$kmeans.k <- floor(1/ioSubEnv$subset.m.perc)
    if(ioSubEnv$kmeans.k < 2){stop("# --Eslib8-- Two cluster at least are expected!")}

  } else if( ioSubEnv$subset.select == "spectralFuzzyCmeansOrder" ){

    if(ioSubEnv$subset.k.perc != -1){
      ioSubEnv$cmeans.k <- 2*ioSubEnv$subset.k
      if(ioSubEnv$cmeans.k < 2){stop("# --Eslib10-- Two cluster at least are expected!")}
    }

  } else if( ioSubEnv$subset.select == "spectralFuzzyCmeansSample" ){

    ioSubEnv$cmeans.nbVarSample = 20
    ioSubEnv$cmeans.k <- floor(1/ioSubEnv$subset.m.perc)
    if(ioSubEnv$cmeans.k < 2){stop("# --Eslib11-- Two cluster at least are expected!")}
  }

  # Set a vector for the speactral bipartitionning
  ioSubEnv$specBi.allClusters <- c()  # Should then have name of variables as names, and nbr clusters as values
  ioSubEnv$specBi.nextCluster <- -1

#   # Create the stats summary file if needed
#   if( scaleNet.createStatSummary(inFilePath = ioSubEnv$statSummary.filePath) != 0 ){
#     warning("# --Wlib1--> Appending to existing stats summary file already...")}

}
