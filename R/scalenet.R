#' scalenet
#'
#' scalenet reconstruct network from local prediction based on subset of variables selected
#' with different clustering methods. The default subset selection correspond to the core
#' scaleNet approach. The other subset selection methods are left for comparison.
#'
#' @param argInData input data file path (a tab separated file) or a data.frame, with features in columns and observations in rows
#' @param argOutDir output directory file, will be created
#' @param argEigenPerc percentage of eigen vectors to use,
#' if -1 this percentage is set with an internal elbow-like heuristic method
#' @param argSubsetType variable subset selection, subset of variables can be selected from:
#' (1) eigen vector elements ("spectral"), (2) spectral k-means clusters ("spectralKmeans"),
#' (3) spectral fuzzy c-means clusters ("spectralFuzzyCmeansOrder" or "spectralFuzzyCmeansSample")
#' (4) spectral bipartition clusters ("spectralBipartition" )
#' @param argReconsMeth reconstruction method, (1) bayesian hill climbing ("bayes_hc"),
#' (2) aracne ("aracne")
#' @param argReconsParam reconstruction method parameters, (1) [bayes_hc] "-s score -r restartNumber",
#' (2) [aracne] "-p epsilon -w estimator"
#' @param argNbSamples number of observation to sample from the original dataset, if more than the
#' original number of observations, the maximum number of available samples is considered
#' @param argNumSeed random seed to sample observations from the original dataset
#' @param argNbCPU number of CPU to be used for local recontructions
#' @param argVerbose display details
#'
#' @export

# Error code from 20,000 to 20,099

# scalenet(argInData = "~/Projects/Projects_largeScale/data/benchmark/andes/input/rawData/andes_20000/andes_20000_0001.txt",
#          argOutDir = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck_spectral",
#          argVarPerc = 0.2, argReconsMeth = "aracne", argReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001)),
#          argSubsetType = "spectral", argPresFreqThresh = c(0.3, 0.8), argVerbose = TRUE)

# scalenet(argInData = "~/Projects/Projects_largeScale/data/benchmark/andes/input/rawData/andes_20000/andes_20000_0001.txt",
#          argOutDir = "~/Projects/Projects_largeScale/data/benchmark/andes/output/test_pck_spectral",
#          argVarPerc = 0.2, argReconsMeth = "bayes_hc", argReconsParam = list(bayes_hc = list(score="bde", restart=21)),
#          argSubsetType = "spectral", argPresFreqThresh = c(0.3, 0.8), argVerbose = TRUE)

# scalenet(argInData = "~/Projects/Projects_largeScale/package/ScaleNet_V1.2/tests/pop2mat.txt",
#          argOutDir = "~/Projects/Projects_largeScale/package/ScaleNet_V1.2/tests/pop2mat",
#          argVarPerc = 0.2, argReconsMeth = "bayes_hc", argReconsParam = list(bayes_hc = list(score="bde", restart=21)),
#          argSubsetType = "spectral", argPresFreqThresh = c(0.3, 0.8), argVerbose = TRUE)

scalenet <- function(argInData, argOutDir, argEigenPerc = -1, argVarPerc, argSubsetType = "spectral",
                     argReconsMeth = "bayes_hc",
                     argReconsParam = list(bayes_hc = list(score="bde", restart=20)),
                     argPresFreqThresh = c(0.5, 1), argNbSamples = 50000,
                     argNumSeed = 6196, argNbCPU = 2, argDiscretize = FALSE, argVerbose = FALSE) {

  # ----------------------------------------
  # Initialize the global variables (subEnv variable)
  # ----------------------------------------
  # An environment variable to avoid passing big table to subfunctions
  subEnv <- new.env()
#   # --> Path to directory that containes reconstruction methods
#   subEnv$reconsMeth.dirPath <- file.path(scaleNet.project.dirPath, "analyses", "reconsMeth")
  # Load/Set variables
  setEnvironment( inData = argInData, outDir = argOutDir, eigenPerc = argEigenPerc,
                           varPerc = argVarPerc, subsetType = argSubsetType,
                           reconsMeth = argReconsMeth, reconsParam = argReconsParam,
                           presFreqThresh = argPresFreqThresh,
                           nbSamples = argNbSamples, numSeed = argNumSeed, nbCPU = argNbCPU,
                           disc = argDiscretize, verbose = argVerbose, ioSubEnv = subEnv )
  # -----------------
  #  Recall the main parameters
  # -----------------
  recallParameters(ioSubEnv = subEnv)

  # ----------------------------------------
  # Compute the LAPLACIAN MATRIX
  # --> [decompose/mat] Lrw=D^-1.L=I-D^-1.W,
  # --> where D=diag(di), di=sum(Coli),
  # --> [affinity.mat] W = MIM
  # ----------------------------------------
  cat("\n# Compute the normalized Laplacian, based on MIM...\n")
  LaplacianRW(ioSubEnv = subEnv)

  # ----------------------------------------
  # Compute the eigenvectors/values
  # ----------------------------------------
  cat("\n# Compute the eigen vectors / values...\n")
  computeEigenVectVal(ioSubEnv = subEnv)

  # ----------------------------------------
  # Plot/Save the eigen vector element distribution
  # ----------------------------------------
  cat("\n# Save the ordered eigen values...\n")
  plotSaveEigenVal(ioSubEnv = subEnv)

  if(subEnv$subset.select %in% c("spectral", "spectralKmeans", "spectralFuzzyCmeansOrder", "spectralFuzzyCmeansSample")){

    # Plot/save the ordered elements of the eigen vectors
    cat("\n# Save the ordered eigen vector element...\n")
    plotSaveEigenVect(ioSubEnv = subEnv )
  }

  # ----------------------------------------
  # Reconstruct all the subgraphs
  # --> ie. two subgraphs per eigen vector
  # ----| G+,m: the m variables corresponding to the m highest eigen vector elements
  # ----| G-,m: the m variables corresponding to the m lowest eigen vector elements
  # ----------------------------------------
  cat("\n# Reconstruct all the G+,m / G-,m subgraphs...\n")
  # --> sub input datasets
  createSubInputData(ioSubEnv = subEnv)
  # --> output subgraphs
  createSubOutputGraphs(ioSubEnv = subEnv)

  # ----------------------------------------
  # Gather the subgraphs
  # ----------------------------------------
  cat("\n# Gather all the G+,m / G-,m or Gclust,m subgraphs...\n")

  # In case of spectralFuzzyCmeansSample, cmeans.nbVarSample globalNet are generated
  tmp.nbVarSamples <- ifelse((subEnv$cmeans.nbVarSample != -1), subEnv$cmeans.nbVarSample, 1)

  # iVarSpl = 1
  for(iVarSpl in seq_len(tmp.nbVarSamples)){
    scaleNet.gatherSubGraphs(ioSubEnv = subEnv, iVarSpl = iVarSpl)
  }

  # ----------------------------------------
  # Convert into a scalenet edgeList/adjMat
  # ----------------------------------------
  cat("\n# Convert into a scaleNet edgeList...\n")

  # In case of spectralFuzzyCmeansSample, cmeans.nbVarSample globalNet are generated
  tmp.nbVarSamples <- ifelse((subEnv$cmeans.nbVarSample != -1), subEnv$cmeans.nbVarSample, 1)

  # iVarSpl = 1
  for(iVarSpl in seq_len(tmp.nbVarSamples)){
    scaleNet.convertToScaleNetFormat(ioSubEnv = subEnv, iVarSpl = iVarSpl)
  }
}
