#' recallParameters
#'
#' recallParameters displays set parameters
#'
#' @param ioSubEnv a global environment variable

recallParameters <- function(ioSubEnv){

  cat( "\n# --------\n# Inputs:\n# ----\n"
       , "# Input data file path -->", ioSubEnv$inputData.filePath, "\n"
       , "# Output directory path -->", ioSubEnv$output.dirPath, "\n"
       , "# Nbr. Variables -->", ncol(ioSubEnv$allData), "\n"
       , "# Subset selection -->", ioSubEnv$subset.select, "\n"
       , "# Perc. of eigenvectors -->", ioSubEnv$subset.k.perc, "\n"
       , "# Perc. of variables -->", ioSubEnv$subset.m.perc, "\n"
       , "# Nbr. kmeans clusters -->", ioSubEnv$kmeans.k, "\n"
       , "# Nbr.cmeans clusters -->", ioSubEnv$cmeans.k, "\n"
       , "# Nbr. CPU -->", ioSubEnv$nbCPU, "\n# --\n"
       , "# MI estimator -->", ioSubEnv$mi.estimator, "\n"
       , "# Laplacian graph -->", ioSubEnv$similarityType, "\n# --\n"
       , "# Stat. Summary file path -->", ioSubEnv$statSummary.filePath, "\n# --\n"
       , "# reconsMeth location -->", ioSubEnv$reconsMeth.dirPath, "\n# --------\n" )

}
