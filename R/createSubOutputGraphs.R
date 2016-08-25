#' createSubInputData
#'
#'
#' @param ioSubEnv a global environment variable

createSubOutputGraphs <- function( ioSubEnv ){

    # Set the directory that will contain all the reconstructed subgraphs
  tmp.dirPath <- file.path(ioSubEnv$output.dirPath, "subGraphs")
  if(!dir.exists(tmp.dirPath)){dir.create(tmp.dirPath)}
  if(!dir.exists(file.path(tmp.dirPath, "subOutput"))){dir.create(file.path(tmp.dirPath, "subOutput"))}

  # --!-- Time
  time.multiReconstruction.output.start <- proc.time()

  if(ioSubEnv$subset.select %in% c("spectral", "random")){

    # Create all command lines
    # --> define the input sub data set file template
    subInData.filePath.template = file.path(tmp.dirPath, "subInput", paste("eigenVector", "XXX_IEIGEN_XXX", "_", ioSubEnv$subset.m,"XXX_SIGN_XXX.tsv", sep=''))
    subOutData.filePath.template = gsub( "subInput", "subOutput", gsub(".tsv", "", subInData.filePath.template))

    # Build a call template
    cmd.template <- ""
    if(ioSubEnv$recons.method == "bayes_hc"){
      cmd.template <- paste( "rMethod.hc(argInData = '", subInData.filePath.template, "', argOutDir = '", subOutData.filePath.template,
                             "', argScore = '", ioSubEnv$recons.param[["score"]], "', argRestart = ",
                             ioSubEnv$recons.param[["restart"]], ", argVerbose = ", as.character(ioSubEnv$verbose), ")", sep = "")

    } else if(ioSubEnv$recons.method == "aracne"){
      cmd.template <- paste( "rMethod.aracne(argInData = '", subInData.filePath.template, "', argOutDir = '", subOutData.filePath.template,
                             "', argEstimator = '", ioSubEnv$recons.param[["estimator"]], "', argEpsilon = ",
                             ioSubEnv$recons.param[["epsilon"]], ", argVerbose = ", as.character(ioSubEnv$verbose), ")", sep = "")
    }

    # --> set with the eigen value number and sign
    cmd.template.split <- unlist(strsplit(cmd.template, split = "XXX_IEIGEN_XXX"))
    cmd.iEigen.paste <- paste(cmd.template.split[1], seq_len(ncol(ioSubEnv$eigen.vectors))[-1],
                              cmd.template.split[2], seq_len(ncol(ioSubEnv$eigen.vectors))[-1],
                              cmd.template.split[3], sep='')

    cmd.all <- c(gsub("XXX_SIGN_XXX", "pos", cmd.iEigen.paste), gsub("XXX_SIGN_XXX", "neg", cmd.iEigen.paste))

    # Call the command in parallel
#     doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
#     tmpOut <- foreach(i=1:length(cmd.all)) %dopar%{ eval(parse(text = cmd.all[i])) }
    for(i in 1:length(cmd.all)) { eval(parse(text = cmd.all[i])) }

  } else if(ioSubEnv$subset.select %in% c("spectralKmeans", "spectralFuzzyCmeansOrder",
                                          "spectralFuzzyCmeansSample", "spectralBipartition")){

    # Create all command lines
    # --> define the input sub data set file template
    subInData.filePath.template <- file.path(tmp.dirPath, "subInput", paste(ioSubEnv$subset.select, "XXX_ICLUST_XXX", "_", ioSubEnv$subset.m,".tsv", sep=''))
    subOutData.filePath.template <- gsub( "subInput", "subOutput", gsub(".tsv", "", subInData.filePath.template))

    # Build a call template
    cmd.template <- ""
    if(ioSubEnv$recons.method == "bayes_hc"){
      cmd.template <- paste( "rMethod.hc(argInData = '", subInData.filePath.template, "', argOutDir = '", subOutData.filePath.template, "', argScore = '", ioSubEnv$recons.param[["bayes_hc"]][["score"]], "', argRestart = ", ioSubEnv$recons.param[["bayes_hc"]][["restart"]], ", argVerbose = ", as.character(ioSubEnv$verbose), ")", sep = "")

    } else if(ioSubEnv$recons.method == "aracne"){
      cmd.template <- paste( "rMethod.aracne(argInData = '", subInData.filePath.template, "', argOutDir = '", subOutData.filePath.template, "', argEstimator = '", ioSubEnv$recons.param[["aracne"]][["estimator"]], "', argEpsilon = ", ioSubEnv$recons.param[["aracne"]][["epsilon"]], ", argVerbose = ", as.character(ioSubEnv$verbose), ")", sep = "")
    }

#     cmd.template <- paste( "cd ", file.path(ioSubEnv$reconsMeth.dirPath, ioSubEnv$recons.method), ";", " ./", ioSubEnv$recons.script, " -i ",
#                            subInData.filePath.template, " -o ", subOutData.filePath.template, " ", ioSubEnv$recons.param, sep = "")
    # --> set with the cluster number
    cmd.template.split <- unlist(strsplit(cmd.template, split = "XXX_ICLUST_XXX"))

    cmd.iClust.paste <- ""
    if(ioSubEnv$subset.select == "spectralKmeans"){

      cmd.iClust.paste <- paste(cmd.template.split[1], seq_len(ioSubEnv$kmeans.k),
                                cmd.template.split[2], seq_len(ioSubEnv$kmeans.k),
                                cmd.template.split[3], sep='')

    } else if(ioSubEnv$subset.select == "spectralFuzzyCmeansOrder"){

      cmd.iClust.paste <- paste(cmd.template.split[1], seq_len(ioSubEnv$cmeans.k),
                                cmd.template.split[2], seq_len(ioSubEnv$cmeans.k),
                                cmd.template.split[3], sep='')

    } else if(ioSubEnv$subset.select == "spectralBipartition"){

      cmd.iClust.paste <- paste(cmd.template.split[1], unique(ioSubEnv$specBi.allClusters),
                                cmd.template.split[2], unique(ioSubEnv$specBi.allClusters),
                                cmd.template.split[3], sep='')

    } else if(ioSubEnv$subset.select == "spectralFuzzyCmeansSample"){

      # Create all command lines
      # --> define the input sub data set file template
      subInData.filePath.template <- file.path(tmp.dirPath, "subInput",
                                               paste(ioSubEnv$subset.select, "XXX_ICLUST_XXX", "_", ioSubEnv$subset.m,"_XXX_IVARSPL_XXX.tsv", sep=''))
      subOutData.filePath.template <- gsub( "subInput", "subOutput", gsub(".tsv", "", subInData.filePath.template))

      # Build a call template
      cmd.template <- ""
      if(ioSubEnv$recons.method == "bayes_hc"){
        cmd.template <- paste( "try(rMethod.hc(argInData = '", subInData.filePath.template, "', argOutDir = '", subOutData.filePath.template, "', argScore = '", ioSubEnv$recons.param[["bayes_hc"]][["score"]], "', argRestart = ", ioSubEnv$recons.param[["bayes_hc"]][["restart"]], ", argVerbose = ", as.character(ioSubEnv$verbose), "))", sep = "")

      } else if(ioSubEnv$recons.method == "aracne"){
        cmd.template <- paste( "try(rMethod.aracne(argInData = '", subInData.filePath.template, "', argOutDir = '", subOutData.filePath.template, "', argEstimator = '", ioSubEnv$recons.param[["aracne"]][["estimator"]], "', argEpsilon = ", ioSubEnv$recons.param[["aracne"]][["epsilon"]], ", argVerbose = ", as.character(ioSubEnv$verbose), "))", sep = "")
      }
#       cmd.template <- paste( "cd ", file.path(ioSubEnv$reconsMeth.dirPath, ioSubEnv$recons.method), ";", " ./", ioSubEnv$recons.script, " -i ",
#                              subInData.filePath.template, " -o ", subOutData.filePath.template, " ", ioSubEnv$recons.param, sep = "")

      cmd.iClust.paste <- c()
      # iNbVarSpl = 1
      for(iNbVarSpl in seq_len(ioSubEnv$cmeans.nbVarSample)){

        # Set the number of the assignment
        tmp.cmd.template <- gsub("XXX_IVARSPL_XXX", iNbVarSpl, cmd.template)

        # --> set with the cluster number
        tmp.cmd.template.split <- unlist(strsplit(tmp.cmd.template, split = "XXX_ICLUST_XXX"))

        cmd.iClust.paste <- c( cmd.iClust.paste,
                               paste(tmp.cmd.template.split[1], seq_len(ioSubEnv$cmeans.k),
                                     tmp.cmd.template.split[2], seq_len(ioSubEnv$cmeans.k),
                                     tmp.cmd.template.split[3], sep=''))
      }
    }

    cmd.all <- cmd.iClust.paste

    # Call the command in parallel
    # i = 1
    doParallel::registerDoParallel(cores = ioSubEnv$nbCPU)
    tmpOut <- foreach(i=1:length(cmd.all)) %dopar%{
      print(cmd.all[i])
      eval(parse(text = cmd.all[i]))
    }

#     for(i in c(1:length(cmd.all))){
#       print(cmd.all[i])
#       eval(parse(text = cmd.all[i]))
#     }

  } else { stop("# Not implemented yet...") }

  # --!-- Time
  ioSubEnv$outputs.time <- (proc.time() - time.multiReconstruction.output.start)
  ioSubEnv$total.time <- ioSubEnv$total.time + ioSubEnv$outputs.time
}
