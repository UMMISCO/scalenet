#
# # load some results
# load("~/Research/workspace_r/predomics/test/db5_mh_2101_gu3.9_bin/terga_db=lgc_db5__sparsity=1_to_30_population=100_convergence=10.rda")
# # get a population of models
# pop <- c(res.terga$models$k_1, res.terga$models$k_2, res.terga$models$k_3, res.terga$models$k_4, res.terga$models$k_5,
#          res.terga$models$k_6, res.terga$models$k_7, res.terga$models$k_8, res.terga$models$k_9, res.terga$models$k_10)
# predomics:::printPop(pop)
#
# # TODO export listOfModelsToDenseCoefMatrix
# pop.2.mat <- predomics:::listOfModelsToDenseCoefMatrix(clf=clf.terga, X=X, y=y, list.models=pop, rm.empty=TRUE, order.row=TRUE)
# image(pop.2.mat)
# write.table(t(pop.2.mat), file="pop2mat.txt", quote=FALSE, row.names = FALSE, sep="\t")


# Launch network reconstruction
#install.packages("~/Desktop/scalenet/scalenet_V1.2.tar.gz", repos = NULL, type = "source")
# source("https://bioconductor.org/biocLite.R")
# biocLite("minet")
library(scalenet)


# scalenet(argInData=pop2mat,
#          argOutDir="/Users/eprifti/Research/workspace_r/scalenet/tests/scalent_results",
#          argVarPerc = 0.2, argReconsMeth = "bayes_hc", argReconsParam = list(bayes_hc = list(score="bde", restart=21)),
#          argSubsetType = "spectral", argPresFreqThresh = c(0.3, 0.8), argVerbose = TRUE)
#
# scalenet(argInData=pop2mat,
#          argOutDir="/Users/eprifti/Research/workspace_r/scalenet/tests/scalent_results",
#          argVarPerc = 0.2, argReconsMeth = "aracne", argReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001)),
#          argSubsetType = "spectral", argPresFreqThresh = c(0.3, 0.8), argVerbose = TRUE)

tmp <- scs(workspaceDir = "scalent_results",
     #argInData = "/Users/eprifti/Research/akkersisters/pooled_cohort_baseline_MGS.txt",
     argInData = "pop2mat.txt",
     argReconsMeth = c("aracne", "bayes_hc"),
     argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"), bayes_hc = list(ort = "y", eweight = "ecorr")),
     argEmbReconsParam = list(aracne = list(estimator="mi.mm", epsilon=0.001), bayes_hc = list(score="bde", restart=21), varPerc = 0.2),
     argPresFreqThresh = c(0.3, 0.8), clean.workspace = FALSE, argDiscretize = TRUE, argVerbose = TRUE)

# pop2mat <- read.delim("pop2mat.txt")
# load("~/Research/akkersisters/pc.bc.dat.big.rda")
# scs(workspaceDir="~/Research/workspace_r/scalenet/tests/",
#     argInData = t(pc.bc.dat.big),
#     argReconsMeth = c("aracne", "bayes_hc"),
#     argReconsMethInfo = list(aracne = list(ort = "n", eweight = "epresenceScore"),
#                              bayes_hc = list(ort = "y", eweight = "ecorr")),
#     argEmbReconsParam = list(aracne = list(estimator = "mi.mm", epsilon = 0.001),
#                              bayes_hc = list(score = "bde", restart = 21), varPerc = 0.2),
#     argPresFreqThresh = c(0.5, 1), clean.workspace = FALSE,
#     argDiscretize = FALSE, argVerbose = TRUE)



# visualize output

#-------------------------------------------------------------------------------------------
# 2. scalenet NETWORK
#-------------------------------------------------------------------------------------------
# flipEdges <- function(var.name = c(node1, node2)){
#   if(length(var.name)!=2) stop("Two variable names are needed")
#   if(var.name[1] > var.name[2]) {
#     tmp <- var.name[1]
#     var.name[1] <- var.name[2]
#     var.name[2] <- tmp
#   }
#   return(var.name)
# }

# node annotation
load("mgs_taxo.rda"); taxo <- mgs_taxo; rm(mgs_taxo)
size <- as.numeric(as.character(taxo$size)); plot(size)


fname <- "scalent_results/consensusGraph/pop2mat_consensusNetpresFreq1/edgesList.txt"
# upper.perc <- 3
# upper.fix <- 307

# load the edge information for spectral3off2 network
edges <- read.delim(fname, as.is = TRUE); dim(edges) # 121278 edges and 7 columns
edges.raw <- edges # save it
edges <- edges[!is.na(edges$eorientScore),]; dim(edges) # 6389 edges
colnames(edges)[1:2] <- c("from","to")
edges$from <- gsub("_",":",edges$from)
edges$to <- gsub("_",":",edges$to)
rownames(edges) <- paste(edges$from, edges$to, sep=" => ")

# take upper.perc% of the links
#ind.signif <- getUpperPercIndex(v=edges$info, upper.fix = upper.fix)
#ind.signif <- getUpperPercIndex(v=edges$info, upper.perc = upper.perc)
#edges <- edges[ind.signif,]; dim(edges) # 307 edges

# # flip node names alphabetically
# edges.flipped <- edges
# flipped <- rep(FALSE, nrow(edges))
# for(i in 1:nrow(edges.flipped)){
#   res <- flipEdges(var.name = edges.flipped[i,c("from","to")])
#   if(paste(res,collapse = " => ") == paste(edges.flipped[i,c("from","to")],collapse = " => ")) flipped[i] <- TRUE
#   edges.flipped[i,c("from","to")] <- res
# }
# ind <- edges.flipped$infOrt !=1 & flipped
# edges.flipped$infOrt[ind] <- edges.flipped$infOrt[ind]*-1
# #require(plyr)
# #edges.flipped[,c("node1","node2")] <- adply(.data = edges.flipped[,c("node1", "node2")],  .margins = 1, .fun = function(x) data.frame(flipEdges(x)))
# # order nodes alphabetically
# edges.flipped <- edges.flipped[order(edges.flipped$from, edges.flipped$to),]
# rownames(edges.flipped) <- paste(edges.flipped$from, edges.flipped$to, sep=" => ")



# BUILD NETWORK
#-------------------------------------------------------------------------------------------
# ANNOTATION of the edges
edges.annot <- taxo[unique(c(edges$from,edges$to)),]; edges.annot <- data.frame(rownames(edges.annot),edges.annot);
colnames(edges.annot)[match("name",colnames(edges.annot))] <- "name_long"; colnames(edges.annot)[1] <- "name"

# Build the igraph netowrk
require(igraph)
# create the igraph object
gD <- graph.data.frame(d = edges, directed = TRUE, vertices = edges.annot)

# Calculate degree for all nodes
degAll <- igraph::degree(gD, v = V(gD), mode = "all")
gD <- set.vertex.attribute(gD, "degree", index = V(gD)$name, value = degAll)
# Calculate betweenness for all nodes
betAll <- igraph::betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll)); rm(betAll)
# Add new node/edge attributes based on the calculated node properties/similarities
gD <- set.vertex.attribute(gD, "betweenness", index = V(gD)$name, value = betAll.norm)

# Calculate edge properties and add to the network
E(gD)$infOrt[E(gD)$infOrt==1] <- 0; E(gD)$infOrt[E(gD)$infOrt==-2] <- 1
#Calculate Dice similarities between all pairs of nodes
dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
# The following function will transform a square matrix to an edge driven one and add values to each edge
F1 <- function(x) {data.frame(dice = dsAll[which(V(gD)$name == as.character(x$from)), which(V(gD)$name == as.character(x$to))])}
library(plyr)
edges.ext <- ddply(edges, .variables=c("from", "to"), function(x) data.frame(F1(x))); dim(edges.ext)

gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
E(gD)[as.character(edges.ext$from) %--% as.character(edges.ext$to)]$similarity <- as.numeric(edges.ext$dice)
# Edge color
E(gD)$color <- c("#DC143C","#A6A6A6")[as.factor(factor(sign(E(gD)$eorient), levels=c('-1','1')))]

# fix orientation coding
# 1 backward
# 2 forward
# 3 bidirected (on laisse 0)
E(gD)$eorient[E(gD)$eorient==1] <- 2
E(gD)$eorient[E(gD)$eorient==-1] <- 1


# #check things out
# node <- "akker"
# dataSet[dataSet$node1 == node,]; dataSet[dataSet$node2 == node,]
# node <- "GU:107"
# dataSet.ext[dataSet.ext$node1 == node,]; dataSet.ext[dataSet.ext$node2 == node,]

# Check the attributes
# Print number of nodes and edges
print(paste("There are",vcount(gD),"nodes and",ecount(gD),"edges"))
# [1] "There are 765 nodes and 1597 edges"
summary(gD)
# IGRAPH DN-- 443 307 --
#   + attr: name (v/c), size (v/n), NA_pc (v/n), BHit_pc (v/n), BHit (v/c), BH_ali (v/n), BH_id (v/n),
# | NA_pc_ (v/n), Taxo_pc (v/n), Taxo (v/c), Taxo_ali (v/n), Taxo_id (v/n), Taxo_level (v/c), annot
# | (v/c), species (v/c), genus (v/c), family (v/c), order (v/c), class (v/c), phylum (v/c),
# | superkingdom (v/c), name_long (v/c), color (v/c), degree (v/n), betweenness (v/n), type (e/c), ui
# | (e/l), info (e/n), cplx (e/n), Nxy_ui (e/l), confidence (e/n), infOrt (e/n), trueOrt (e/l), isOrt
# | (e/l), isOrtOk (e/l), essential (e/l), sign (e/c), diff (e/n), similarity (e/n), color (e/c)

l <- layout.fruchterman.reingold(gD)
# l <- layout.auto(gD)
# l <- layout.random(gD)
# l <- layout.circle(gD)
# l <- layout.sphere(gD)
# l <- layout.kamada.kawai(gD)
# l <- layout.reingold.tilford(gD)
# l <- layout.lgl(gD)
# l <- layout.graphopt(gD)

pdf(file=paste("igraph.pdf",sep=""),width = 10,height = 10)
plot(gD,
     vertex.label = V(gD)$name_long,
     vertex.size=log10(V(gD)$size)*3,
     edge.arrow.size=.4,
     asp=TRUE,
     rescale=TRUE,
     layout=l,
     edge.arrow.mode = E(gD)$eorient,
     vertex.label.cex = 0.7,
     vertex.label.dist=0)
#plot(gD, vertex.label = NA, vertex.size=log10(V(gD)$size), edge.arrow.size=.4, asp=TRUE, rescale=TRUE, layout=l)
dev.off()
save(gD, edges, edges.raw, l, file=paste("graph_data.rda",sep=""))
#save(gD, edges, edges.raw, l, file=paste(fname,"graph_data_",upper.fix,"_fix.rda",sep=""))
