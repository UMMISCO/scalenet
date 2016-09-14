# Matrix to edge
adj2pairs <- function(mat, symetric=TRUE, silent=FALSE){
  if(is.vector(mat)){
    tmp <- as.data.frame(matrix(NA,nrow = length(mat),ncol = length(mat)))
    rownames(tmp) <- mat; colnames(tmp) <- mat
    mat <- tmp
  }
  square = TRUE
  if(ncol(mat)!=nrow(mat)) {warning("This is not a square matrix!"); square = FALSE}
  edges <- c() # initialize empty
  if(square & symetric){ # no need for diagonal and only the upper triangle is computed
    print("Running edge translation for a symetric square matrix")
    for(row in 1:nrow(mat)){ # for each row
      if(!silent & round(row/nrow(mat)*100) %% 5==0) print(round(row/nrow(mat)*100))
      for(col in 1:ncol(mat)){ # for each column
        if(col>row){
          #print(c(rownames(mat)[row], colnames(mat)[col]))
          edges <- rbind(edges,c(rownames(mat)[row], colnames(mat)[col]))
        }
      }
    }
  }else{ # non symetric. Even the diagonal is computed
    print("Running edge translation for a non symetric and/or non square matrix")
    for(row in 1:nrow(mat)){
      if(!silent & round(row/nrow(mat)*100) %% 5==0) print(round(row/nrow(mat)*100))
      for(col in 1:ncol(mat)){
        #if(col>row){
        edges <- rbind(edges,c(rownames(mat)[row], colnames(mat)[col]))
        #}
      }
    }
  }
  colnames(edges) <- c("node1","node2"); rownames(edges) <- paste("edge",(1:nrow(edges)),sep="")
  return(as.data.frame(edges))
}


# edge to matrix
pairs2adj <- function(edges, symetric=TRUE, silent=FALSE){
  all.nodes <- unique(c(as.character(edges$node1),as.character(edges$node2)))
  rows <- unique(as.character(edges$node1))
  cols <- unique(as.character(edges$node2))
  square <- TRUE # by default
  if(length(cols)!=length(rows)) {
    square <- FALSE
    stop("The matrix will not be square. To be impemented for non square matrixes.")
  }else {
    mat <- as.data.frame(matrix(NA, nrow = length(all.nodes), ncol = length(all.nodes)))
    rownames(mat) <- all.nodes; colnames(mat) <- all.nodes
    # fill the matrix
    for(i in 1:nrow(edges)){
      x <- as.character(edges$node1)[i]
      y <- as.character(edges$node2)[i]
      if(symetric) {
        mat[y,x] <- mat[x,y] <- edges$dist[i]
      }else{
        mat[x,y] <- edges$dist[i]
      }
    }
    warning("The diagonal is empty. Please treat after depending on the method.")
  }
  return(mat)
}

edgeDistance <- function(x,y, method="spearman"){
  if (method %in% c("spearman", "pearson", "jsd") == FALSE)    # error corrected elc
    stop("Please provide a valid mapping type spearman, pearson, jsd")
  if(method=="spearman"){
    res <- cor(rank(x),rank(y), method="pearson")
  }
  if(method=="pearson"){
    res <- cor(x,y, method="pearson")
  }
  if(method=="jsd"){
    pseudocount=1e-11
    xp <- x; xp[x==0] <- pseudocount; yp <- y; yp[y==0] <- pseudocount
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    res <- JSD(xp,yp)
  }
  return(res)
}
