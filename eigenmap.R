
EigenMap <- function(W,reduced.dim=10,bEigs){
# W: the affinity matrix
# ReducedDim: the dimensionality of the reduced subspace
# bEigs: whether to use eigs to speed up


  maxMatrixSize <- 1600
  eigvectorRatio <- 0.1

  W.size = dim(W);
  if(W.size[1] != w.size[2]) stop("W must square matrix!")
  nSmp <- W.size[2]
  reduced.dim <- min(reduced.dim+1,W.size[1])
  D.mhalf <- rowSums(W)^-0.5
  library(Matrix)
  
  D.mhalfMatrix <- Diagonal(length(D.mhalf),D.mhalf)
  W <- D.mhalfMatrix%*%W%*%D.mhalfMatrix
  ind.min <- which(W<t(W))
  W[ind] <- t(W)[ind]

  dim.Matrix <- nrow(W)
  if (missing(bEigs)){
	  if(dim.Matrix>maxMatrixSize&&reduced.dim<dim.Matrix*eigvectorRatio){
		bEigs <- 1
	  } else {
		bEigs <- 0   
	  }
}
  eig.out <- eigen(W)
  if (bEigs) {
    Y <- eig.out$vectors$vectors[,1:reduced.dim]
    eigvalue <- eig.out$values[1:reduced.dim]
  } else {
    Y <- eig.out$vectors
    eigvalue <- eig.out$values
    sort.out <- sort.int(-eigvalue,index.return = TRUE)
    eigvalue <- eigvalue[sort.out$ix]
    Y <- Y[,sort.out$ix]
    if (reduced.dim < length(eigvalue)){
      Y <- Y[,1:reduced.dim]
      eigvalue <- eigvalue[1:reduced.dim]
    }
  }

  eigIdx <- which(abs(eigvalue) > 1e-6)
  
  eigvalue <- eigvalue[eigIdx]
  Y <- Y[,eigIdx]
  nGotDim <- length(eigvalue)

  idx <- 1
  while (abs(eigvalue[idx] - 1) < 1e-12) {
    idx <- idx + 1
    if (idx > nGotDim)
      break
  }
  idx = idx - 1;
  if (idx > 1) {
    u <- matrix(0, nrow(Y), idx)
    d.m <- as.matrix(1/D.mhalf)
    cc <- 1/norm(d.m, type= "F" )
    u[,1] <- cc/D.mhalf
    b.Done <- 0
    for (i in 1:idx) {
      if ( abs(t(Y[,i])%*%u[,1]-1) < 1e-14 ) {
        Y[,i] <- Y[,1]
        Y[,i] <- u[,1]
        b.Done <- 1
      }          
    }
    if (!b.Done) {
      for ( i in 2:idx ) {
        u[,i] <- Y[,i]
        for ( j in 1:i-1 ) {
          u[,i] <- u[,i] - (t(u[,j])%*%Y[,i])%*%u[,j]
        }
        u[,i] <- u[,i]/norm(as.matrix(u[,i]),type = "F")
      }
      Y[,1:idx] <- u
    }
  }
  Y <- D.mhalfMatrix%*%Y
  Y <- Y[,-1]
  eigvalue <- eigvalue[-1]
  return (list(Y=Y, eigvalue = eigvalue))
  }