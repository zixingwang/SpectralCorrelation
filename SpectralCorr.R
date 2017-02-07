 SpectralCorr <- function(sample,k){
# %                   
# %         sample        -  the sample matrix with n samples, p genes. nXp
# %         options: Struct value in Matlab. The fields in options that can be set:
# %                  
# %               
# %            k         -   The parameter needed under 'KNN' NeighborMode.


W = constructW(sample,k);

#%%%%%%%%%%%%%%%%%%%%%%%%%%%DECLUSTERING%%%%%%%%%%%%%
distcorr <- function(X,y) {
  n <- nrow(X)
  distX <- as.matrix(dist(X))
  disty <- as.matrix(dist(y))
  colMX <- matrix(colMeans(distX),n,n)
  colMy <- matrix(colMeans(disty),n,n)
  distX <- distX - colMX - t(colMX) + matrix(mean(colMX),n,n)
  disty <- disty - colMy - t(colMy) + matrix(mean(colMy),n,n)
  dCov <- sum(distX*disty)/n^2
  dVarX <- sum(distX*distX)/n^2
  dVary <- sum(disty*disty)/n^2
  return(sqrt(dCov/sqrt(dVarX*dVary)))
}

##########
Y = Eigenmap(W,k);
Dist=matrix(0,ncol(sample));
for (indx in 1:ncol(sample)){
    Dist[indx]=distcorr(Y,sample[:,indx]);    
}

Feature_order= Dist[order(-Dist)];
return Feature_order
}