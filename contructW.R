constructW <- function(fea,options.k=5){

# Usage: construct the similarity matrix using KNN neighboring:
#                'KNN'            -  k = 0
#                                       Complete graph
#                                    k > 0
#                                      Put an edge between two nodes if and
#                                      only if they are among the k nearst
#                                      neighbors of each other. You are
#                                      required to provide the parameter k in
#                                      the options. Default k=5.
#    For more details about the different ways to construct the W, please
#    refer:
#    Zixing Wang, Peng Qiu, Wenlong Xu, Yin Liu, "Spectral feature selection and its application in high dimensional gene expression studies" ACM BCB Sep. 2014.

if (missing(options.k))
   stop("options.k must be a number ")		
nSmp = size(fea,1);
maxM = 62500000; #500M
bSpeed  = 1;
BlockSize = floor(maxM/(nSmp*3));	
G = matrix(0,nSmp*(options.k+1),3)
for (i in 1:ceiling(nSmp/BlockSize)){
	if (i == ceiling(nSmp/BlockSize)){
		smpIdx = (i-1)*BlockSize+1:nSmp
		dist = fea[smpIdx,:]*t(fea)
		if bSpeed {
			nSmpNow = length(smpIdx)
			dump = matrix(0,nSmpNow,options.k+1)
			idx = dump;
			for (j in 1:options.k+1){
				idx[:,j] = max.col(dist)
				temp = (idx[:,j]-1)%*%nSmpNow + matrix(1:nSmpNow,nSmpNow,1)
				dist[temp] = 0;
			}
		} else {
			idx = t(apply(-dist,1,order))
			idx = idx[:,1:options.k+1];
			dump = -dump[:,1:options.k+1];				
		}
		
		G[(i-1)*BlockSize*(options.k+1)+1:(nSmp*(options.k+1)),1] = apply(t(smpIdx),2,rep,options.k+1) 
		G[(i-1)*BlockSize*(options.k+1)+1:(nSmp*(options.k+1)),2] = c(idx);
		G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = c(dump);
		
	}else{
		smpIdx = (i-1)*BlockSize+1:i*BlockSize;
		dist = fea[smpIdx,:] %*% t(fea)                
		if bSpeed{
			nSmpNow = length(smpIdx);
			dump = matrix(0,nSmpNow,options.k+1);
			idx = dump;
			for (j in 1:options.k+1){
				idx[:,j] = max.col(dist)
				temp = (idx[:,j]-1)%*%nSmpNow + matrix(1:nSmpNow,nSmpNow,1)
				dist[temp] = 0;
				}
			} else {
				idx = t(apply(-dist,1,order))
				idx = idx[:,1:options.k+1];
				dump = -dump[:,1:options.k+1];						
			}
		G[(i-1)*BlockSize*(options.k+1)+1:(nSmp*(options.k+1)),1] = apply(t(smpIdx),2,rep,options.k+1) 
		G[(i-1)*BlockSize*(options.k+1)+1:(nSmp*(options.k+1)),2] = c(idx);
		G((i-1)*BlockSize*(options.k+1)+1:nSmp*(options.k+1),3) = c(dump)			
	}		
} 
W <- matrix(0,nSmp,nSmp,sparse=TRUE)
W[cbind(G[:,1],G[,2])] <- G[:,3]	
return W	
}

  