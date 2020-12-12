SampCells=function(fcs,ToUse,n,N,cf,numCore){
#Inputs:
	#fcs: The full path name for the file names
	#ToUse: The indices of functional markers to use for clustering
	#n: The number of cells to take from each sample (e.g. 5000)
	#N: The total number of cells (e.g. 30,000)
  # should be <= n x no. of fcs files
  #numCore: number of cores for foreach
#Outputs
	#$PatMark: Functional marker exp in each iteration cluster per patient
	#$PatProp: Frequency data--> prop of cells in each cluster for each sample

  MN=pData(parameters(read.FCS(fcs[1])))[,2]	
  cl=makeCluster(numCore)
  registerDoParallel(cl)

#creating a custom combine for S
custom_s=function(L1,L2){
PatMark=cbind(L1$PatMark,L2$PatMark)
PatProp=cbind(L1$PatProp,L2$PatProp)
ConsFeatByMark=rbind(L1$ConsFeatByMark,L2$ConsFeatByMark)
return(list(PatMark=PatMark,PatProp=PatProp,ConsFeatByMark=ConsFeatByMark))
#return(list(PatMark=PatMark,PatProp=PatProp))
}
custom_i=function(L1,L2){
	SubDM=rbind(L1$SubDM,L2$SubDM)
	AllMarker=rbind(L1$AllMarker,L2$AllMarker)
	return(list(SubDM=SubDM,AllMarker=AllMarker))
}

##############################################
#Step 1: Define data matrix used in clustering
##############################################
DM_Build=foreach(i=fcs,.combine=custom_i,.packages=c('flowCore','plyr')) %dopar% {
  #arcsinh transformation
  cytoftrans=arcsinhTransform(transformationId='cytofTransform',a=0,b=(1/cf),c=0)
  frame=read.FCS(i)
  translist=transformList(colnames(exprs(frame)),cytoftrans)
  frame=transform(frame,translist)
  X=exprs(frame)
  colnames(X)=MN
  #sample n cells from the file
  NCU=min(n,nrow(X))
  SampInds=sample(1:nrow(X),NCU,replace=FALSE)
  AllMarker=X[SampInds,1:length(MN)]
  SubDM=X[SampInds,ToUse]
  colnames(SubDM)=MN[ToUse]
  colnames(AllMarker)=MN
  return(list(SubDM=SubDM,AllMarker=AllMarker))
}
AllMarker=DM_Build$AllMarker
#return only the sub data matrix
ClDM=DM_Build$SubDM
#now subsample N number of cells
SampInds=sample(1:nrow(ClDM),N,replace=FALSE)
ClDM=ClDM[SampInds,]
stopCluster(cl)
ClDM
}
