Clust=function(i,k,cf,fcs,numCore){
#Description: Metaclustering code to partition cells into k metaclusters
#Inputs:
	# i: The number of iterations to do (e.g. 50)
	# k: The number of metaclusters (e.g. 50)
  # cf: co-factor for asinh transformation
	# fcs: The full path for the fcs files (e.g. c('File1.fcs','File2.fcs'))
	# numCore: number of cores for foreach
  MN=pData(parameters(read.FCS(fcs[1])))[,2]
  cl=makeCluster(numCore)
  registerDoParallel(cl)
  #custom combine for within-file clustering
  custom=function(LL1,LL2){
    Assn=c(LL1$Assn,LL2$Assn)
    Centers=rbind(LL1$Centers,LL2$Centers)
    UClus=c(LL1$UClus,LL2$UClus)
    NumCells=c(LL1$NumCells,LL2$NumCells)
    XSample=rbind(LL1$XSample,LL2$XSample)
    return(list(XSample=XSample,Assn=Assn,Centers=Centers,UClus=UClus,NumCells=NumCells))
  }
  DM_Cluster=foreach(i=fcs,.combine=custom,.packages=c('flowCore','plyr','FlowSOM')) %dopar% {
    #arcsinh transformation
    cytoftrans=arcsinhTransform(transformationId='cytofTransform',a=0,b=(1/cf),c=0)
    frame=read.FCS(i)
    translist=transformList(colnames(exprs(frame)),cytoftrans)
    frame=transform(frame,translist)
    X=scale(exprs(frame))
    colnames(X)=MN
    # within file clustering
    NumCl=floor(sqrt((nrow(X)/2)))
    somdim=floor(sqrt(NumCl))
    som=SOM(X, xdim=somdim, ydim=somdim)
    # collect results
    Assn=c(som[["mapping"]][,1])
    Assn2=as.factor(Assn)
    Centers=som$codes
    UClus=nrow(Centers)
    NumCells=nrow(X)
    X=data.frame(X,Assn2)
    XSample=ddply(X,.(Assn2),numcolwise(median))
    XSample=as.matrix(XSample[,-1])
    return(list(XSample=XSample,Assn=Assn,Centers=Centers,UClus=UClus,NumCells=NumCells))
  }
  stopCluster(cl)
  
  # begin metaclustering
  cl=makeCluster(numCore)
  registerDoParallel(cl)
  custom_s=function(L1,L2){
    FMed=cbind(L1$FMed,L2$FMed)
    FProp=cbind(L1$FProp,L2$FProp)
    IterNumClus=c(L1$IterNumClus,L2$IterNumClus)
    ClusMed=rbind(L1$ClusMed,L2$ClusMed)
    return(list(FMed=FMed,FProp=FProp,IterNumClus=IterNumClus,ClusMed=ClusMed))
  }
  SecondDM=DM_Cluster$Centers
  SecondDM2=scale(SecondDM)
  XSampFull=DM_Cluster$XSample
  print('SOM clustering done, metaclustering now')
  
  BootMeta=foreach(s=1:i,.combine='custom_s',.packages=c('foreach','flowCore','plyr')) %dopar%{
    #actual metaclustering
    MetaClust=kmeans(SecondDM2,centers=k)
    MetaAssn=MetaClust$cluster
    MetaCenter=MetaClust$centers
    #calculate the median marker expression for each cluster
    toAssn=as.factor(MetaAssn)
    XX=data.frame(XSampFull,toAssn)
    names(XX)[ncol(XX)]='Assn'
    ClusMed=ddply(XX,.(toAssn),numcolwise(median))
    ClusMed=as.matrix(ClusMed[,-1])
    IterNumClus=k
    #mapping metacluster labels back to single cells
    #Individual file cluster centers were assigned to metaclusters
    #Individual file clusters are comprised of cells in each file
    UClus=DM_Cluster$UClus
    Assn=DM_Cluster$Assn
    NumCells=DM_Cluster$NumCells
    #Assn keeps track of within-file cluster assignment for each cell across files
    #We will update sample-specific entries of Assn to metacluster labels
    #Index book keeping
    old=0
    cellOld=0
    startEndMat=matrix(0,nrow=length(UClus),ncol=2)
    for(i in 1:length(UClus)){
      Start=old+1
      End=Start+(UClus[i]-1)
      old=End
      startCell=cellOld+1
      endCell=startCell+(NumCells[i]-1)
      cellOld=endCell
      startEndMat[i,1]=startCell
      startEndMat[i,2]=endCell
      RelAssn=Assn[startCell:endCell]
      clIndexes=Start:End
      Convert=clIndexes[RelAssn]
      Assn[startCell:endCell]=Convert
    }
    #CellMetaLab is the converted version ot cell-to-metacluster labels
    CellMetaLab=rep(0,length(Assn))
    c_i=function(L1,L2){
      calcMed=rbind(L1$calcMed,L2$calcMed)
      propVec=rbind(L1$propVec,L2$propVec)
      CellMetaLab=cbind(L1$CellMetaLab,L2$CellMetaLab)
      return(list(calcMed=calcMed,propVec=propVec,CellMetaLab=CellMetaLab))
    }
    c_j=function(L1,L2){
      calcMed=c(L1$calcMed,L2$calcMed)
      propVec=c(L1$propVec,L2$propVec)
      return(list(calcMed=calcMed,propVec=propVec))
    }
    Base2PFeat=foreach(i=1:length(fcs),.combine=c_i,.packages='foreach') %do% {
      startCell=startEndMat[i,1]
      endCell=startEndMat[i,2]
      GetAssn=Assn[startCell:endCell]
      WorkWith=MetaAssn[GetAssn]
      CellMetaLab[startCell:endCell]=WorkWith
      ClsMeta=c(1:max(MetaAssn))
      #begin calculating frequency and functional marker features
      foreach(j=ClsMeta,.combine=c_j) %do% {
        relInds=which(WorkWith==j)
        stage2Assn=GetAssn[relInds]
        tabStage2Assn=table(stage2Assn)
        #keeps track of how many cells in each file were assigned to each individual cluster
        Weights=t(as.matrix(tabStage2Assn))
        IndexWeight=as.numeric(names(tabStage2Assn))
        if(length(unique(stage2Assn))==0){
          calcMed=rep(0,ncol(XSampFull))
        }
        else if(length(unique(stage2Assn))==1){
          calcMed=as.matrix(XSampFull[IndexWeight,])
          calcMed=calcMed[,1]
        }else{
          #get weighted mean
          WeightMed=Weights%*%as.matrix(XSampFull[IndexWeight,])
          calcMed=WeightMed/sum(Weights)
          calcMed=calcMed[1,]
        }
        #get names for medians
        Name4Med=paste(j,colnames(XSampFull),sep='_')
        names(calcMed)=Name4Med
        #calculate frequencies
        propVec=length(relInds)/NumCells[i]
        names(propVec)=paste(j,'--prop',sep='')
        return(list(calcMed=calcMed,propVec=propVec))
      } #j
    } #i
    FMed=Base2PFeat$calcMed
    FProp=Base2PFeat$propVec
    return(list(FMed=FMed,FProp=FProp,IterNumClus=IterNumClus,ClusMed=ClusMed))
  } #s
  stopCluster(cl)
  BootMeta
}

