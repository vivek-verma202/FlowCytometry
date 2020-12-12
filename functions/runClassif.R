runClassif=function(FuncDF,Y,FPV,IterNumClus,propTrain,numPerm,sampID,numCore){
    # inputs:
    # FuncDF: The sample x feature per iteration matrix returned by runRepMetaclust.R
    # Y: the binary response vector
    # FPV: The number of features to use per view.
    # IterNumClus: The number of clusters per iteration
    # propTrain: the proportion of the data to train with
    # numPerm: Number of classification trials to do
    # sampID: Your input will either be 0 or a vector of subject IDs for each sample.
    # Input a vector of subjectIDs for each sample if subject has more than 1 FCS.
    # numCore: The number of cores to use for parallelization

    # Output:
    # StoreAUCs: the vector of AUCs

    # Parameter Specification  ##
    # Parameters for feature selection
    N=10
    #Use the maximal number of clustering iterations provided
    m=length(IterNumClus)

    StoreAUCs=c()
    for(ss in 1:numPerm){
        # Random Forest
        #set up parallelization
        cl=makeCluster(numCore)
        registerDoParallel(cl)

        rfCustom=function(L1,L2){
            ansVec=cbind(L1$ansVec,L2$ansVec)
            modImp=rbind(L1$modImp,L2$modImp)
            return(list(ansVec=ansVec,modImp=modImp))
        }

        kNN=function(Mat,NN){
            k=NN
            dist_mat <- as.matrix(dist(Mat, method = "euclidean", upper = TRUE, diag=TRUE))
            nrst <- lapply(1:nrow(dist_mat), function(i) k.nearest.neighbors(i, dist_mat, k = k))
            w <- matrix(nrow = dim(dist_mat), ncol=dim(dist_mat)) ## all NA right now
            w[is.na(w)] <- 0 ## populate with 0
            for(i in 1:length(nrst)) for(j in nrst[[i]]) w[i,j] = 1
            Adj2=w
            Net=graph.adjacency(Adj2,mode='undirected')
            FinalAdj=get.adjacency(Net,type='both',sparse=FALSE)
            diag(FinalAdj)=0
            FinalAdj
        }

        Laplacian=function(Adj){
            Ds=rowSums(Adj)
            D=matrix(0,nrow=nrow(Adj),ncol=nrow(Adj))
            diag(D)=Ds
            L=D-Adj
            Out=list()
            Out[[1]]=D
            Out[[2]]=L
            Out
        }

        GetTopFeat=function(X,NumKeep,NN){
            #Step 1 is to create kNN graph
            NNGraph=kNN(X,NN)
            #Step 2 is to build similarity matrix
            S=matrix(0,nrow=nrow(NNGraph),ncol=nrow(NNGraph))
            for(ss in 1:nrow(S)){
                Edges=which(NNGraph[ss,]==1)
                for(j in 1:length(Edges)){
                    if(ss<Edges[j]){
                        x1=X[ss,]
                        x2=X[Edges[j],]
                        diff=norm(as.matrix(x1-x2),type='F')
                        diff=(diff^2)/5
                        S[ss,Edges[j]]=exp(-diff)
                        S[Edges[j],ss]=exp(-diff)
                    }
                }
            }
            #Step 3: Compute Laplacian
            GetLaplace=Laplacian(S)
            L=GetLaplace[[2]]
            D=GetLaplace[[1]]
            #Generate matrix of ones
            Ones=as.matrix(rep(1,nrow(S)),ncol=1)
            #Step 4: Compute score for each feature
            FScore=c()
            for(r in 1:ncol(X)){
                #get their feature vector
                fr=as.matrix(X[,r],ncol=1)
                Num=t(fr)%*%D%*%Ones
                Denom=t(Ones)%*%D%*%Ones
                Subtract=c(Num/Denom)*Ones
                TfR=fr-Subtract
                #Compute Laplacian Score
                Num2=t(TfR)%*%L%*%TfR
                Denom2=t(TfR)%*%D%*%TfR
                Lr=Num2/Denom2
                FScore=c(FScore,Lr)
            }
            #get indices for the top scoring features
            TopFeat=order(FScore,decreasing=TRUE)[1:NumKeep]
            #print(TopFeat)
            Out=X[,TopFeat]
            Out
        } ##function end

        MultiviewFS=function(NumIt,DataMat,ClusNum,FPV,NN){
            ConcatPCA=rep(0,nrow(DataMat))
            for(i in 1:NumIt){
                if(i==1){
                    low=1
                }
                else{low=sum(ClusNum[1:(i-1)])+1
                }
                high=sum(ClusNum[1:i])
                IterMat=DataMat[,low:high]
                #prevent from asking for more columns than we have
                FPV2=min(ncol(IterMat),FPV)
                Y1=GetTopFeat(IterMat,FPV2,NN)
                if(ClusNum[i]==0){
                    ConcatPCA=ConcatPCA
                }
                else{
                    ConcatPCA=cbind(ConcatPCA,Y1)}
            }
            ConcatPCA=ConcatPCA[,-1]
            ConcatPCA
        }

        #number of CV iterations to do
        itSeq=1:500
        X2=MultiviewFS(m,FuncDF,IterNumClus,FPV,N)
        RFIter=foreach(s=itSeq,.combine=rfCustom,.packages='randomForest') %dopar% {
            UniqueSamp=unique(sampID)
            repeat{
                trainGet=sample(UniqueSamp,length(UniqueSamp)*propTrain)
                print(trainGet)
                #get all sample IDs that is an element of this
                relSampInds=which(is.element(sampID,trainGet))

                train=relSampInds
                if(length(unique(Y[train]))==2)
                    break
            }
            test=seq(length(Y))[-train]
            ansVec=rep(NA,length(Y))
            mod=randomForest(X2[train,],Y[train])
            predVec=predict(mod,X2[test,],type='prob')[,2]
            ansVec[test]=predVec
            return(list(ansVec=ansVec))
        }
        stopCluster(cl)
        ansMat=RFIter$ansVec
        ans=as.numeric(rowMedians(ansMat,na.rm=TRUE))
        #calculate AUC
        pred=prediction(ans,Y)
        r1=performance(pred,measure='auc')
        r=r1@y.values[[1]]
        print('trial #')
        print(ss)
        print(r)
        StoreAUCs=c(StoreAUCs,r)
    } #for ss
    StoreAUCs
} #function end
