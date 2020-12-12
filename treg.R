##%######################################################%##
#                                                          #
####                00. load dependencies               ####
#                                                          #
##%######################################################%##
# getwd()
setwd("/home/vverma3/repo/FlowCytometry")

# load libraries
load.lib <- c("flowCore",
              "FlowSOM",
              "foreach",
              "doParallel",
              "iterators",
              "Rtsne",
              "plyr",
              "igraph",
              "randomForest",
              "matrixStats",
              "ROCR",
              "FastKNN",
              "miscTools",
              "ggplot2",
              "reshape2",
              "viridis",
              "pROC",
              "scales",
              "uwot",
              "FNN")
sapply(load.lib,require,character=T)
# >>>> load functions
R.utils::sourceDirectory(path = "./functions/")

##%######################################################%##
#                                                          #
####                    01. load data                   ####
#                                                          #
##%######################################################%##

## FCS files
fcs_nk    <- list.files(path='./data/fcs/nk',pattern='.fcs',full.names=T)
fcs_bmd   <- list.files(path='./data/fcs/bmd',pattern='.fcs',full.names=T)
fcs_cyto  <- list.files(path='./data/fcs/cyto',pattern='.fcs',full.names=T)
fcs_chemo <- list.files(path='./data/fcs/chemo',pattern='.fcs',full.names=T)
fcs_treg  <- list.files(path='./data/fcs/treg',pattern='.fcs',full.names=T)
# marker names:
MN_nk    <- pData(parameters(read.FCS(fcs_nk[1])))[,2]
MN_bmd   <- pData(parameters(read.FCS(fcs_bmd[1])))[,2]
MN_cyto  <- pData(parameters(read.FCS(fcs_cyto[1])))[,2]
MN_chemo <- pData(parameters(read.FCS(fcs_chemo[1])))[,2]
MN_treg  <- pData(parameters(read.FCS(fcs_treg[1])))[,2]

#Indices corresponding to the columns of FCS files for clustering
ToUse_nk    <- c(1:length(MN_nk))
ToUse_bmd   <- c(1:length(MN_bmd))
ToUse_cyto  <- c(1:length(MN_cyto))
ToUse_chemo <- c(1:length(MN_chemo))
ToUse_treg  <- c(1:length(MN_treg))

# meta dataframe:
ID     <- gsub("\\.fcs$","",gsub("./data/fcs/nk/","",fcs_nk))
meta1  <- data.frame(cbind(ID = ID,fcs = fcs_nk,
                           panel = rep("nk",length(ID))))
ID     <- gsub("\\.fcs$","",gsub("./data/fcs/bmd/","",fcs_bmd))
meta2  <- data.frame(cbind(ID = ID,fcs = fcs_bmd,
                           panel = rep("bmd",length(ID))))
ID     <- gsub("\\.fcs$","",gsub("./data/fcs/cyto/","",fcs_cyto))
meta3  <- data.frame(cbind(ID = ID,fcs = fcs_cyto,
                           panel = rep("cyto",length(ID))))
ID     <- gsub("\\.fcs$","",gsub("./data/fcs/chemo/","",fcs_chemo))
meta4  <- data.frame(cbind(ID = ID,fcs = fcs_chemo,
                           panel = rep("chemo",length(ID))))
ID     <- gsub("\\.fcs$","",gsub("./data/fcs/treg/","",fcs_treg))
meta5  <- data.frame(cbind(ID = ID,fcs = fcs_treg,
                           panel = rep("treg",length(ID))))
meta   <- data.frame(rbind(meta1,meta2,meta3,meta4,meta5))
df     <- read.csv("./data/FM_pheno.csv")
meta   <- merge(meta,df,by="ID",all.x = T)
meta$FM <- relevel(factor(meta$FM), ref = "Control")
saveRDS(meta,"./data/meta.RDS")
meta <- readRDS("./data/meta.RDS")

##%######################################################%##
#                                                          #
####                02: run clustering                  ####
#                                                          #
##%######################################################%##
rm(list = (c("df","load.lib","ID","meta1","meta2","meta3","meta4","meta5")))
gc()

Build_nk    <- RunCl_Flow(S=100,K=20,FileNames=fcs_nk,
                       doCPF='auto',MN=MN_nk,transformInds=ToUse_nk,
                       ToUse=ToUse_nk,numCore=90)
Build_bmd   <- RunCl_Flow(S=100,K=20,FileNames=fcs_bmd,
                       doCPF='auto',MN=MN_bmd,transformInds=ToUse_bmd,
                       ToUse=ToUse_bmd,numCore=90)
Build_cyto  <- RunCl_Flow(S=100,K=20,FileNames=fcs_cyto,
                        doCPF='auto',MN=MN_cyto,transformInds=ToUse_cyto,
                        ToUse=ToUse_cyto,numCore=90)
Build_chemo <- RunCl_Flow(S=100,K=20,FileNames=fcs_chemo,
                         doCPF='auto',MN=MN_chemo,transformInds=ToUse_chemo,
                         ToUse=ToUse_chemo,numCore=90)
Build_treg  <- RunCl_Flow(S=100,K=20,FileNames=fcs_treg,
                         doCPF='auto',MN=MN_treg,transformInds=ToUse_treg,
                         ToUse=ToUse_treg,numCore=90)
gc()
#extract features
FuncDF=Build[[1]]
FreqDF=Build[[2]]

##%######################################################%##
#                                                          #
####                 03. Visualizations                 ####
#                                                          #
##%######################################################%##

# step 1: Sample cells across all FCS files
CM_nk    <- SampCells(fcs = fcs_nk,ToUse = ToUse_nk,n=5000,N=50000,
                   cf=200,numCore=90); gc()
CM_bmd   <- SampCells(fcs = fcs_bmd,ToUse = ToUse_bmd,n=5000,N=50000,
                   cf=200,numCore=90); gc()
CM_cyto  <- SampCells(fcs = fcs_cyto,ToUse = ToUse_cyto,n=5000,N=50000,
                   cf=200,numCore=90); gc()
CM_chemo <- SampCells(fcs = fcs_chemo,ToUse = ToUse_chemo,n=5000,N=50000,
                   cf=200,numCore=90); gc()
CM_treg  <- SampCells(fcs = fcs_treg,ToUse = ToUse_treg,n=5000,N=50000,
                   cf=200,numCore=90); gc()
# step 2: run UMAP
umap_nk    <- tumap(CM_nk, n_threads = 90,  verbose = T); gc()
umap_bmd   <- tumap(CM_bmd, n_threads = 90,  verbose = T); gc()
umap_cyto  <- tumap(CM_cyto, n_threads = 90,  verbose = T); gc()
umap_chemo <- tumap(CM_chemo, n_threads = 90,  verbose = T); gc()
umap_treg  <- tumap(CM_treg, n_threads = 90,  verbose = T); gc()
  # check if umap looks good
ggplot(data.frame(umap_treg))+aes(x = X1, y = X2, color = runif(50000))+ 
  geom_point(size = 1, alpha = 0.3) + scale_color_viridis() +
  theme_bw() 

# step 3: meatcluster profile & pvalues
vizAtlas(CellMat = CM_nk,Build = Build_nk,Y = meta$FM[meta$panel=="nk"],
         ToUse = ToUse_nk,SampsToUse=NULL,numCore = 90,
         layout = umap_nk,outDir = "./plots/nk"); gc()
vizAtlas(CellMat = CM_bmd,Build = Build_bmd,Y = meta$FM[meta$panel=="bmd"],
         ToUse = ToUse_bmd,SampsToUse=NULL,numCore = 90,
         layout = umap_bmd,outDir = "./plots/bmd"); gc()
vizAtlas(CellMat = CM_cyto,Build = Build_cyto,Y = meta$FM[meta$panel=="cyto"],
         ToUse = ToUse_cyto,SampsToUse=NULL,numCore = 90,
         layout = umap_cyto,outDir = "./plots/cyto"); gc()
vizAtlas(CellMat = CM_chemo,Build = Build_chemo,Y = meta$FM[meta$panel=="chemo"],
         ToUse = ToUse_chemo,SampsToUse=NULL,numCore = 90,
         layout = umap_chemo,outDir = "./plots/chemo"); gc()
vizAtlas(CellMat = CM_treg,Build = Build_treg,Y = meta$FM[meta$panel=="treg"],
         ToUse = ToUse_treg,SampsToUse=NULL,numCore = 90,
         layout = umap_treg,outDir = "./plots/treg"); gc()
# step 4: color by differential frequency
vizAtlas_Freq_Directional(CellMat = CM,Build = Build,
                          FreqMat = FreqDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap,
                          outDir = "./plots/treg",numIter = 100)

# step 5: Make a functional pval maps
vizAtlas_Func(CellMat = CM,Build = Build,
              FuncMat = FuncDF,Y = meta$FM,
              ToUse = ToUse,
              SampsToUse=NULL,numCore = 90,
              layout = umap,outDir = "./plots/treg",
              FuncNames = c("CD45RA","Ki.67","CXCR5","CXCR3",
                            "CCR6","CRTH2"),
              numIter = 100)

# step 6: Make a functional map, color by direction
vizAtlas_Func_Directional(CellMat = CM,Build = Build,
                          FuncMat = FuncDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap,outDir = "./plots/treg",
                          FuncNames = c("CD45RA","Ki.67","CXCR5","CXCR3",
                                        "CCR6","CRTH2"),
                          numIter = 100)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#purpose: generate a distribution of AUC
IterNumClus=Build$IterNumClus
#rename to not perturb RF
ForC=1:ncol(FuncDF)
NewName=paste(colnames(FuncDF),ForC,sep='_')
colnames(FuncDF)=NewName
Y=meta$FM

#do baseline distribution first
Base=c()
for(i in 1:100){
  #choose an iteration to get features from and get relavent columns
  rndIter=sample(1:length(IterNumClus),1)
  NumClus=IterNumClus[1]
  start=((rndIter-1)*NumClus)+1
  end=NumClus*rndIter
  SubDF=FuncDF[,start:end]
  AUC=runClassif(SubDF,meta$FM,IterNumClus[1],
                 IterNumClus[1],0.5,1,as.character(meta$ID),35)
  Base=c(Base,AUC)
}

#get repeated distribution
Boot=runClassif(FuncDF,meta$FM,40,IterNumClus,0.5,100,as.character(meta$ID),90)


#plotting
AllBL=rbind(Base,Boot)
rownames(AllBL)=c('Baseline','BootStrap')
DF=melt(AllBL)
names(DF)=c('Method','perm','AUC')

ValVec=c('black','#FF0266')
p26 <- ggplot(DF, aes(x=Method, y=AUC,color=Method)) +
  geom_boxplot(notch=FALSE,lwd=1)+theme(text = element_text(size=15))+xlab('')+ylab('')+scale_color_manual(values = ValVec)+ggtitle('')
p26=p26+geom_point(size=.7,position=position_jitterdodge(),alpha=.5)+theme_classic()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size=14))      
p26=p26+theme(legend.position='bottom')+ggtitle('')+ylab('AUC')
ggsave('OutDir/Preg_Dist.pdf',p26,width=4,height=4)




