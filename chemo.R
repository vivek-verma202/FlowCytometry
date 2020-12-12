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
FileNames=list.files(path='./data/fcs/chemo',pattern='.fcs',full.names=TRUE)
# marker names:
MN = pData(parameters(read.FCS(FileNames[1])))[,2]
#Indices corresponding to the columns of FCS files for clustering
ToUse=c(1:10)

# meta dataframe:
df   <- read.csv("./data/FM_pheno.csv")
ID   <- gsub("\\.fcs$","",gsub("./data/fcs/chemo/","",FileNames))
meta <- data.frame(cbind(ID,FileNames))
meta <- merge(meta,df,by="ID",all.x = T)
meta$FM <- factor(meta$FM)
meta$FM <- relevel(meta$FM, ref = "Control")
saveRDS(meta,"./data/chemo_meta.RDS")

##%######################################################%##
#                                                          #
####                02: run clustering                  ####
#                                                          #
##%######################################################%##
rm(list = (c("df","load.lib","ID")));gc()
Build <- RunCl_Flow(S=100,K=20,FileNames=FileNames,
                    doCPF='auto',MN= MN,transformInds=ToUse,
                    ToUse=ToUse,numCore=90)
saveRDS(Build,"./data/chemo_build.RDS")
#extract features
FuncDF=Build[[1]]
FreqDF=Build[[2]]

##%######################################################%##
#                                                          #
####                 03. Visualizations                 ####
#                                                          #
##%######################################################%##

# step 1: Sample cells across all FCS files
CM <- SampCells(fcs = FileNames,ToUse = ToUse,n=5000,N=50000,
                cf=200,numCore=90)
# step 2: run UMAP
umap <- tumap(CM, n_threads = 90,  verbose = T)
# check if umap looks good
df <- data.frame(umap)
# testing ggplot2  
ggplot(df)+aes(x = X1, y = X2)+ 
  geom_point(size = 1, alpha = 0.3, color = "steelblue") +
  theme_bw() + theme(legend.position = "none")
saveRDS(umap,"./data/fcs/chemo/chemo_umap.RDS")

# step 3: per-cell differentiation score & visualizations
vizAtlas(CellMat = CM,Build = Build,Y = meta$FM,
         ToUse = ToUse,SampsToUse=NULL,numCore = 90,
         layout = umap,outDir = "./plots/chemo")

# step 4: color by differential frequency
vizAtlas_Freq_Directional(CellMat = CM,Build = Build,
                          FreqMat = FreqDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap,
                          outDir = "./plots/chemo",numIter = 100)

# step 5: Make a functional pval maps
vizAtlas_Func(CellMat = CM,Build = Build,
              FuncMat = FuncDF,Y = meta$FM,
              ToUse = ToUse,
              SampsToUse=NULL,numCore = 90,
              layout = umap,outDir = "./plots/chemo",
              FuncNames = c("CD45RA","Ki.67","CXCR5","CXCR3",
                            "CCR6","CRTH2"),
              numIter = 100)

# step 6: Make a functional map, color by direction
vizAtlas_Func_Directional(CellMat = CM,Build = Build,
                          FuncMat = FuncDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap,outDir = "./plots/chemo",
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




