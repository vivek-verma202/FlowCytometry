##%######################################################%##
#                                                          #
####                00. load dependencies               ####
#                                                          #
##%######################################################%##
# getwd()
# "/home/vverma3/repo/FlowCytometry"

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
FileNames=list.files(path='./data/fcs/cyto',pattern='.fcs',full.names=TRUE)
# marker names:
MN = pData(parameters(read.FCS(FileNames[1])))[,2]
#Indices corresponding to the columns of FCS files for clustering
ToUse=c(1:10)

# meta dataframe:
df   <- read.csv("./data/FM_pheno.csv")
ID   <- gsub("\\.fcs$","",gsub("./data/fcs/cyto/","",FileNames))
meta <- data.frame(cbind(ID,FileNames))
meta <- merge(meta,df,by="ID",all.x = T)
meta$FM <- factor(meta$FM)
meta$FM <- relevel(meta$FM, ref = "Control")
saveRDS(meta,"./data/cyto_meta.RDS")

##%######################################################%##
#                                                          #
####                02: run clustering                  ####
#                                                          #
##%######################################################%##
rm(list = (c("df","load.lib","ID")));gc()
Build <- RunCl_Flow(S=100,K=20,FileNames=FileNames,
                    doCPF='auto',MN= MN,transformInds=ToUse,
                    ToUse=ToUse,numCore=90)

saveRDS(Build,"./data/cyto_build.RDS")
#extract features
FuncDF=Build[[1]]
FreqDF=Build[[2]]
#functional features


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
saveRDS(umap,"./data/fcs/cyto/cyto_umap.RDS")

# step 3: per-cell differentiation score & visualizations
vizAtlas(CellMat = CM,Build = Build,Y = meta$FM,
         ToUse = ToUse,SampsToUse=NULL,numCore = 90,
         layout = umap,outDir = "./plots/cyto")

# step 4: color by differential frequency
vizAtlas_Freq_Directional(CellMat = CM,Build = Build,
                          FreqMat = FreqDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap,
                          outDir = "./plots/cyto",numIter = 100)

# step 5: Make a functional pval maps
vizAtlas_Func(CellMat = CM,Build = Build,
              FuncMat = FuncDF,Y = meta$FM,
              ToUse = ToUse,
              SampsToUse=NULL,numCore = 90,
              layout = umap,outDir = "./plots/cyto",
              FuncNames = c("CD45RA","Ki.67","Foxp3","IL.4",
                            "IL.2","IL.17A","IFN.g"),
              numIter = 100)

# step 6: Make a functional map, color by direction
vizAtlas_Func_Directional(CellMat = CM,Build = Build,
                          FuncMat = FuncDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap,outDir = "./plots/cyto",
                          FuncNames = c("CD45RA","Ki.67","Foxp3","IL.4",
                                        "IL.2","IL.17A","IFN.g"),
                          numIter = 100)




