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
              "reticulate",
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
FileNames=list.files(path='./data/fcs/bmd',pattern='.fcs',full.names=TRUE)
# marker names:
MN = pData(parameters(read.FCS(FileNames[1])))[,2]
#Indices corresponding to the columns of FCS files for clustering
ToUse=c(1:14)

# meta dataframe:
df   <- read.csv("./data/FM_pheno.csv")
ID   <- gsub("\\.fcs$","",gsub("./data/fcs/bmd/","",FileNames))
meta <- data.frame(cbind(ID,FileNames))
meta <- merge(meta,df,by="ID",all.x = T)
meta$FM <- factor(meta$FM)
meta$FM <- relevel(meta$FM, ref = "Control")
saveRDS(meta,"./data/bmd_meta.RDS")

##%######################################################%##
#                                                          #
####                02: run clustering                  ####
#                                                          #
##%######################################################%##
rm(list = (c("df","load.lib","ID")));gc()
Build <- RunCl_Flow(S=100,K=50,FileNames=FileNames,
                    doCPF='auto',MN= MN,transformInds=ToUse,
                    ToUse=ToUse,numCore=90)

saveRDS(Build,"./data/bmd_build.RDS")
#frequency features
FreqDF=Build[[2]]
#functional features
FuncDF=Build[[1]]

##%######################################################%##
#                                                          #
####                 03. Visualizations                 ####
#                                                          #
##%######################################################%##

# step 1: Sample cells across all FCS files
CM <- SampCells(fcs = FileNames,ToUse = ToUse,n=5000,N=50000,
                cf=200,numCore=90)
# step 2: run UMAP
umap <- umap::umap(CM, n_threads = 0)
# check if umap looks good
df <- data.frame(umap$layout)
# testing ggplot2  
ggplot(df)+aes(x = X1, y = X2)+ 
  geom_point(size = 1, alpha = 0.3, color = "magenta") + theme_bw() + theme(legend.position = "none")
saveRDS(umap,"./data/fcs/bmd/bmd_umap.RDS")

# step 3: per-cell differentiation score & visualizations
vizAtlas(CellMat = CM,Build = Build,Y = meta$FM,
         ToUse = ToUse,SampsToUse=NULL,numCore = 90,
         layout = umap$layout,outDir = "./plots/bmd")

# step 4: color by differential frequency
vizAtlas_Freq_Directional(CellMat = CM,Build = Build,
                          FreqMat = FreqDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap$layout,
                          outDir = "./plots/bmd",numIter = 100)

# step 5: Make a functional pval maps
vizAtlas_Func(CellMat = CM,Build = Build,
              FuncMat = FuncDF,Y = meta$FM,
              ToUse = ToUse,
              SampsToUse=NULL,numCore = 90,
              layout = umap$layout,outDir = "./plots/bmd",
              FuncNames = c("CD123","CD19","CD14","CD11c",
                            "CD21","CD38","CD56","CD16","CD10",
                            "CD27","CD282","IgD","CD141","HLADR"),
              numIter = 100)

# step 6: Make a functional map, color by direction
vizAtlas_Func_Directional(CellMat = CM,Build = Build,
                          FuncMat = FuncDF,Y = meta$FM,
                          ToUse = ToUse,
                          SampsToUse=NULL,numCore = 90,
                          layout = umap$layout,outDir = "./plots/bmd",
                          FuncNames = c("CD123","CD19","CD14","CD11c",
                                        "CD21","CD38","CD56","CD16","CD10",
                                        "CD27","CD282","IgD","CD141","HLADR"),
                          numIter = 100)















