setwd("file")
library(Seurat)
library(tidyverse)
library(harmony)
library(ggplot2)
library(monocle)
library(patchwork)
library(readxl)
library(openxlsx)
library(scuttle)
library(RColorBrewer)
library(ggsci)
library(clusterProfiler)
library(enrichR)
library(enrichR)
library(farver)
library(RColorBrewer)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(devtools)
library(scRNAtoolVis)
library(EnhancedVolcano)

#hy ：
HYCT<-readRDS("SMT_new_gene.rds")
DTCT<-readRDS("DTCT.rds")
dt<-DTCT
hy<-HYCT
dt_to_hy <- read.table('dt_to_hy_2.txt') 
dt_to_hy <- dt_to_hy[,c(1,2)]
names(dt_to_hy) <- c('dt','hy')
hy_to_dt <- read.table('hy_to_dt_2.txt') 
hy_to_dt <- hy_to_dt[,c(2,1)]
names(hy_to_dt) <- c('hy','dt') 
hy_to_dt_all <- rbind(dt_to_hy,hy_to_dt)
hy_to_dt_all <- hy_to_dt_all[!duplicated(hy_to_dt_all$dt),]
matrix_dt <- dt@assays$RNA@counts

hy_to_dt_all$dt <- stringr::str_replace_all(hy_to_dt_all$dt,c('evm.model.'=''))
rownames(matrix_dt) <- stringr::str_replace_all(rownames(matrix_dt),c('EVM%20prediction%20'='','-'='_'))

index <- intersect(hy_to_dt_all$dt,rownames(matrix_dt))
rownames(hy_to_dt_all) <- hy_to_dt_all$dt
hy_to_dt_all <- hy_to_dt_all[index,]
matrix_dt_new <- matrix_dt[index,]
rownames(matrix_dt_new) <- hy_to_dt_all$hy
matrix_dt_new <- matrix_dt_new[!duplicated(rownames(matrix_dt_new)), ]

dim(matrix_dt_new)
dim(hy)


matrix_hy<-hy@assays$RNA@counts
data<-list(matrix_dt_new,matrix_hy)
name<-c("DTCT","HYCT")
names(data)<-name
scRNA.list<-list()
for (i in name) {
  scRNA.list[[i]] <- data[[i]]
  scRNA.list[[i]]<-CreateSeuratObject(scRNA.list[[i]],project = i)}
scRNA<-merge(scRNA.list[[1]],scRNA.list[2:length(scRNA.list)])
#table(scRNA$orig.ident)
#scRNA<-SCTransform(scRNA)
#scRNA<-RunPCA(scRNA,assay = "SCT")
#scRNA<-RunUMAP(scRNA,dims = 1:30)
#scRNA<-FindNeighbors(scRNA,dims = 1:30)%>%FindClusters()
#DimPlot_scCustom(scRNA,label = F,group.by = "orig.ident")+ggtitle("")
#saveRDS(scRNA,file = "haircellDTvsHY_nonCCA.rds")

scRNAlist<-SplitObject(scRNA,split.by = "orig.ident")
scRNAlist<-parallel::mclapply(scRNAlist,FUN = function(x) SCTransform(x),mc.cores = 8)
scRNA.features<-SelectIntegrationFeatures(scRNAlist,nfeatures = 2000)
scRNAlist<-PrepSCTIntegration(scRNAlist,anchor.features = scRNA.features)
scRNA.anchors<-FindIntegrationAnchors(object.list = scRNAlist,normalization.method = "SCT",anchor.features = scRNA.features)
scRNA.sct.int<-IntegrateData(scRNA.anchors,normalization.method = "SCT")
scRNA<-RunPCA(scRNA.sct.int,npcs = 30,verbose = T)
scRNA<-scRNA%>%RunUMAP(dims = 1:30)
scRNA<-FindNeighbors(scRNA)
scRNA<-FindClusters(scRNA,resolution = c(0.1,0.2,0.3,0.4,0.5,0.6))
library(scCustomize)
DimPlot_scCustom(scRNA,label = T,group.by = "integrated_snn_res.0.3",pt.size = 1,label.size = 6)
Idents(scRNA)<-"integrated_snn_res.0.3"
p2 = DimPlot(scRNA,label = T,split.by = "orig.ident",pt.size = 1,label.size = 6)
p1/p2
saveRDS(scRNA,file = "crossspeciesDT_vs_HY.rds")