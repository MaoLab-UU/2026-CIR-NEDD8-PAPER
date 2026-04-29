library(tidyverse)
library(readxl)
library(Seurat)
library(UCell)
library(glmGamPoi)
library(patchwork)
library(CellChat)
library(splitstackshape)
library(patchwork)
library(hexbin)
library(jsonlite)
library(stringr)
library(reshape2)
library(infercnv)
library(reticulate)
library(gridExtra)
library(DESeq2)
library(ComplexHeatmap)
library(geneconverter)
library(scales)
library(SingleCellExperiment)
library(slingshot)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(limma)
library(org.Mm.eg.db)
library(gprofiler2)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 32000*1024^2)
lapply(list.files("C:/bioinfo/R tools",full.names = T),source) %>% invisible()

obj <- read_rds('1.general.ident.assigned.rds')


## subset cluster
obj <- obj %>% subset(celltype.general.1 == 'T cell')

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL
## re normalize and cluster

obj <- SCTransform(obj,vars.to.regress = c("mitoRatio"),clip.range = c(-5,5))

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')

### seperate CD4 CD8

FeaturePlot(obj,c('Cd8a','Cd4','Cd3d','Cd3e','Col6a1','Krt8','Itgam','Itgax','Cd14')) | 
DimPlot(obj,label=T)

DotPlot(obj,list(Tcell = c('Cd8a','Cd4','Cd3d','Cd3e'),
            stromal =  c('Col6a1','Krt8'),
            myeloid = c('Itgam','Itgax','Cd14'),
            `NK/T`=  c('Ncr1','Klrb1c'),
             gdT = c('Trdc','Il2ra','Il7r'),
            Naive = c('Cd44','Tcf7','Sell','Lef1')))+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))





c <- FindMarkers(obj,ident.1 = 2)

### assign cell identity

idents <- read_excel('Step2 Lymphoid cell/T cell Prelim not shown/step2 lymphoid cell identity.xlsx')

obj$lymphoid.cell.type <- 'NA'

for(n in 1:nrow(idents)){
  id <- idents$SCT.res.2[n]
  obj$lymphoid.cell.type[which(obj$SCT_snn_res.2 == id)] = idents$ident[n]
}

DimPlot(obj,group.by = 'lymphoid.cell.type',raster = F,label = T)

DotPlot(obj,c('Cd8a','Cd4','Cd3d','Cd3e',
              'Col6a1','Krt8',
              'Itgam','Itgax','Cd14',
              'Ncr1','Nkg7',
              'Trdc','Il2ra','Il7r'),group.by = 'lymphoid.cell.type')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) + coord_flip()

obj <- obj %>% subset(lymphoid.cell.type %in% c('Epi contaminant','Myeloid contaminant') == F)



### re run the clustering

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL
## re normalize and cluster

obj <- SCTransform(obj,vars.to.regress = c("mitoRatio"),clip.range = c(-5,5))

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')

DimPlot(obj,group.by = 'lymphoid.cell.type')

## save process
#write_rds(obj,'Step2 Lymphoid cell/Tcell.rds')
FeaturePlot(obj,c('Cd3d','Cd3e','Cd4','Cd8a','Tcf7','Gzmb'))


