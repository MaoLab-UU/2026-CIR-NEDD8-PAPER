library(dplyr)
library(Seurat)
library(readxl)
library(readr)
library(stringr)
library(scales)

lapply(c(list.files("../R tools",full.names = T,pattern = '.R'),
         list.files("../R tools",full.names = T,pattern = '.r')),source) %>% invisible()

options(future.globals.maxSize= 32000*1024^2)

#### create obj

coh1 <- readRDS('1863-counts_cells_cohort1.rds')
meta1 <- read_delim('1872-BIOKEY_metaData_cohort1_web.csv',
                    delim = ',') %>% as.data.frame()
obj <- CreateSeuratObject(coh1,meta.data = meta1,
                          min.cells=3,min.features = 200)

#tnbc <- subset(obj,BC_type == 'TNBC' & expansion != 'n/a' & timepoint=='Pre')
obj <- subset(obj,expansion != 'n/a')
rm(coh1)
rm(meta1)
rm(requirements)

#### normalization

obj <- QC(obj,species = 'hs')

obj <- QC_plot(obj,filtering = T)

obj <- NormalizeData(obj)
obj <- cell_cycle_score(obj,species = 'hs')

### SCTransform

obj <- SCTransform(obj,
                   vars.to.regress = c(
                                       "nFeature_RNA", 
                                       "HbRatio",
                                       "mitoRatio", 
                                       "HspRatio",
                                       'cc.difference'),
                   verbose = TRUE)



obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:20, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')

write_rds(obj,'all_cells_all_type_cancer.rds')
# Fig 1.A
DimPlot(obj,group.by = 'cellType')




