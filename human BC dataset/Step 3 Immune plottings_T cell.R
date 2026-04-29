library(tidyverse)
library(ggplot2)
library(Seurat)
library(readxl)
library(readr)
library(stringr)
library(scales)
library(splitstackshape)
options(future.globals.maxSize= 32000*1024^2)

######## Tumor bulk DEG analysis
obj <- read_rds('all_cells_all_type_cancer.rds')

t.cell <- obj  %>% subset(cellType == 'T_cell')

DefaultAssay(t.cell) <- 'RNA'

t.cell[['SCT']] <- NULL

t.cell <- NormalizeData(t.cell)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

t.cell <- CellCycleScoring(t.cell, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

t.cell <- SCTransform(t.cell,vars.to.regress = c('nFeature_RNA',
                                                 "S.Score", "G2M.Score",
                                                 'patient_id'),clip.range = c(-5,5))


t.cell <- t.cell %>%
  RunPCA(assay = "SCT", npcs = 20)%>%
  FindNeighbors( dims = 1:20, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')

Idents(t.cell) <- 'SCT_snn_res.0.5'

write_rds(t.cell,'t.cell.rds')

#### 2026.2.14.stopped

t.cell <- read_rds('t.cell.rds')

DimPlot(t.cell,label = T)+
DotPlot(t.cell,c('KLRB1','CD8A','CD4','CD3D'))


### Extract CD8 T cell

CD8.T <- t.cell %>% subset(SCT_snn_res.0.5 %in% c(1,4,5,9,10))

dt <- read_tsv('NEDD8_phenotype_label.txt')

CD8.T$N8_status <- lapply(CD8.T$patient_id,function(x) dt$N8_status[which(dt$patient_id == x)]) %>% unlist()





###  CD8T marker expression difference
###  
## Fig S2B

CD8.T$group <- paste(CD8.T$timepoint,CD8.T$expansion)
CD8.T$group <- factor(CD8.T$group,levels = c('Pre NE','Pre E','On NE','On E'))
DotPlot(CD8.T,group.by = 'group',
        list(Cytotoxic = c('GZMA','GZMB','GZMK','PRF1','IFNG','TNF'),
             Exhaustion = c('PDCD1','CTLA4','TOX','TIGIT','LAG3'),
             `IFN res` = c('STAT1','STAT3','ISG15','ISG20','IRF6','IRF8')))+
  #Naive = c('LEF1','SELL','TCF7','CCR7') )
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.4,'cm'))

###  # Fig 1D
CD8.T$group <- paste(CD8.T$timepoint,CD8.T$N8_status)
CD8.T$group <- factor(CD8.T$group,levels = c('On N8_high','On N8_low','Pre N8_high','Pre N8_low'))
DotPlot(CD8.T,group.by = 'group',
        list(Cytotoxic = c('GZMA','GZMB','GZMK','PRF1','IFNG','TNF'),
             Exhaustion = c('PDCD1','CTLA4','TOX','TIGIT','LAG3'),
             `IFN res` = c('STAT1','STAT3','ISG15','ISG20','IRF6','IRF8')))+
             #Naive = c('LEF1','SELL','TCF7','CCR7') )
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.4,'cm'))




################ CD8T cell marker patient average expression
#Fig S2a

genes <-    c(c('GZMA','GZMB','PRF1','IFNG'),
              c('PDCD1','CTLA4','TOX','TIGIT','LAG3'),
              c('STAT1','STAT3','ISG15','ISG20','IRF8') ) 

exp <- get_average_expression(CD8.T,
                              features =  genes ,
                              assay = 'SCT',
                              group.by = 'group',
                              individual_labels='patient_id',
                              reshape_for_ggplot = T,
                              remove_outliers_by_group.by = T)

exp$gene <- factor(exp$gene,levels =genes)

exp$group <- factor(exp$group,levels = c('Pre N8-high','Pre N8-low','On N8-high','On N8-low'))


library(ggh4x)
ggplot(exp,aes(x=group,y=value,colour = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = group))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45,hjust=-0.1),
        axis.text.x.bottom = element_blank(),
        axis.ticks.length.x = unit(0,'cm'))+
  scale_color_manual(values = c('#b20f0f','#b20f0f','#0F99B2','#0F99B2'))+
  scale_shape_manual(values = c(16, 17, 16, 17))+  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  facet_nested_wrap(~gene,ncol = 5,scales = 'free_y')+
  geom_signif(comparisons = list(c('Pre N8-high','Pre N8-low'),c('On N8-high','On N8-low')))
ylab('Normalized expression')



#### Fig 1C T cell DEG
#### 
Idents(CD8.T) <- 'group'
mks <- FindMarkers(CD8.T,ident.1 = 'Pre N8_low',ident.2 = 'Pre N8_high',)

FindMks_Volcano(mks,avg_lfc.hold = 0.4,p_adj.hold = 0.01)

mks <- mks %>% filter(p_val_adj < 0.01 & pct.1 >0.05 & pct.2 >0.05)

highlights <- read_tsv("DEG/cd8/highlight_key.txt.txt",col_names = F) %>% unlist

FindMks_Volcano(mks,avg_lfc.hold = 0.4,p_adj.hold = 0.01,top_n_plot = 0,
                gene.highlight =highlights)+ylim(0,100)




