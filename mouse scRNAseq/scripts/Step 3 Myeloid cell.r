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

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL
## subset cluster
obj <- obj %>% subset(celltype.general.1 =='Macrophage')

## re normalize and cluster
obj <- SCTransform(obj, 
                   vars.to.regress = c("mitoRatio"),
                   clip.range = c(-5,5),
                   verbose = TRUE)


obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:20, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')

Idents(obj) <- 'SCT_snn_res.0.8'
DimPlot(obj,label = T)

FeaturePlot(obj,'Ptprc')

FeaturePlot(obj,c('Krt8','Cd3e','Cd3d','Cd8a','Cd4','Ncr1','Itgam'))

DotPlot(obj,c('Krt8','Cd3e','Cd3d','Cd8a','Cd4','Ncr1','Itgam'))

## remove the Krt8
obj <- obj %>% subset(SCT_snn_res.0.8 != 12)

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')

Idents(obj) <- 'SCT_snn_res.2'

## save process
#write_rds(obj,'Step3 Myeloid cell/Macro.rds')

# c <- FindMarkers(obj,ident.1 = 6)

obj <- read_rds('./Step3 Myeloid cell/Macrophage/Macro.rds')

obj$group <-factor(obj$group,
                   levels =c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))

DotPlot(obj,c('Krt8','Itga6','Krt19','Cd3e','Cd3d','Cd8a','Cd4','Ncr1','Itgam','klrd1c','Cd14'))

DimPlot(obj,label = T)
obj <- obj %>% subset(SCT_snn_res.2 != '30')

FeaturePlot(obj,c('Krt8','Cd3e','Cd3d','Ncr1','Itgam','klrb1c','Cd14'))

## function signature score

markers<- read_excel('Step3 Myeloid cell/Macrophage/macrophage markers.xlsx',sheet = 1)

obj <- calculate_sig(obj,markers,h2m_conversion = F,split.by = ',')

signature <- obj[[2]]
obj<-obj[[1]]


avg_heatmap(obj,signature%>%names(),group_by = 'group')
# DotPlot(obj,IFN_res,group.by = 'group')+
#   scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))



macro  <- list(
Chemotaxis = c('Ccl2','Ccl7','Cxcl2','Cxcl8','Cxcl9','Cxcl10'),
# M1 = c('Nos2','Macro','Tlr2','Cd80','Cd86',
#        'Csf2','Tnf','Il1b','Il6','Tlr4','Cxcl2','Ifng','Il1r1')
# 
# M2 = c('Mrc1','Pparg','Arg1','Cd163','Clec10a',
#        'Clec7a','Pdcd1lg2','Retnla','Ccl22','Il10','Fcrg3a',
#        'Il4','Irf4','Pdgfb','Stat6','Chil3')

lipid = c('Apoe','Lipa','Fabp4','Fabp5','Trem2','Lpl'),
#Trem2, Lipa, Lpl, Ctsb, Ctsl, Fabp4, Fabp5, Lgals1, Lgals3, Cd9, and Cd36
#anti_tumor <- c('Foxp1','Tnf','Il6','Il1b','Nlrp3'),

#M1 = c('Ly6c2','Tnf','Stat1','Nos2','Cd86','Il12a','Il6','Ccr2'),
#M2 = c('Tgfb1','Il13ra1','Prkaa1','Clec10a','Il4ra','Arg1','Retnla','Cxc3cr1','Il10','Cd163','Igf1','Mrc1'),

M2 = c('Selenop','Slc40a1','Pltp','F13a1','Fuca2','Arg1','Mrc1','Cd163','Tgfb1','Retnlg'), #'Ccl20'

M1 = c('Pld4','Il1b','Foxp1','Tnf','Tlr4','Cd80','Cd86','Nos2'),

IFN_res = c('Isg20','Isg15','Ifnar1','Stat1','Irf1','Irf2','Irf8'),
exhaustion = c('Entpd1','Havcr2','Pdcd1lg2','Cd274'))


DotPlot(obj,macro,group.by = 'group')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  theme(axis.text=element_text(angle = -90),axis.text.y = element_text(hjust = 0.5))
  

DotPlot(obj,macro,group.by = 'sample_id')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  geom_hline(yintercept = c(4.5,8.5,11.5),linetype='dashed')+
  theme(axis.text.x=element_text(angle = -90),axis.text.y = element_text(hjust = 0.5))



###
FeaturePlot(obj[,sample(colnames(obj),10000)],names(signature))&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

VlnPlot(obj,'REACTOME_FATTY_ACID_METABOLISM',group.by = 'group',pt.size = 0)
VlnPlot(obj,names(signature),group.by = 'group',pt.size = 0)

### each group downsmaple 10k cells
sub_obj <- lapply(unique(obj$group), function(x) obj%>% subset(group==x) %>% colnames() %>% sample(10000)) %>% unlist() 
sub_obj <- obj[,sub_obj]


FeaturePlot(sub_obj,names(signature),split.by = 'group')&
  scale_colour_gradientn(colors = c('#5dade2','#EEEEEF','#e25d5d'))

FeaturePlot(sub_obj,names(signature))&
  scale_colour_gradientn(colors = c('#5dade2','#EEEEEF','#e25d5d'))

### avg signature score of all samples ( each sample down to 2000 cells)

sub_obj2 <- lapply(unique(obj$sample_id), function(x) obj%>% subset(sample_id==x) %>% colnames() %>% sample(2000)) %>% unlist() 
sub_obj2 <- obj[,sub_obj2]



avg_heatmap(obj ,names(signature),
            group_by = 'group',
            label = F,
            x.order = c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))+
  theme_bw()

DotPlot(obj,c('Arg1','Vegfa','Tgfb1','Nos2','Tnf','Ifng','Apoe'),group.by = 'group')
FeaturePlot(obj,c('Ly6c1','Ly6g','H2-Aa','H2-Eb1','H2-Ab1'))


avg_heatmap(obj ,names(signature),
            group_by = 'sample_id',
            label = F)+
  theme_bw()+
  geom_vline(xintercept = 4.5)+
  geom_vline(xintercept = 8.5)+
  geom_vline(xintercept = 11.5)

VlnPlot(obj,c('Apoc1', 'Apoe', 'Acp5','Fabp5'),group.by = 'group')




table(obj[['sample_id']])
DotPlot(obj,names(signature),group.by = 'group')


## DEG
Idents(sub_obj2) <- 'group'
DEG.IgG <- FindMarkers(sub_obj2,ident.1 = 'N8KO_IgG',ident.2 = 'WT_IgG')
DEG2RNK(DEG.IgG,log2fc_hold = 0.1,p_hold = 0.1,name = 'Macro_DEG_N8KO_IgG-vs-WT_IgG')

DEG.aPD1 <- FindMarkers(sub_obj2,ident.1 = 'N8KO_aPD1',ident.2 = 'WT_aPD1')
DEG2RNK(DEG.aPD1,log2fc_hold = 0.1,p_hold = 0.1,name = 'Macro_DEG_N8KO_aPD1-vs-WT_aPD1')


## before GSEA check if you have run easyBuild()
## Use the full path


GSEA_batch.rnk(
    GSEA_installed_path = "C:/bioinfo/irineos_mice_scRNA_final_version/GSEA_4.3.2",
    DEG_path='C:/bioinfo/irineos_mice_scRNA_final_version/Step3 Myeloid cell/Macrophage/GSEA/RNKs',
    species = 'mouse',
    gene_sets =c(`hallmark gene sets`='mh.all.v2024.1.Mm.symbols.gmt',
                 `reactome gene sets` = 'm2.cp.reactome.v2024.1.Mm.symbols.gmt',
                 `GOBP gene sets`='m5.go.bp.v2024.1.Mm.symbols.gmt',
                 `wiki pathways` = 'm2.cp.reactome.v2024.1.Mm.symbols.gmt'),
    symbol_chip='Mouse_Gene_Symbol_Remapping_MSigDB.v2024.1.Mm.chip',
    out_dir='auto',
    GSEA_plots_number=30,
    collapse='Remap_Only'
)

### GSEA bubble plot
GSEA_bubble('C:/bioinfo/irineos_mice_scRNA_final_version/Step3 Myeloid cell/Macrophage/GSEA',GSEA_fdr_hold = 1)
  




##------------------------------ end 


### dendritic cell
obj <- read_rds('1.general.ident.assigned.rds')

## subset cluster
obj <- obj %>% subset(celltype.general.1 =='DCs')

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL


## re normalize and cluster
obj <- SCTransform(obj, 
                   vars.to.regress = c("mitoRatio"),
                   clip.range = c(-5,5),
                   verbose = TRUE)


obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:20, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')


DimPlot(obj,label = T)


FeaturePlot(obj,c('C1qa','C1qb','Krt8','Cd3d','Cd3e',
                  'Nkg7','Ncr1','Adgre1',
                  'Col1a1'))

obj <- obj%>%subset(SCT_snn_res.2 %in% c( 19,21) == F)


table(obj[['sample_id']])

pca <- prcomp(t(log1p(assays(sce)$counts)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

## DEG
Idents(sub_obj2) <- 'group'
DEG.IgG <- FindMarkers(sub_obj2,ident.1 = 'N8KO_IgG',ident.2 = 'WT_IgG')
DEG2RNK(DEG.IgG,log2fc_hold = 0.1,p_hold = 0.1,name = 'DCs_DEG_N8KO_IgG-vs-WT_IgG.rnk')

DEG.aPD1 <- FindMarkers(sub_obj2,ident.1 = 'N8KO_aPD1',ident.2 = 'WT_aPD1')
DEG2RNK(DEG.aPD1,log2fc_hold = 0.1,p_hold = 0.1,name = 'DCs_DEG_N8KO_aPD1-vs-WT_aPD1.rnk')

# unable to retrieve enough DEG

write_rds(obj,'DC.rds')

