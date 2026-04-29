library(tidyverse)
library(readxl)
library(Seurat)
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
library(SingleR)
library(scp)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 32000*1024^2)
lapply(list.files("C:/bioinfo/R tools",full.names = T),source) %>% invisible()

### read object
obj <- read_rds('all_cells_raw.rds')
table(obj$sample_id)

obj$group <- str_replace_all(obj$group,'_',' ')
obj$group <- str_replace_all(obj$group,'a','α')
obj$group <- factor(obj$group,levels = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'))


### check cell cycle difference
DefaultAssay(obj) <- 'RNA'
obj <- obj %>%FindVariableFeatures() %>% ScaleData()%>% RunPCA()

DimPlot(obj,split.by = 'Phase')

FeaturePlot(obj,'cc.difference')&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

### SCTransform
obj <- SCTransform(obj,
                     vars.to.regress = c("nCount_RNA",
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

### check sample variance
DimPlot(obj,split.by = 'sample_id',ncol = 5,raster=FALSE,group.by = 'SCT_snn_res.0.8')
DimPlot(obj,split.by = 'treatment',raster=FALSE,group.by = 'SCT_snn_res.0.8')
DimPlot(obj,split.by = 'genotype',raster=FALSE,group.by = 'SCT_snn_res.0.8')

DimPlot(obj,raster=FALSE,label = T)



#write_rds(obj,'all_cells_sctransformed_regressed.rds')

### 1.calculate signature score
# read obj
#obj <- readRDS('all_cells_sctransformed_regressed.rds')



# read marker file
Markers <- read_xlsx('Step1 cell general identity/Markers_expressed.xlsx',skip = 1) 

# calculate signature
obj <- calculate_sig(obj,Markers,h2m_conversion=F)
signatures.mm <- obj[[2]]
obj <- obj[[1]]



### general signature score plot
p<- FeaturePlot(obj[,sample(obj%>%colnames,10000)],names(signatures.mm),raster = F)&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

pdf('General Signatures.pdf',
    width = length(names(signatures.mm)) ^ 0.5 %>% floor()*6,
    height = length(names(signatures.mm)) ^ 0.5 %>% floor()*6)
print(p)
dev.off()

### plot genes of each cell type signature set
for(ct in names(signatures.mm)){
  plot_genes <- signatures.mm[[ct]]
  p<-FeaturePlot(obj[,sample(obj%>%colnames,10000)],plot_genes,raster = F,
                 ncol = length(plot_genes) ^ 0.5 %>% floor())&
    scale_color_gradientn(colors = c('#EEEEEF','#E0191B'))
  
  pdf(paste0(ct,'_detailed_markers.pdf'),
      width = length(plot_genes) ^ 0.5 %>% floor()*6,
      height = length(plot_genes) ^ 0.5 %>% ceiling()*6)
  print(p)
  dev.off()
}

# dotplot 
DotPlot(obj,colnames(Markers),group.by =  'SCT_snn_res.2')+
  scale_color_gradientn(colors = c('#EEEEEF','#E0191B'))+
  coord_flip()

### 2.assign cell identity

idents <- read_excel('Step1 cell general identity/step1 general cell identity.xlsx')

obj$celltype.general.1 <- 'NA'

for(n in 1:nrow(idents)){
  id <- idents$SCT.res.2[n]
  obj$celltype.general.1[which(obj$SCT_snn_res.2 == id)] = idents$ident[n]
}

DimPlot(obj,group.by = 'celltype.general.1',raster = F,label = T)



plot_names <- colnames(Markers)[! colnames(Markers) %in%  c('Neutrophil','B cell') ]
avg_heatmap(obj,plot_names,group_by = 'celltype.general.1')+
  ggtitle('y axis cell signature score, x axis cell type')+
  theme(axis.text.x = element_text(angle=90, hjust=1))


###### 3.cell typist on immune cells
## only need to run once to bulid environment
# easyBuild()

model_path  <- '../Celltypist_models'
use_condaenv('celltypist')
celltypist = import('celltypist') 
scanpy= import('scanpy')
pandas= import('pandas')
numpy= import('numpy')

#Model_download(model_path )           ######### install models, only need to run once
#Model_list(model_path)                ######### list available models


sub_obj <- obj %>% subset(celltype.general.1 %in% c('Epithelial cell','Fibroblast','Endothelial cell') == F) 
DefaultAssay(sub_obj) <- 'RNA'
sub_obj[['SCT']] <- NULL
sub_obj <- sub_obj[,sample(colnames(sub_obj),30000)]
sub_obj <- NormalizeData(sub_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#rm(obj)
predicted_data <- Runcelltypist(sub_obj,model='Adult_Mouse_Gut',majority_voting = F)

# d1 <- DimPlot(predicted_data,group.by = 'celltype.general.1')
# d2 <- DimPlot(predicted_data,group.by = "typist_prediction")
# d1+d2

### typist prediction percentage, ignore rare cell sub population in plots
predicted.cell.pct <- table(predicted_data$typist_prediction) %>% prop.table() * 100 
plot.cell.types <- names(predicted.cell.pct)[predicted.cell.pct>0.5]
predicted_data <- predicted_data %>% subset(typist_prediction %in% plot.cell.types)

plot.list <- lapply(predicted_data$typist_prediction %>% unique(),
                    function(x){ predicted_data %>% 
                        DimPlot(cells.highlight = colnames(predicted_data)[predicted_data$typist_prediction == x]) + 
                        # xlim(c(-15,15)) + 
                        # ylim(c(-17,17)) +
                        ggtitle(x)+
                        theme(legend.position="none")} )

do.call("grid.arrange", c(plot.list, 
                          ncol= length(plot.cell.types) ^ 0.5 %>% ceiling(),
                          top='Immune cell Typist prediction (cell population with pct > 0.5%)'))



### save progress
write_rds(obj,'1.general.ident.assigned.rds')




# ### check IRAK3 expression
obj<-readRDS('1.general.ident.assigned.rds')

DotPlot(obj,list(Epithelial = c('Krt8','Itga6'),
                 `T cell` = c('Cd3d','Cd3e','Cd4','Cd8a'),
                 NK = c('Klrb1c','Ncr1'),
                 DCs = c('Batf3','Xcr1','Itgax','Clec9a'),
                 `Macrophage/Monocyte` = c('C1qa','Adgre1','Lyz2','Ly6c1'),
                 Neutrophil = c('S100a8','S100a9','Ly6g'),
                 Endothelial = c('Pecam1','Cdh5'),
                 Fibroblast= c('Col1a2','Col1a1','Dcn'),
                 `Mast cell` = c('Tpsab1','Mcpt4')),group.by = 'celltype.general.1')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))


DimPlot(obj,group.by = 'celltype.general.1',raster = F)


## general cell percentage in tumor

meta <- read_xlsx('metadata.xlsx')
meta$sample_id <- meta$Sample
meta$group <- paste0(meta$Genotype,'_',meta$Treatment)

group_tb <- obj@meta.data[c('celltype.general.1','sample_id')] %>% group_by(sample_id)%>%
  table() %>% t() %>% as.data.frame()
group_tb <- merge(group_tb,meta,by='sample_id') 
group_tb$Freq <- group_tb$Freq/group_tb$`Total cells after quality control` 
write_csv2(group_tb,'T_percent_in_total_tumor.csv')

### make box plot
group_tb <- group_tb[,c('sample_id','group','Freq','celltype.general.1')]
plot_dt <- group_tb %>% melt()
box.plot.list <- lapply( split(plot_dt,plot_dt[['celltype.general.1']]), function(x) 
  ggplot(x,aes(x= factor(group,
                         levels =c('N8KO_aPD1','WT_aPD1','N8KO_IgG','WT_IgG') ),
               y=value ,color=group))+
    geom_boxplot(outliers = F)+
    geom_point(position=position_jitter())+
    theme(axis.text.x=element_text(angle=90))+
    geom_signif(comparisons = list(c('N8KO_aPD1','WT_aPD1'),c('N8KO_IgG','WT_IgG')),
                y_position = c(max(x$value)*1.05,max(x$value)*1.1),color='black')+
    theme_bw()+
    ylab('Cell Frequency')+
    xlab('Group')+
    ylim(0,max(x$value)*1.2))


legend_plot.list(box.plot.list,'box',name = 'general cell pct',width.factor = 5,height.factor = 3)

