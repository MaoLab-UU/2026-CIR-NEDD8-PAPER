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
library(colorspace)
library(ggridges)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 32000*1024^2)
lapply(list.files("C:/bioinfo/R tools",full.names = T),source) %>% invisible()

##
obj <- read_rds('Step2 Lymphoid cell/Tcell.rds')

obj$group <- factor(obj$group,
                    levels =c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))

## 1.cd4 cd8 percentage in T cell population
{
group_tb <- obj@meta.data[c('lymphoid.cell.type','sample_id')] %>% group_by(sample_id)%>%
  table() %>% prop.table(margin=2) %>% t() %>% as.data.frame()

meta <- read_xlsx('metadata.xlsx')
meta$sample_id <- meta$Sample
meta$group <- paste0(meta$Genotype,'_',meta$Treatment)



group_tb <- merge(group_tb,meta,by='sample_id') 

write_csv2(group_tb,'T_percent_in_T cell population.csv')

  ### make box plot
group_tb <- group_tb[,c('sample_id','group','Freq','lymphoid.cell.type')]
plot_dt <- group_tb %>% melt()
box.plot.list <- lapply( split(plot_dt,plot_dt[['lymphoid.cell.type']]), function(x) 
  ggplot(x,aes(x= factor(group,
                         levels =c('N8KO_aPD1','WT_aPD1','N8KO_IgG','WT_IgG') ),
               y=value,color=group))+
    geom_boxplot(outliers = F)+
    geom_point(position=position_jitter())+
    theme(axis.text.x=element_text(angle=90))+
    geom_signif(comparisons = list(c('N8KO_aPD1','WT_aPD1'),c('N8KO_IgG','WT_IgG')),
                y_position = c(max(x$value)*1.05,max(x$value)*1.1),color='black')+
    theme_bw()+
    ylab('Cell Frequency')+
    xlab('Group')+
    ylim(0,max(x$value)*1.2))

name <- 'T_percent_in_T cell population'
legend_plot.list(box.plot.list,'box',name = name,width.factor = 5,height.factor = 3 )
}

## 2.cd4 cd8 percentage in total tumor 
{
meta <- read_xlsx('metadata.xlsx')
meta$sample_id <- meta$Sample
meta$group <- paste0(meta$Genotype,'_',meta$Treatment)

group_tb <- obj@meta.data[c('lymphoid.cell.type','sample_id')] %>% group_by(sample_id)%>%
  table() %>% t() %>% as.data.frame()
group_tb <- merge(group_tb,meta,by='sample_id') 
group_tb$Freq <- group_tb$Freq/group_tb$`Total cells after quality control` 
write_csv2(group_tb,'T_percent_in_total_tumor.csv')

  ### make box plot
group_tb <- group_tb[,c('sample_id','group','Freq','lymphoid.cell.type')]
plot_dt <- group_tb %>% melt()
box.plot.list <- lapply( split(plot_dt,plot_dt[['lymphoid.cell.type']]), function(x) 
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

name <- 'T_percent_in_total_tumor'
legend_plot.list(box.plot.list,'box',name = name,width.factor = 5,height.factor = 3)

}

##########
DimPlot(obj,group.by = 'lymphoid.cell.type')
cd4 <- obj %>% subset(lymphoid.cell.type == 'CD4 T')
cd8 <- obj %>% subset(lymphoid.cell.type == 'CD8 T')
rm(obj)
# write_rds(cd4,'Step2 Lymphoid cell/CD4/CD4.rds')
# write_rds(cd8,'Step2 Lymphoid cell/CD8/CD8.rds')


####
####  CD8 analysis
obj <- read_rds('Step2 Lymphoid cell/CD8/cd8.rds')


DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

## re normalize and cluster

obj <- SCTransform(obj,vars.to.regress = c("mitoRatio"), clip.range = c(-5,5))

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) 

DimPlot(obj,label = T,group.by = 'SCT_snn_res.0.8')
DimPlot(obj,label = T,group.by = 'SCT_snn_res.2')
DimPlot(obj,group.by = 'SCT_snn_res.0.8',split.by = 'sample_id',ncol = 5)

FeaturePlot(obj,c('Cd8a','Cd8b1','Cd3d','Cd3e'))

### Add cell markers
# read marker file
Markers <- read_xlsx('Step2 Lymphoid cell/CD8/Makers_human.xlsx',skip = 1) 

# calculate signature
obj <- calculate_sig(obj,Markers,h2m_conversion=T)
signatures.mm <- obj[[2]]
obj <- obj[[1]]


FeaturePlot(obj,signatures.mm%>%names())&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

# DotPlot(obj,signatures.mm$`Cytokine/Cytokine receptor`)&
#    scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))


avg_heatmap(obj,signatures.mm%>%names(),group_by = 'SCT_snn_res.0.8',label = F)
avg_heatmap(obj,signatures.mm%>%names(),group_by = 'SCT_snn_res.2')


avg_heatmap(obj,signatures.mm%>%names(),group_by = 'group',x.order = c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

avg_heatmap(obj,signatures.mm%>%names(),group_by = 'sample_id',label = F)+  theme_bw()+
  geom_vline(xintercept = 4.5)+
  geom_vline(xintercept = 8.5)+
  geom_vline(xintercept = 11.5)



Idents(obj) <- 'SCT_snn_res.0.8'

### marker dot plot
DotPlot(obj,list(
              naive= c('Ccr7','Tcf7','Il7r','Sell','Lef1'), # naive
              `Effector/Actication` = c('Gzmk','Gzmb','Prf1','Cd44','Cd69','Ifng','Eomes','Il2ra','Il2rb'),# effector
              #M = c('Klrd1','Klrb1','Prdm1','Cxcr3','Cxcr4'),
              `TR` = c('Itga1','Itgae'), # memory
              Exhasution = c('Pdcd1','Lag3','Havcr2','Ctla4'), ## exhasution
              Proliferation = c('Mki67','Tuba1b','Top2a') ### proliferation
              ),dot.scale = 5)+ 
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))



DotPlot(obj,list(
  naive= c('Ccr7','Tcf7','Il7r','Sell','Lef1'), # naive
  `Effector/Actication` = c('Gzmk','Gzmb','Prf1','Cd44','Cd69','Ifng','Eomes','Il2ra','Il2rb'),# effector
  #M = c('Klrd1','Klrb1','Prdm1','Cxcr3','Cxcr4'),
  `TR` = c('Itga1','Itgae'), # memory
  Exhasution = c('Pdcd1','Lag3','Havcr2','Ctla4'), ## exhasution
  Proliferation = c('Mki67','Tuba1b','Top2a') ### proliferation
),dot.scale = 5,group.by = 'group')+ 
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))




DotPlot(obj,list(
  naive= c('Ccr7','Tcf7','Il7r','Sell','Lef1'), # naive
  `Effector/Actication` = c('Gzmk','Gzmb','Prf1','Cd44','Cd69','Ifng','Eomes','Il2ra','Il2rb'),# effector
  #M = c('Klrd1','Klrb1','Prdm1','Cxcr3','Cxcr4'),
  Idents = c('Cd3d','Cd3e','Cd4','Cd8a','C8b1'),
  `TR` = c('Itga1','Itgae'), # memory
  Exhasution = c('Pdcd1','Lag3','Havcr2','Ctla4'), ## exhasution
  TCR = c(),
  Proliferation = c('Mki67','Tuba1b','Top2a') ### proliferation
),dot.scale = 5,group.by = 'sample_id')+ 
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  geom_hline(yintercept = c(4.5,8.5,11.5),linetype='dashed')


### slingshot lineage

dim<- 'SCT_snn_res.0.8'
sce <- as.SingleCellExperiment(obj,assay = 'SCT')
sce <- slingshot(sce, clusterLabels = dim, reducedDim = "PCA",
                 allow.breaks = T)

#lnes <- getLineages(reducedDim(sce,"PCA"),sce[[dim]], start.clus = "11")  

sce.new <- embedCurves(sce,obj@reductions$umap@cell.embeddings)


### lineage plot
plot_dt <- obj@reductions$umap@cell.embeddings
plot_dt <- cbind(plot_dt,obj[[dim]]) %>% as.data.frame()
colnames(plot_dt)[3] <- 'SCT.cluster'
  
colors <- plot_dt$SCT.cluster %>% unique() %>% length() %>% rainbow()
plot_dt$color <-  colors[plot_dt$SCT.cluster] 

legend <- unique(plot_dt$SCT.cluster)[order(unique(plot_dt$SCT.cluster))]
fill   <-   unique(plot_dt$color)[order(legend)]

c <- SlingshotDataSet(sce.new)
n_line <- c@lineages %>% length()


plot(x=plot_dt$umap_1,y=plot_dt$umap_2, type = "p",col = NULL, bg = alpha('black',0.2),pch = 21,cex = 0.4)+
lines(SlingshotDataSet(sce.new), lwd = 3, col = viridis_pal()(n_line)) +
legend("bottomright",
           legend = c@lineages %>% names(),
           fill = viridis_pal()(n_line),
           title = paste0('trajectories ',dim),
           cex = 0.5,
           pt.cex = 0.5)





### avg time 

avg.ptime <- averagePseudotime(sce %>% slingPseudotime)

obj@meta.data$avg.ptime <- avg.ptime %>% scale()

FeaturePlot(obj,'avg.ptime')&
  scale_color_gradientn(colors = viridis_pal()(30))

###individual lineage plot

lin_plot.list <- lapply(1:length(c@curves), function(lin){
  obj[[paste0('slingPseudotime_',lin)]] <- sce[[paste0('slingPseudotime_',lin)]]
  
  cell.names <- colnames(sce)[which(sce[[paste0('slingPseudotime_',lin)]] != 'NA')]

  return(FeaturePlot(obj,paste0('slingPseudotime_',lin),cells = cell.names)&
    scale_color_gradientn(colors = viridis_pal()(30)))
  
})

legend_plot.list(lin_plot.list,extension = '',name='CD8 trajectories each plot individually')

### Assign subsets

idents <- read_excel('Step2 Lymphoid cell/CD8/CD8T subcluster.xlsx')

obj$CD8.subset <- 'NA'

for(n in 1:nrow(idents)){
  id <- idents$SCT.res.0.8[n]
  obj$CD8.subset[which(obj$SCT_snn_res.0.8 == id)] = idents$ident[n]
}

DimPlot(obj,group.by = 'CD8.subset',raster = F,label = T)+theme(legend.position = 'none')

DimPlot(obj,split.by = 'sample_id',raster = F,label = F,ncol=5)

### plot markers and signature for each cell type

avg_heatmap(obj,signatures.mm%>%names(),group_by = 'CD8.subset')+
  theme_bw ()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

DotPlot(obj,list(
  naive= c('Ccr7','Tcf7','Il7r','Sell','Lef1'), # naive
  `Effector/Actication` = c('Gzmk','Gzmb','Prf1','Cd44','Cd69','Ifng','Eomes','Il2ra','Il2rb'),# effector
  #M = c('Klrd1','Klrb1','Prdm1','Cxcr3','Cxcr4'),
  `TR` = c('Itga1','Itgae'), # memory
  Exhasution = c('Pdcd1','Lag3','Havcr2','Ctla4'), ## exhasution
  Proliferation = c('Mki67','Tuba1b','Top2a') ### proliferation
  
),dot.scale = 5,group.by = 'sample_id')+ 
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))


FeaturePlot(obj,c('Exhaustion','Cytotoxicity',
                  'Naïve','Activation:Effector function',
                  'IFN Response','Senescence','Stress response'))&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

### cell subset percentage

######### percentage plot in CD8 subset
c<-table(obj@meta.data[,c('CD8.subset','sample_id')]) %>% prop.table(margin = 2) %>% as.data.frame()
meta <- read_xlsx("metadata.xlsx")
c$group <- NA
for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
}
c$Freq <- c$Freq* 100

pct_plot(c)


# ggplot(c,aes(x=CD8.subset,y=Freq* 100,colour=group))+
#   geom_boxplot(outliers = F)+
#   geom_point(position = position_jitterdodge(jitter.width=0.1) )+
#   stat_signif(comparisons = list(c('N8KO aPD1','N8KO IgG')))+
#   theme_bw()+
#   ylab('Cell percent %')

######### percentage plot in total tumor cells

c<-table(obj@meta.data[,c('CD8.subset','sample_id')]) %>% as.data.frame()

meta <- read_xlsx("metadata.xlsx")
c$group <- NA
c$total_number <- NA

for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
  c$total_number[i] <- meta$`Total cells after quality control`[which(meta$Sample==c$sample_id[i])]
}
c$Freq <- c$Freq/c$total_number * 100 

pct_plot(c)

# 
# ggplot(c,aes(x=CD8.subset,y=(Freq/total_number * 100),colour=group))+
#   geom_boxplot(outliers = F)+
#   geom_point(position = position_jitterdodge(jitter.width=0.1) )+
#   theme_bw()+ 
#   ylab('Cell percent %')




write_rds(obj,'Step2 Lymphoid cell/CD8/cd8.rds')










###########
########### process cd4 T cell
obj <- read_rds('Step2 Lymphoid cell/CD4/CD4.rds')

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

## re normalize and cluster

obj <- SCTransform(obj,vars.to.regress = c("mitoRatio"), clip.range = c(-5,5))

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) 

Idents(obj) <- 'SCT_snn_res.0.8'
DimPlot(obj,label = T,group.by = 'SCT_snn_res.0.8')+
DimPlot(obj,label = T,group.by = 'SCT_snn_res.2')
DimPlot(obj,group.by = 'SCT_snn_res.0.8',split.by = 'sample_id',ncol = 5)


Idents(obj) <- 'SCT_snn_res.2'

DotPlot(obj,c('Cd4','Cd3e','Cd3d','Il2ra','Foxp3','Pdcd1','Cd8a'))&
  scale_color_gradientn(colors = c('#EEEEEF','#E0191B'))

DotPlot(obj,c('Nkg7','Prf1'))

### remove SCT_snn_res.2 cluster 15 it is CD8
obj <- obj %>% subset(SCT_snn_res.2 != 15)

FeaturePlot(obj,c('Cd4','Cd3e','Cd3d','Il2ra','Foxp3','Pdcd1','Cd8a','Cd8b1'))&
  scale_color_gradientn(colors = c('#EEEEEF','#E0191B'))


Markers <- read_xlsx('Step2 Lymphoid cell/CD4/Markers_human.xlsx',skip = 1) 

# calculate signature
obj <- calculate_sig(obj,Markers,h2m_conversion=T)
signatures.mm <- obj[[2]]
obj <- obj[[1]]


FeaturePlot(obj,signatures.mm%>%names())&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

DotPlot(obj,signatures.mm$`Cytokine/Cytokine receptor`)&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))



avg_heatmap(obj,signatures.mm%>%names(),group_by = 'SCT_snn_res.2')


avg_heatmap(obj,signatures.mm%>%names(),group_by = 'group',x.order = c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

avg_heatmap(obj,signatures.mm%>%names(),group_by = 'sample_id',label = F)+  theme_bw()+
  geom_vline(xintercept = 4.5)+
  geom_vline(xintercept = 8.5)+
  geom_vline(xintercept = 11.5)






#### sub cluster T helper

Idents(obj) <- 'SCT_snn_res.2'


DotPlot(obj ,list(
  # Thf
  Thf = c('Bcl6','Cxcr5','Ica1','Gng4','Nmb','Ebi3','Il21'),
  #Th1
  Th1 = c('Ifng','Tbx21','Gzmb','Nr4a1','Ccl3','Ccl4'),
  #Th17
  Th17 = c('Il17a','Il17f','Rora',#'Klrb1',
           'Ccr6','Ccr4','Irf4','Batf'),
  # memory
  memory = c('Il7r','Cd69','Gpr183')
  
  ),group.by = 'SCT_snn_res.2')&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))





### slingshot lineage for Treg

sub_obj <- obj %>% subset(SCT_snn_res.2 %in% c(0:3,5,7,8,11:14))
# 
# DefaultAssay(sub_obj) <- 'RNA'
# sub_obj[['SCT']] <- NULL
# 
# sub_obj  <- SCTransform(sub_obj ,vars.to.regress = c("mitoRatio"), clip.range = c(-5,5))%>%
#   RunPCA(assay = "SCT", npcs = 50)%>%
#   RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')%>%
#   FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
#   FindClusters(resolution = c(0.8,2)) 
# 
# DimPlot(sub_obj)



dim<- 'SCT_snn_res.2'
sce <- as.SingleCellExperiment(sub_obj,assay = 'SCT')
sce <- slingshot(sce, clusterLabels = dim, reducedDim = "PCA",
                 allow.breaks = T)

#lnes <- getLineages(reducedDim(sce,"PCA"),sce[[dim]], start.clus = "11")  

sce.new <- embedCurves(sce,sub_obj@reductions$umap@cell.embeddings)


### lineage plot
plot_dt <- sub_obj@reductions$umap@cell.embeddings
plot_dt <- cbind(plot_dt,sub_obj[[dim]]) %>% as.data.frame()
colnames(plot_dt)[3] <- 'SCT.cluster'

colors <- plot_dt$SCT.cluster %>% unique() %>% length() %>% rainbow()
plot_dt$color <-  colors[plot_dt$SCT.cluster] 

legend <- unique(plot_dt$SCT.cluster)[order(unique(plot_dt$SCT.cluster))]
fill   <-   unique(plot_dt$color)[order(legend)]

c <- SlingshotDataSet(sce.new)
n_line <- c@lineages %>% length()


plot(x=plot_dt$umap_1,y=plot_dt$umap_2, type = "p",col = NULL, bg = alpha('black',0.2),pch = 21,cex = 0.4)+
  lines(SlingshotDataSet(sce.new), lwd = 3, col = viridis_pal()(n_line)) +
  legend("topright",
         legend = c@lineages %>% names(),
         fill = viridis_pal()(n_line),
         title = paste0('trajectories ',dim),
         cex = 0.5,
         pt.cex = 0.5)





### avg time 

avg.ptime <- averagePseudotime(sce %>% slingPseudotime)

sub_obj@meta.data$avg.ptime <- -avg.ptime %>% scale()

FeaturePlot(sub_obj,'avg.ptime')&
  scale_color_gradientn(colors = viridis_pal()(30))

###individual lineage plot

lin_plot.list <- lapply(1:length(c@curves), function(lin){
  sub_obj[[paste0('slingPseudotime_',lin)]] <- -sce[[paste0('slingPseudotime_',lin)]]
  
  cell.names <- colnames(sce)[which(sce[[paste0('slingPseudotime_',lin)]] != 'NA')]
  
  return(FeaturePlot(sub_obj,paste0('slingPseudotime_',lin),cells = cell.names)&
           scale_color_gradientn(colors = viridis_pal()(30)))
  
})

legend_plot.list(lin_plot.list,extension = '',name='CD4 Treg trajectories each plot individually')


### assign subsets

idents <- read_excel('Step2 Lymphoid cell/CD4/CD4T subcluster.xlsx')

obj$CD4.subset <- 'NA'

for(n in 1:nrow(idents)){
  id <- idents$SCT.res.2[n]
  obj$CD4.subset[which(obj$SCT_snn_res.2 == id)] = idents$ident[n]
}

DimPlot(obj,group.by = 'CD4.subset',raster = F,label = T) + theme(legend.position = 'none')

## plot signature score and cell markers for each cell type

avg_heatmap(obj,signatures.mm %>% names,group_by = 'CD4.subset' )+theme(axis.text.x=element_text(angle=90))



DotPlot(obj ,list(
    # Thf
    Thf = c('Bcl6','Cxcr5','Ebi3','Il21'),
    
    Th17 = c('Il17a','Il17f'),
    #Th1
    Cytotoxic = c('Ifng','Gzmb','Gzmk','Prf1','Nkg7'),
  
    # memory
    memory = c('Il7r','Cd69','Gpr183','Cd44'),
    
    Naive = c('Tcf7','Lef1','Sell'),
    
    Treg = c('Pdcd1','Foxp3','Il2ra'),
    
    `IFN response` = c('Stat1','Isg15','Ifit1')
    
  ),group.by = 'CD4.subset')&
    scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))
  




### cell subset percentage

######### percentage plot in CD4 subset
c<-table(obj@meta.data[,c('CD4.subset','sample_id')]) %>% prop.table(margin = 2) %>% as.data.frame()
meta <- read_xlsx("metadata.xlsx")
c$group <- NA
for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
}

c$Freq <- c$Freq *100
# ggplot(c,aes(x=CD4.subset,y=Freq* 100,colour=group))+
#   geom_boxplot(outliers = F)+
#   geom_point(position = position_jitterdodge(jitter.width=0.1) )+
#   theme_bw()+
#   ylab('Cell percent %')

pct_plot(c,'CD4.subset') + ylab('Cell percent %')

######### percentage plot in total tumor cells

c<-table(obj@meta.data[,c('CD4.subset','sample_id')]) %>% as.data.frame()

meta <- read_xlsx("metadata.xlsx")
c$group <- NA
c$total_number <- NA

for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
  c$total_number[i] <- meta$`Total cells after quality control`[which(meta$Sample==c$sample_id[i])]
}

c$Freq <- c$Freq/c$total_number *100
# ggplot(c,aes(x=CD4.subset,y=Freq* 100,colour=group))+
#   geom_boxplot(outliers = F)+
#   geom_point(position = position_jitterdodge(jitter.width=0.1) )+
#   theme_bw()+
#   ylab('Cell percent %')

pct_plot(c,'CD4.subset') + ylab('Cell percent %')


# ggplot(c,aes(x=CD4.subset,y=(Freq/total_number * 100),colour=group))+
#   geom_boxplot(outliers = F)+
#   geom_point(position = position_jitterdodge(jitter.width=0.1) )+
#   theme_bw()+ 
#   #geom_signif(comparisons = list(c('N8KO aPD1','WT aPD1')))+
#   ylab('Cell percent %')

 
write_rds(obj,'cd4.rds')




###### NK cell

obj <- read_rds('1.general.ident.assigned.rds')

obj <- obj%>%subset(celltype.general.1=='NK cell')

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

## re normalize and cluster
obj <- SCTransform(obj,vars.to.regress = c("mitoRatio"), clip.range = c(-5,5))

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) 

Idents(obj) <- 'SCT_snn_res.0.8'
DimPlot(obj,label = T,group.by = 'SCT_snn_res.0.8')

### remove containmination
FeaturePlot(obj,c('Cd3d','Cd3e','Krt8','Cd8a','Cd4','Ncr1','Cd14'))

obj <- obj%>%subset( SCT_snn_res.0.8 %in% c(10,11,5,12) == F)

### reclustering

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

## re normalize and cluster
obj <- SCTransform(obj,vars.to.regress = c("mitoRatio"), clip.range = c(-5,5))

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  RunUMAP(dims = 1:50,seed.use = 123,reduction='pca')%>%
  FindNeighbors( dims = 1:50, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) 

DimPlot(obj,label = T,group.by = 'SCT_snn_res.0.8')
## check containmination
FeaturePlot(obj,c('Cd3d','Cd3e','Krt8','Cd8a','Cd4','Ncr1','Cd14'))

## check batch effect
DimPlot(obj,label = F,group.by = 'SCT_snn_res.0.8',split.by = 'sample_id',ncol=5)

#write_rds(obj,'NK.rds')



##### resume process
obj <- read_rds('Step2 Lymphoid cell/NK/NK.rds')

FeaturePlot(obj,'Klrb1c')

Idents(obj) <- 'SCT_snn_res.0.8'
DimPlot(obj,label = T)

DotPlot(obj,c('Il7r','Rora','Gata3','Eomes','Ncr1','Krt83'))&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))


# determine ILC
# obj$nk.subset <- paste0('NK_',obj$SCT_snn_res.0.8)
# obj$nk.subset[which(obj$SCT_snn_res.0.8==3)] <- 'ILC'

# Idents(obj) <- 'nk.subset'
# DimPlot(obj,label = T)
# 
# DotPlot(obj,list(cytotoxic = c('Gzma','Gzmb','Gzmk','Ifng','Prf1'),
#                  exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox')),group.by = 'nk.subset')&
#   scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))
# 
# 
# obj <- obj %>% subset(nk.subset != 'ILC')
### calculate signature

 Markers <- read_xlsx('Step2 Lymphoid cell/NK/NK signature.xlsx',skip = 1,sheet = 'GSEA') 
 obj <- calculate_sig(obj,Markers,h2m_conversion=F,split.by = ',')

Markers <- read_xlsx('Step2 Lymphoid cell/NK/NK signature.xlsx',skip = 0,sheet = 'human_GSEA') 
obj <- calculate_sig(obj,Markers,h2m_conversion=T,split.by = ',')

Markers <- read_xlsx('Step2 Lymphoid cell/NK/NK signature.xlsx',skip = 1,sheet = 'human') 
obj <- calculate_sig(obj,Markers,h2m_conversion=T,split.by = NULL)

# calculate signature

signatures.mm <- obj[[2]]
obj <- obj[[1]]

###
DotPlot(obj,list(cytotoxic = c('Gzma','Gzmb','Gzmk','Ifng','Prf1'),
                 exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox','Cd96'),
                 migration = c('Ccl3', 'Ccl4', 'Ccl5','Cxcr4','Itgam'),
                 Intratumoral = c('Itga1')),group.by = 'SCT_snn_res.0.8')&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

DotPlot(obj,signatures.mm %>% names,group.by = 'SCT_snn_res.0.8')&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))
####


avg_heatmap(obj,names(signatures.mm),group_by = 'SCT_snn_res.0.8')


avg_heatmap(obj,names(signatures.mm),group_by = 'sample_id')

avg_heatmap(obj,names(signatures.mm),group_by = 'group',x.order = c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))

FeaturePlot(obj,names(signatures.mm))&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))


## check NK cell number
table(obj$sample_id)

## downsample based on the lowest number of cell
cell_downsample <- data.frame(cell.name = obj@meta.data %>% rownames(),sample_id = obj@meta.data$sample_id) %>%
                   group_by(sample_id) %>% 
                   sample_n(min(table(obj$sample_id)))
sub_obj <- obj[,cell_downsample$cell.name]
table(sub_obj$sample_id)


avg_heatmap(sub_obj,names(signatures.mm),group_by = 'group',x.order = c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))


VlnPlot(sub_obj,names(signatures.mm),group.by = 'group')




###
# 
# sub_obj <- obj

sub_obj$group <- factor(sub_obj$group,levels = c('WT_IgG','WT_aPD1','N8KO_IgG','N8KO_aPD1'))

DotPlot(sub_obj,list(cytotoxic = c('Gzma','Gzmb','Gzmc','Ifng','Prf1'),
                     exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox','Cd96'),
                     migration = c('Ccl3', 'Ccl4', 'Ccl5','Cxcr4'),
                     Intratumoral = c('Itga1','Itga2','Itgam'),
                     TFs = c('Eomes','Stat1','Bcl2','Isg15','Tcf7'),
                     proliferation = c('Mki67','Top2a')),group.by = 'group')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))


####


DotPlot(sub_obj,list(cytotoxic = c('Gzma','Gzmb','Gzmc','Ifng','Prf1'),
                     exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox','Cd96'),
                     migration = c('Ccl3', 'Ccl4', 'Ccl5','Cxcr4'),
                     Intratumoral = c('Itga1','Itga2','Itgam'),
                     TFs = c('Eomes','Stat1','Bcl2','Isg15','Tcf7'),
                     proliferation = c('Mki67','Top2a')),group.by = 'sample_id')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  geom_hline(yintercept = c(4.5,8.5,11.5),linetype ='dashed')

###

DotPlot(sub_obj,list(cytotoxic = c('Gzma','Gzmb','Gzmk','Ifng','Prf1'),
                     exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox','Cd96'),
                     migration = c('Ccl3', 'Ccl4', 'Ccl5','Cxcr4'),
                     Intratumoral = c('Itga1','Itga2','Itgam'),
                     TFs = c('Eomes','Rora','Il7r','Tbx21','Gata3')),group.by = 'SCT_snn_res.0.8')&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

DotPlot(sub_obj,list(cytotoxic = c('Gzma','Gzmb','Gzmk','Ifng','Prf1'),
                     exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox','Cd96'),
                     migration = c('Ccl3', 'Ccl4', 'Ccl5','Cxcr4'),
                     Intratumoral = c('Itga1','Itga2','Itgam'),
                     TFs = c('Eomes','Rora','Il7r','Tbx21','Gata3')),group.by = 'group')&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

### CD49b- tissue resident NK and ILC  possess memory potential
### CD49a tumor infiltration level, quick express after entering TME
### CD49b circulating NK cells
### Itgam mature cNK


## NK markers dotplot

Markers <- read_xlsx('Step2 Lymphoid cell/NK/NK signature.xlsx',skip = 0,sheet = 'mouse') 

Markers <- sapply(Markers %>% as.list(),function(x) x[is.na(x)==F])

DotPlot(obj,Markers)&
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

#NK cluster density - downsampled by sample
dim_dt <- sub_obj@reductions$umap@cell.embeddings %>%as.data.frame()
dim_dt$group <- sub_obj$group

ggplot(dim_dt,aes(x=umap_1,y=umap_2))+
  facet_wrap(~group)+
  geom_density_2d_filled()

#NK cluster density - downsampled by group
dim_dt <- obj@reductions$umap@cell.embeddings %>%as.data.frame()
dim_dt$group <- obj$group

dim_dt <- dim_dt %>% group_by(group) %>%sample_n(min(table(dim_dt$group)))
ggplot(dim_dt,aes(x=umap_1,y=umap_2))+
  facet_wrap(~group)+
  geom_point(alpha = 0.7,size = 0.3)+
  geom_density_2d_filled(alpha = 0.8)+
  theme_base ()




##DEG
Idents(obj) <- 'group'
mks <- FindMarkers(obj,ident.1 = 'N8KO_aPD1',ident.2 = 'WT_aPD1')
DEG2RNK(mks,log2fc_hold = 0.15,p_hold = 0.05,name ='N8KO_aPD1 vs WT_aPD1')

mks <- FindMarkers(obj,ident.1 = 'N8KO_IgG',ident.2 = 'WT_IgG')
DEG2RNK(mks,log2fc_hold = 0.15,p_hold = 0.05,name ='N8KO_IgG vs WT_IgG')

Idents(obj) <- 'SCT_snn_res.0.8'
mks <- FindMarkers(obj,ident.1 = '3',min.pct = 0.15)




FeaturePlot(obj,c('Il2ra','Il7r','Cxcr6'))
FeaturePlot(obj,c('Itga2'))
FeaturePlot(obj,c('Itga1'))
FeaturePlot(obj,c('Itgam'))


