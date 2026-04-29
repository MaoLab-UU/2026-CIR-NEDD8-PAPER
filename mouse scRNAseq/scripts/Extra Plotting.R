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
library(ggh4x)
library(cowplot)
library(ggsignif)
library(ggrepel)
library(ggnewscale)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 32000*1024^2)
lapply(list.files("C:/bioinfo/R tools",full.names = T),source) %>% invisible()





########  general cell idents UMAP

obj <-  read_rds('cleaned.cells.rds')


###########





axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm"),
)


p<- DimPlot(obj,group.by = 'celltype.general.1',
        label = T,
        label.size = 4,
        repel=TRUE,
        raster = F)+
  xlab('UMAP1')+
  ylab('UMAP2')+
  theme(
  aspect.ratio = 1,
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line(arrow = arrow(angle = 45,length = unit(0.25,'cm'))),
  axis.title = element_text(hjust = 0.05),
  legend.key.height = unit(0.8, 'cm'),   
  legend.text = element_text(size=16) )+ 
  
  
  guides( x = axis, y = axis) +
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  ggtitle(label = NULL)


ggsave('all cell UMAP.png',p,width = 6)


### general cell markers

obj$celltype.general.1 <- factor(obj$celltype.general.1,levels = c("T cell",'NK cell',
                                                                   'Macrophage','DCs',
                                                                   'Epithelial cell',
                                                                   'Fibroblast',
                                                                   'Endothelial cell',
                                                                   'Mast cell'
                                                                   ))

p1 <- DotPlot(obj,group.by = 'celltype.general.1',features  = list(
  Im = c('Ptprc'),
  `T cell` = c('Cd3d','Cd3e','Cd4','Cd8a'),
  `NK cell` = c('Klrd1','Ncr1','Klrb1c'),
  `Macrophage` = c('Adgre1','C1qa','Lyz2'),
  DCs = c('Xcr1','Batf3','Ccr7','Itgax'),
  Epi = c('Krt8','Itga6'),
  Fibro = c('Col1a2','Col1a1'),
  Endo = c('Pecam1','Cdh5'),
  Mast = c('Tpsab1','Mcpt4')
  
))+scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,vjust = 0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggsave('general marker 1.png',p1,width = 8,height = 4)




### cell percentage, general cell percentage in tumor

c<-table(obj@meta.data[,c('celltype.general.1','sample_id')]) %>% prop.table(margin = 2) %>% as.data.frame()
meta <- read_xlsx("metadata.xlsx")
c$group <- NA
for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
}

c$Freq <- c$Freq *100
c$group <- factor(c$group,levels = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'))




p <- pct_plot(c,'celltype.general.1',compare = list(c('WT αPD1','N8KO αPD1'),
                                                   c('WT IgG','N8KO IgG'))) + 
  ylab('Cell percentage in tumor %')+
  scale_color_manual(values = c('#8ecbe8','#928ee8','#e88ed1','#e8a08e'))+
  theme(
    aspect.ratio = 0.3,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    legend.position = 'top',
    axis.text.x  = element_text(size = 12))+
  xlab(NULL)



ggsave('all cell pct.png',p,width = 14,height = 4)


############# Macrophage markers
obj <- read_rds('./Step3 Myeloid cell/Macrophage/Macro.rds')




macro  <- list(
  
  `Recruit` = c('Cxcl9','Cxcl10'), # immune attraction
  
  `Lipid` = c('Apoe','Lipa','Fabp5','Lpl','Trem2'), # lipid metabolism
  
  `MHC-II` = c('Cd74','H2-Ab1','H2-Aa','H2-Eb1','H2-DMa','H2-DMb1'), #,'H2-DMb2'
#  `MHC-I` = c('B2m','H2-Q4','H2-Q6','H2-Q7','H2-T24','H2-M3'),#,'H2-T23'

  IFN_res = c('Ifnar1','Stat1','Irf1','Isg15','Isg20'),
  #  IFN_res = c('Ifnar1','Stat1','Irf1','Irf7', 'Ifi47','Ifit3','Ifit2','Isg15','Isg20'),
  
  `co-stim` = c('Cd86','Cd80'),
  
  IC = c('Pdcd1lg2','Cd274'))


p1 <- 
  DotPlot(obj,macro,group.by = 'group')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,vjust = 0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggsave('Macrophage marker 1.png',p1,width = 8,height = 2.2)






### heat map marcrophage



signatures <- c( 
"IFN-TAMs" ,                                                   
"Inflam-TAMs",                                                  
"Reg-TAMs",                                                   
"REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",             
"REACTOME_FATTY_ACID_METABOLISM",                             
"GOBP_MACROPHAGE_ACTIVATION",                                  
"REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",                   
"REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION")



avg_heatmap(obj,signatures,group_by = 'group',wrap_y = 24,
            x.order = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'),
            y.order = signatures
            )+
        theme_light() +
        theme(axis.text.y = element_text(size=9),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size = 10,angle = 90))+
        labs(fill = "Scaled Score")

ggsave('macrophage signatures.png',width = 5.4,height = 5)




############ CD4

obj <- read_rds('Step2 Lymphoid cell/CD4/cd4.rds')

obj$CD4.subset <- str_replace_all(obj$CD4.subset,' - ','\n')


### functional signature plots in general cd4



signatures <- c( 
  'Exhaustion',	
  'IFN response',
  'Treg signature',
  'Costimulatory molecules'

)

### signature heatmap

avg_heatmap(obj,signatures,group_by = 'group',wrap_y = 26,
            x.order = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'),
            y.order = signatures
)+
  theme_light() +
  theme(axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5))+
  labs(fill = "Scaled Score")

ggsave('CD4 signatures.png',width = 5,height = 4)


### cd4 GSEA + volcano
Idents(obj) <- 'group'
DEG.IgG <- FindMarkers(obj,ident.1 = 'N8KO IgG',ident.2 = 'WT IgG')
DEG2RNK(DEG.IgG,log2fc_hold = 0.1,p_hold = 0.1,name = 'CD4_DEG_N8KO_IgG-vs-WT_IgG')
write.csv(DEG.IgG,'CD4.DEG.IgG.csv')


DEG.aPD1 <- FindMarkers(obj,ident.1 = 'N8KO αPD1',ident.2 = 'WT αPD1')
DEG2RNK(DEG.aPD1,log2fc_hold = 0.1,p_hold = 0.1,name = 'CD4_DEG_N8KO_aPD1-vs-WT_aPD1')
write.csv(DEG.aPD1,'CD4.DEG.aPD1.csv')

volcano('C:/bioinfo/irineos_mice_scRNA_final_version/Step2 Lymphoid cell/CD4/DEG',
        width = 8,
        height = 4,
        tops = 30,log2_fc_hold = 0.35,p_hold = 0.01)


## before GSEA check if you have run easyBuild()
## Use the full path

GSEA_batch.rnk(
  GSEA_installed_path = "C:/bioinfo/irineos_mice_scRNA_final_version/GSEA_4.3.2",
  DEG_path='C:/bioinfo/irineos_mice_scRNA_final_version/Step2 Lymphoid cell/CD4/DEG',
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


### CD4 subset marker expression


p <- DotPlot(obj ,c(
  'Bcl6','Cxcr5','Ebi3','Il21',
  'Ifng','Gzmb','Gzmk','Prf1','Nkg7',
  'Il7r','Cd69','Gpr183','Cd44',
  'Tcf7','Lef1','Sell','Ccr7',
  'Pdcd1','Foxp3','Il2ra',
  'Stat1','Isg15','Ifit1')
  
  #### 'Il10'    'Tgfb'
  
,group.by = 'CD4.subset')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,vjust = 0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggsave('CD4 marker.png',p,width = 8,height = 2.5)





### Treg featureplot


p.l <- lapply(c('Cd3d','Cd3e','Cd4','Il2ra','Pdcd1','Foxp3'), 
            function(x) {FeaturePlot(obj,x)+
              scale_color_continuous(
                low ='#EEEEEF',
                high ='#E0191B',
                name = 'Expression',
                breaks = c(min(obj@assays$SCT@data[x,]),max(obj@assays$SCT@data[x,])),
                labels = c("min", "max")  # Set the legend labels to "min" and "max"
              )+
              theme(plot.margin = margin(0,0,0,0, "cm"),
                    axis.line=element_blank(),axis.text.x=element_blank(),
                    axis.text.y=element_blank(),axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    legend.title = element_text(size=10,vjust = 2))}
            )



p <- ggarrange(plotlist=p.l , widths = c(1,1), common.legend = T,legend="right",
          hjust = 0,
          vjust = 0,font.label = list(size = 10, color = "black", face = "bold", family = NULL) )

ggsave('Cd4 marker features.png',p,width = 8,height = 6)





### CD4 Treg pct per sample

# obj$subset.show <- 'CD4 other'
# obj$subset.show[grep('Treg',obj$CD4.subset)] <- 'CD4 Treg'


p<-DimPlot(obj,group.by = 'CD4.subset')+
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size=12),
    legend.key.spacing.y = unit(0.2,'cm'))+ 
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  ggtitle(label = NULL)

ggsave('Cd4 Treg dimplot.png',p,width = 5,height = 5)


### subset pct plot


c<-table(obj@meta.data[,c('CD4.subset','sample_id')]) %>% prop.table(margin = 2) %>% as.data.frame()

meta <- read_xlsx("metadata.xlsx")
c$group <- NA
for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
}

c$Freq <- c$Freq *100
c$group <- factor(c$group,levels = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'))



p <- pct_plot(c,group.by = 'CD4.subset',compare = list(c('WT αPD1','N8KO αPD1'),
                                                         c('WT IgG','N8KO IgG')))+

  ylab('CD4 subset percentage in CD4 %')+
    
  scale_color_manual(values = c('#8ecbe8','#928ee8','#e88ed1','#e8a08e'))+
  
  theme(
    aspect.ratio = 0.45,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = 'top',
    axis.text.x  = element_text(size = 10),
    legend.margin = margin(0,0,0,0))+
  xlab(NULL)



ggsave('CD4 cell pct in CD4.png',p,width = 10,height = 4)



### CD8


obj <- read_rds('Step2 Lymphoid cell/CD8/cd8.rds')

obj$CD8.subset <- factor(obj$CD8.subset,levels = c(
  "Naive CD8",
  "Early activated CD8",
  "Effector CD8",
  "Exhausted CD8",
  "Central Memory CD8",
  "Proliferating CD8"
))

### signature plot by group




### marker dotplot by cell type
mks <- c('Ccr7','Tcf7','Il7r','Sell','Lef1',
            'Gzmk','Gzmb','Prf1','Cd44','Cd69','Ifng','Eomes','Il2ra','Il2rb',
            'Klrd1','Klrb1','Prdm1','Cxcr3','Cxcr4',
            'Itga1','Itgae',
            'Pdcd1','Lag3','Havcr2','Ctla4',
            'Stat1','Isg15','Ifit1',
            'Mki67','Tuba1b','Top2a')
            

p <- DotPlot(obj,mks,group.by = 'CD8.subset')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,vjust = 0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggsave('CD8 marker by subset.png',p,width = 10,height = 2.5)


### CD8 dimplot

DimPlot(obj,group.by = 'CD8.subset')+
  theme(
  aspect.ratio = 1,
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_blank(),
  axis.title = element_blank(),
  legend.text = element_text(size=12),
  legend.key.spacing.y = unit(0.2,'cm'))+ 
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  ggtitle(label = NULL)

ggsave('Cd8 dimplot.png',width = 5,height = 5)

### subset pct plot

c<-table(obj@meta.data[,c('CD8.subset','sample_id')]) %>% prop.table(margin = 2) %>% as.data.frame()
meta <- read_xlsx("metadata.xlsx")

c$group <- NA
for(i in 1:nrow(c)){
  c$group[i] <- paste(meta$Genotype[which(meta$Sample==c$sample_id[i])],meta$Treatment[which(meta$Sample==c$sample_id[i])])
}

c$Freq <- c$Freq *100
c$group <- factor(c$group,levels = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'))



p <- pct_plot(c,group.by = 'CD8.subset',compare = list(c('WT αPD1','N8KO αPD1'),
                                                       c('WT IgG','N8KO IgG')))+
  
  ylab('CD8 subset percentage in CD8 %')+
  
  scale_color_manual(values = c('#8ecbe8','#928ee8','#e88ed1','#e8a08e'))+
  
  theme(
    aspect.ratio = 0.4,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = 'top',
    axis.text.x  = element_text(size = 10),
    legend.margin = margin(0,0,0,0))+
  xlab(NULL)



ggsave('CD8 cell pct in CD8.png',p,height = 4)




### CD8 DEG

Idents(obj) <- 'group'
DEG.IgG <- FindMarkers(obj,ident.1 = 'N8KO IgG',ident.2 = 'WT IgG')
DEG2RNK(DEG.IgG,log2fc_hold = 0.1,p_hold = 0.1,name = 'CD4_DEG_N8KO_IgG-vs-WT_IgG')
write.csv(DEG.IgG,'CD8.DEG.IgG.csv')


DEG.aPD1 <- FindMarkers(obj,ident.1 = 'N8KO αPD1',ident.2 = 'WT αPD1')
DEG2RNK(DEG.aPD1,log2fc_hold = 0.1,p_hold = 0.1,name = 'CD4_DEG_N8KO_aPD1-vs-WT_aPD1')
write.csv(DEG.aPD1,'CD8.DEG.aPD1.csv')

volcano('C:/bioinfo/irineos_mice_scRNA_final_version/Plotting/CD8/DEG',
        width = 8,
        height = 4,
        tops = 20,log2_fc_hold = 0.35,p_hold = 0.01)

### GSEA

GSEA_batch.rnk(
  GSEA_installed_path = "C:/bioinfo/irineos_mice_scRNA_final_version/GSEA_4.3.2",
  DEG_path='C:/bioinfo/irineos_mice_scRNA_final_version/Plotting/CD8/DEG',
  species = 'mouse',
  gene_sets =c(`hallmark gene sets`='mh.all.v2024.1.Mm.symbols.gmt',
               `reactome gene sets` = 'm2.cp.reactome.v2024.1.Mm.symbols.gmt',
               `GOBP gene sets`='m5.go.bp.v2024.1.Mm.symbols.gmt',
               `wiki pathways` = 'm2.cp.reactome.v2024.1.Mm.symbols.gmt'),
  symbol_chip='Mouse_Gene_Symbol_Remapping_MSigDB.v2024.1.Mm.chip',
  out_dir='auto',
  GSEA_plots_number=10,
  collapse='Remap_Only'
)


GSEA_bubble('C:/bioinfo/irineos_mice_scRNA_final_version/Plotting/CD8/DEG')


### marker dotplot by group
mks <- list(cytotoxicity = c('Gzmb','Gzmk','Prf1','Ifng','Tnf'),
         Exhaustion = c('Pdcd1','Lag3','Havcr2','Ctla4','Tox'),
         `IFN res` = c('Stat1','Isg15','Isg20','Irf4','Irf8'))
         #Cytokines = c('Csf1','Ccl3','Ccl4','Ccl5'))

p <- DotPlot(obj,mks,group.by = 'group',scale.min = 10)+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,vjust = 0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggsave('CD8 marker by group.png',p,width = 6,height = 2)



### CD8 signature heatmap


signatures <- c( 
'Fatty acid metabolism',
'Glycolysis',
'IFN Response',
'MAPK Signaling',
'Stress response',
'Cytotoxicity',
'Exhaustion'
)



p<-avg_heatmap(obj,signatures,group_by = 'group',wrap_y = 26,
            x.order = c('WT IgG','N8KO IgG','WT αPD1','N8KO αPD1'),
            y.order = signatures
)+
  theme_light() +
  theme(axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5))+
  labs(fill = "Scaled Score")

ggsave('CD8 signatures.png',p,width = 4.5,height = 4)




#### nk

obj <- read_rds("Step2 Lymphoid cell/NK/NK.rds")

mks <-  list(Cytotoxic = c('Gzma','Gzmb','Gzmc','Ifng','Prf1'),
           Exhaustion = c('Pdcd1','Ctla4','Havcr2','Lag3','Tox','Cd96','Tigit'),
           Chemotaxis = c('Ccl3', 'Ccl4', 'Ccl5','Cxcr4'),
           Integrins = c('Itga1','Itga2','Itgam'),
           TFs = c('Eomes','Stat1','Isg15','Tcf7'))
  
p<-DotPlot(obj,mks,group.by = 'group')+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B')) +
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90,vjust = 0.4),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=7))

ggsave('NK marker by group.png',p,width = 8.5,height = 2.2)


#### Epi GSEA

obj <- read_rds('Step4 tumor cell/Epi.rds')

DimPlot(obj)

FeaturePlot(obj,c('Ptprc','Cd14','Itgam','H2-Aa','Adgre1'))


Idents(obj) <- 'group'
DEG.IgG <- FindMarkers(obj,ident.1 = 'N8KO IgG',ident.2 = 'WT IgG',min.pct=0.01)
DEG2RNK(DEG.IgG,log2fc_hold = 0.25,p_hold = 0.01,name = 'Plotting/Epi/DEG/IgG/Epi_DEG_N8KO_IgG-vs-WT_IgG')
write.csv(DEG.IgG,'Plotting/Epi/DEG/IgG/Epi.DEG.IgG.csv')


DEG.aPD1 <- FindMarkers(obj,ident.1 = 'N8KO αPD1',ident.2 = 'WT αPD1',min.pct=0.01)
DEG2RNK(DEG.aPD1,log2fc_hold = 0.25,p_hold = 0.01,name = 'Plotting/Epi/DEG/aPD1/Epi_DEG_N8KO_aPD1-vs-WT_aPD1')
write.csv(DEG.aPD1,'Plotting/Epi/DEG/aPD1/Epi.DEG.aPD1.csv')

#### Run string enrichment


### plot terms


terms.aPD1.up <- c(
  'Regulation of cell migration',
  'Inflammatory response',
  'Response to cytokine',
  'Lipid metabolic process',
  'Regulation of localization',
  'Endocytosis',
  'Response to oxidative stress',
  'Response to interferon-beta',
  'Granulocyte chemotaxis',
  'Activation of innate immune response'
)
terms.aPD1.dn <- c(
  'Cell cycle process',
  'Regulation of cell cycle',
  'Cell division',
  'Regulation of cell cycle phase transition',
  'Cellular response to stress',
  'DNA metabolic process',
  'Mitotic nuclear division',
  'Negative regulation of cell cycle'
)



dt <- read_tsv('Plotting/Epi/DEG/aPD1/enrichment.Process.up.apd1.tsv')
dt.2 <- read_tsv('Plotting/Epi/DEG/aPD1/enrichment.Process.down.apd1.tsv')
dt.2$signal <- -dt.2$signal

dt$sig <- 'up'

dt.2$sig <- 'dn'

dt <- dt %>% filter(`term description` %in% terms.aPD1.up)
dt.2 <- dt.2 %>% filter(`term description` %in% terms.aPD1.dn)

dt.all <- rbind(dt,dt.2)




dt.all$`term description` <- factor(dt.all$`term description`,
                                    levels = dt.all$`term description`[order(dt.all$signal,dt.all$signal)])



p <- ggplot(dt.all,aes(x=signal,y=`term description`))+
  
  
  geom_col(aes(colour = `false discovery rate`,
               fill =  `false discovery rate`),
           data = subset(dt.all, sig=='dn'),
           width = 0.2,alpha=0.5)+
  
  geom_point(aes(colour = `false discovery rate`,
                 size=`observed gene count`),
             data = subset(dt.all, sig=='dn'))+ 
  
  scale_colour_gradientn('FDR.Down',trans = 'log10',colours=brewer.pal(9,'Blues')[2:8]%>%rev(),
                         guide = guide_colorbar(reverse = TRUE,order = 2))+
  scale_fill_gradientn('FDR.Down',trans = 'log10',colours=brewer.pal(9,'Blues')[2:8]%>%rev(),
                       guide = guide_colorbar(reverse = TRUE,order = 2))+
 
  ##### add second scale for up
  new_scale("color") +
  new_scale("fill") +

  
  geom_col(aes(colour = `false discovery rate`,
               fill =  `false discovery rate`),
           data = subset(dt.all, sig=='up'),
           width = 0.2,alpha=0.5)+
  
  geom_point(aes(colour = `false discovery rate`,
                 size=`observed gene count`),
             data = subset(dt.all, sig=='up'))+ 
  
  scale_colour_gradientn('FDR.Up',trans = 'log10',colours=brewer.pal(9,'Reds')[2:8]%>%rev(),
                         guide = guide_colorbar(reverse = TRUE,order = 1))+
  scale_fill_gradientn('FDR.Up',trans = 'log10',colours=brewer.pal(9,'Reds')[2:8]%>%rev(),
                       guide = guide_colorbar(reverse = TRUE,order = 1))+


  theme_bw() +
  theme(panel.grid.major.y = element_line(linetype = 'dashed'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=14))+
  guides(size = guide_legend(override.aes = list(color='grey'),order = 3))
  
  #guides(size = guide_legend(override.aes = list(color='grey',
  #fill=guide_colorbar(reverse = TRUE))))
   
  

ggsave('up in N8KO apd1 vs WT apd1.png',p,width = 9,height = 7)





###########



terms.igg.up <- c(
  'Response to stress',
  'Vasculature development',
  'Circulatory system development',
  'Regulation of innate immune response',
  'Angiogenesis',
  'Positive regulation of cytokine production',
  'Response to oxidative stress',
  'Lipid metabolic process'
)

terms.igg.dn <- c(
  'Cellular response to stress',
  'Mitotic cell cycle',
  'DNA replication',
  'Oxidative phosphorylation',
  'DNA repair'
)



dt <- read_tsv('Plotting/Epi/DEG/IgG/enrichment.Process.up.igg.tsv')

dt.2 <- read_tsv('Plotting/Epi/DEG/IgG/enrichment.Process.down.igg.tsv')


dt.2$signal <- -dt.2$signal

dt$sig <- 'up'

dt.2$sig <- 'dn'

dt <- dt %>% filter(`term description` %in% terms.igg.up)
dt.2 <- dt.2 %>% filter(`term description` %in% terms.igg.dn)

dt.all.2 <- rbind(dt,dt.2)


dt.all.2$`term description` <- factor(dt.all.2$`term description`,
                                    levels = dt.all.2$`term description`[order(dt.all.2$signal,dt.all.2$signal)])






p2 <- ggplot(dt.all.2,aes(x=signal,y=`term description`))+
  
  geom_col(aes(colour = `false discovery rate`,
               fill =  `false discovery rate`),
           data = subset(dt.all.2, sig=='dn'),
           width = 0.2,alpha=0.5)+
  
  geom_point(aes(colour = `false discovery rate`,
                 size=`observed gene count`),
             data = subset(dt.all.2, sig=='dn'))+ 
  
  scale_colour_gradientn('FDR.Down',trans = 'log10',colours=brewer.pal(9,'Blues')[2:8]%>%rev(),
                         guide = guide_colorbar(reverse = TRUE,order = 2))+
  scale_fill_gradientn('FDR.Down',trans = 'log10',colours=brewer.pal(9,'Blues')[2:8]%>%rev(),
                       guide = guide_colorbar(reverse = TRUE,order = 2))+
  
  ## add new scale for up
  new_scale("color") +
  new_scale("fill") +
  
  geom_col(aes(colour = `false discovery rate`,
               fill =  `false discovery rate`),
           data = subset(dt.all, sig=='up'),
           width = 0.2,alpha=0.5)+
  
  geom_point(aes(colour = `false discovery rate`,
                 size=`observed gene count`),
             data = subset(dt.all, sig=='up'))+ 
  
  scale_colour_gradientn('FDR.Up',trans = 'log10',colours=brewer.pal(9,'Reds')[2:8]%>%rev(),
                         guide = guide_colorbar(reverse = TRUE,order = 1))+
  scale_fill_gradientn('FDR.Up',trans = 'log10',colours=brewer.pal(9,'Reds')[2:8]%>%rev(),
                       guide = guide_colorbar(reverse = TRUE,order = 1))+
  

  theme_bw() +
  theme(panel.grid.major.y = element_line(linetype = 'dashed'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=14))+
  guides(size = guide_legend(override.aes = list(color='grey'),order = 3))

#guides(size = guide_legend(override.aes = list(color='grey',
#fill=guide_colorbar(reverse = TRUE))))



ggsave('up in N8KO igg vs WT igg.png',p2,width = 8,height = 5)




### inregrated plot




dt.all$cohort <- 'αPD1 Comparison'
dt.all.2$cohort <- 'IgG Comprison'

plot.dt <- rbind(dt.all,dt.all.2)




p.all <- ggplot(plot.dt,aes(x=signal,y=`term description`))+
  facet_wrap(~cohort,scales = 'free')+
  
  geom_col(aes(colour = `false discovery rate`,
               fill =  `false discovery rate`),
           data = subset(plot.dt, sig=='dn'),
           width = 0.2,alpha=0.5)+
  
  geom_point(aes(colour = `false discovery rate`,
                 size=`observed gene count`),
             data = subset(plot.dt, sig=='dn'))+ 
  
  scale_colour_gradientn('FDR.Down',trans = 'log10',colours=brewer.pal(9,'Blues')[2:8]%>%rev(),
                         guide = guide_colorbar(reverse = TRUE,order = 2))+
  scale_fill_gradientn('FDR.Down',trans = 'log10',colours=brewer.pal(9,'Blues')[2:8]%>%rev(),
                       guide = guide_colorbar(reverse = TRUE,order = 2))+

  new_scale("color") +
  new_scale("fill") +
  
  
  geom_col(aes(colour = `false discovery rate`,
               fill =  `false discovery rate`),
           data = subset(plot.dt, sig=='up'),
           width = 0.2,alpha=0.5)+
  
  geom_point(aes(colour = `false discovery rate`,
                 size=`observed gene count`),
             data = subset(plot.dt, sig=='up'))+ 
  
  scale_colour_gradientn('FDR.Up',trans = 'log10',colours=brewer.pal(9,'Reds')[2:8]%>%rev(),
                         guide = guide_colorbar(reverse = TRUE,order = 1))+
  scale_fill_gradientn('FDR.Up',trans = 'log10',colours=brewer.pal(9,'Reds')[2:8]%>%rev(),
                       guide = guide_colorbar(reverse = TRUE,order = 1))+
  
  theme_bw() +
  theme(panel.grid.major.y = element_line(linetype = 'dashed'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=14))+
  guides(size = guide_legend(override.aes = list(color='grey'),order = 3))


ggsave('up in N8KO vs WT integrated .png',p.all,width = 16,height = 7)



### Cell chat

cellchat <- read_rds('cc,rds')

cell.colors <- scPalette(cellchat@net[[1]]$count %>% ncol())
names(cell.colors) <- cellchat@net[[1]]$count %>% colnames()




diff.heatmap <- function(obj,compares,cell.colors,slot='count'){
  
  dt <- obj@net[[compares[2]]][[slot]] - obj@net[[compares[1]]][[slot]]
  
  dt.col <- data.frame(group = colnames(dt))
  rownames(dt.col) <- colnames(dt)
  dt.row <- data.frame(group = rownames(dt))
  rownames(dt.row) <- rownames(dt)
  
  
  col_annotation <- HeatmapAnnotation(df = dt.col, col = list(group = cell.colors), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  row_annotation <- HeatmapAnnotation(df = dt.row, col = list(group = cell.colors), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  p <- Heatmap(dt,cluster_rows = F,
          cluster_columns = F,
          row_names_side = "left",
          col = c('#2760DA','#EEEEEF','#E0191B'),
          row_title = "Sources cell",
          column_title = 'Target cell',
          column_title_side = 'bottom',
          bottom_annotation = col_annotation,
          column_names_gp = gpar(fontsize = 12,fontface="plain"),
          row_names_gp = gpar(fontsize = 12,fontface="plain"),
          left_annotation = row_annotation,
          column_title_gp = gpar(fontsize = 12,fontface="plain"),
          row_title_gp = gpar(fontsize =12,fontface="plain"),
          heatmap_legend_param = list(title =paste0('Diff ',slot,' of interactions\n',compares[[1]],' vs ',compares[[2]]), 
                                      title_gp = gpar(fontsize = 10), 
                                      border = NA, 
                                      legend_height = unit(20, "mm"), 
                                      labels_gp = gpar(fontsize = 8,fontface="plain"),
                                      grid_width = unit(4, "mm")))
  return(p)
}


#


diff.heatmap(cellchat,c('WT IgG','N8KO IgG'),cell.colors,'weight')+

diff.heatmap(cellchat,c('WT αPD1','N8KO αPD1'),cell.colors,'weight')+


diff.heatmap(cellchat,c('N8KO IgG','N8KO αPD1'),cell.colors,'weight')+

diff.heatmap(cellchat,c('WT IgG','WT αPD1'),cell.colors,'weight')







#> Do heatmap based on a merged object

res <- 1200
length <- 480*res/72

png('heatmap_diff_interactions_IgG.png',width = length,height = length*0.80,pointsize=1,res=res)
gg1
dev.off()

png('heatmap_diff_interactions_aPD1.png',width = length,height = length*0.75,pointsize=1,res=res)
gg2
dev.off()





######## diff bubble plot

### DEG analysis


#αPD1



result <- read_rds('cc.result.rds') #'cc.igg.rds'

cellchat <- mergeCellChat(result,add.names =names(result))


compare <- c('N8KO IgG','WT IgG')


cellchat <- subsetCellChat(cellchat,cells.use = rownames(cellchat@meta)[cellchat@meta$datasets %in% c(compare)]   )


cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = compare[[1]],
                                       features.name = paste0(compare[[1]],"new"), 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.01, 
                                       thresh.fc = 0.05,
                                       thresh.p = 0.01, 
                                       group.DE.combined = FALSE) 





net <- netMappingDEG(cellchat, features.name = paste0(compare[[1]],"new"), variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = compare[[1]])
net.down <- subsetCommunication(cellchat, net = net, datasets = compare[[2]])
net.all <- rbind(net.up,net.down)

net.all$receptor.pvalues[is.na(net.all$receptor.pct.2)] <- 1
net.all$receptor.logFC[is.na(net.all$receptor.logFC)] <- 0




net.all$Logfc <- net.all$ligand.logFC + net.all$receptor.logFC
net.all$trend <- ifelse(net.all$Logfc > 0,'Up','Down')
#net.all$pct <- rowMeans(net.all[,c('ligand.pct.1','ligand.pct.2')])/100 * rowMeans(net.all[,c('receptor.pct.1','receptor.pct.2')])/100
net.all$pval <- lapply(1:nrow(net.all),function(x) min(net.all$ligand.pvalues[x], 
                                                       net.all$receptor.pvalues[x])) %>% unlist()
## keep the significant 
net.all <- net.all %>% filter(abs(net.all$Logfc) >0.1 & 
                                0 < net.all$pval & 
                                net.all$pval< 0.01)


plot.dt <- net.all %>% filter(source %in% c('Macrophage') &
                             target %in% c('CD4 Treg cell','CD8 T cell','NK cell'))

plot.dt$cellpair <- paste0(plot.dt$source,' -> ',plot.dt$target)

plot.dt$inter.name <- paste0(plot.dt$ligand,' - ',plot.dt$receptor)

plot.dt$pathway_name <- factor(plot.dt$pathway_name,levels = unique(plot.dt$pathway_name))

plot.dt <- plot.dt[order(plot.dt$pathway_name),]



ggplot(plot.dt,aes(x=cellpair,y=inter.name,colour = Logfc,size = prob))+
  geom_point()+
  ggtitle(paste0('Positive - Up in ',compare[[1]]))+

  scale_color_gradientn(colors  = c('#2760DA','#EEEEEF','#E0191B'))+theme_bw()




pairLR.use <- net.all[,"interaction_name", drop = F]

netVisual_bubble_new(cellchat, 
                        sources.use = 5, 
                        pairLR.use = pairLR.use,
                        targets.use = c(1,3,4), 
                        comparison = c(2, 4),
                        thresh = 0.01,
                        angle.x = 90, 
                        remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))




net.all.sub <- net.all %>% filter(source %in% c('CD4 Treg cell'),
                                  target %in% c('Macrophage','CD4 T cell','CD8 T cell','NK cell'))
pairLR.use <- net.all.sub [,"interaction_name", drop = F]
netVisual_bubble_new(cellchat, 
                     sources.use = 2, 
                     pairLR.use = pairLR.use,
                     targets.use = c(5,1,3,6), 
                     comparison = c(2, 4),
                     thresh = 0.01,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#2760DA','#EEEEEF','#E0191B'))












cellchat <- mergeCellChat(result,add.names =names(result))


compare <- c('N8KO αPD1','WT αPD1')


cellchat <- subsetCellChat(cellchat,cells.use = rownames(cellchat@meta)[cellchat@meta$datasets %in% c(compare)]   )


cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = compare[[1]],
                                       features.name = paste0(compare[[1]],"new"), 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.01, 
                                       thresh.fc = 0.05,
                                       thresh.p = 0.01, 
                                       group.DE.combined = FALSE) 




net <- netMappingDEG(cellchat, features.name = paste0(compare[[1]],"new"), variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = compare[[1]])
net.down <- subsetCommunication(cellchat, net = net, datasets = compare[[2]])
net.all <- rbind(net.up,net.down)

net.all[is.na(net.all)] <- 0

net.all$Logfc <- net.all$ligand.logFC + net.all$receptor.logFC
net.all$trend <- ifelse(net.all$Logfc > 0,'Up','Down')
net.all$pct <- rowMeans(net.all[,c('ligand.pct.1','ligand.pct.2')])/100 * rowMeans(net.all[,c('receptor.pct.1','receptor.pct.2')])/100
net.all$pval <- lapply(1:nrow(net.all),function(x) max(net.all$ligand.pvalues[x], 
                                                       net.all$receptor.pvalues[x])) %>% unlist()
## keep the significant 
net.all <- net.all %>% filter(abs(net.all$Logfc) >0.1 & 
                                net.all$pval< 0.05)

net.all.sub <- net.all %>% filter(source %in% c('Macrophage'),
                              target %in% c('CD4 Treg cell','CD8 T cell','NK cell'))

pairLR.use <- net.all.sub [,"interaction_name", drop = F]
plot.dt <- net.all %>% filter(source %in% c('Macrophage') &
                                target %in% c('CD4 Treg cell','CD8 T cell','NK cell'))

plot.dt$cellpair <- paste0(plot.dt$source,' -> ',plot.dt$target)

plot.dt$inter.name <- paste0(plot.dt$ligand,' - ',plot.dt$receptor)

plot.dt$pathway_name <- factor(plot.dt$pathway_name,levels = unique(plot.dt$pathway_name))

plot.dt <- plot.dt[order(plot.dt$pathway_name),]



ggplot(plot.dt,aes(x=cellpair,y=inter.name,colour = Logfc,size = prob))+
  geom_point()+
  ggtitle(paste0('Positive - Up in ',compare[[1]]))+
  
  scale_color_gradientn(colors  = c('#2760DA','#EEEEEF','#E0191B'))+theme_bw()






netVisual_bubble_new(cellchat, 
                     sources.use = 5, 
                     pairLR.use = pairLR.use,
                     targets.use = c(2,3,6), 
                     comparison = c(1, 3),
                     thresh = 0.01,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#2760DA','#EEEEEF','#E0191B'))







net.all.sub <- net.all %>% filter(source %in% c('DCs'),
                                  target %in% c('CD4 T cell','CD8 T cell','CD4 Treg cell'))

pairLR.use <- net.all.sub [,"interaction_name", drop = F]


netVisual_bubble_new(cellchat, 
                     sources.use = 2, 
                     pairLR.use = pairLR.use,
                     targets.use = c(5,3,6), 
                     comparison = c(1, 3),
                     thresh = 0.05,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))












netVisual_bubble_new(cellchat, 
                     sources.use = 4, 
                     pairLR.use = pairLR.use,
                     targets.use = c(5,3,1), 
                     comparison = c(1, 3),
                     thresh = 0.05,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))








netVisual_bubble_new(cellchat, 
                     sources.use = 4,
                     targets.use = c(2,3,1), 
                     comparison = c(1, 3),
                     max.dataset = 3,
                     thresh = 0.05,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))





















netVisual_bubble_new(cellchat, 
                     sources.use = 2,
                     targets.use = c(5,3,1), 
                     comparison = c(1, 4),
                     max.dataset = 4,
                     thresh = 0.05,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))



netVisual_bubble_new(cellchat, 
                     sources.use = 4,
                     targets.use = c(2,3,1), 
                     comparison = c(2, 3),
                     max.dataset = 3,
                     thresh = 0.1,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))




netVisual_bubble_new(cellchat, 
                     sources.use = 1,
                     targets.use = c(3), 
                     comparison = c(2, 3),
                     max.dataset = 3,
                     thresh = 0.1,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))









netVisual_bubble_new(cellchat, 
                     sources.use = 4,
                     targets.use = c(1,3), 
                     comparison = c(2, 3),
                     max.dataset = 3,
                     thresh = 0.1,
                     angle.x = 90, 
                     remove.isolate = T,return.data = F)  +
  scale_color_gradientn(colors  = c('#EEEEEF','#E0191B'))
















