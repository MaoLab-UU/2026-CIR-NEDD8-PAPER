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
obj <- obj %>% subset(celltype.general.1 =='Epithelial cell')

### re cluster remove contaminants

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

obj <- SCTransform(obj, 
                   vars.to.regress = c("mitoRatio"),
                   clip.range = c(-5,5),
                   verbose = TRUE)

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:20, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')

FeaturePlot(obj,c('Ptprc','Cd3d','Cd3e','Klrb1c','Cd14','Itgax'))

DimPlot(obj,group.by = 'SCT_snn_res.0.8',label = T)

obj <- obj %>% subset(SCT_snn_res.0.8 %in% c(11,12,15,14)==F)

### RE- RE CLUSTER

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

obj <- SCTransform(obj, 
                   vars.to.regress = c("mitoRatio"),
                   clip.range = c(-5,5),
                   verbose = TRUE)

obj <- obj %>%
  RunPCA(assay = "SCT", npcs = 50)%>%
  FindNeighbors( dims = 1:20, reduction = "pca",nn.method = "rann")%>%
  FindClusters(resolution = c(0.8,2)) %>%
  RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')

FeaturePlot(obj,c('Ptprc','Cd3d','Cd3e','Klrb1c','Cd14','Itgax'))




#write_rds(obj,'Epi.rds')

########### Infer CNV

obj.nk <- read_rds("Step2 Lymphoid cell/NK/NK.rds") %>% colnames()

obj.macro <- read_rds("Step3 Myeloid cell/Macrophage/Macro.rds") %>% colnames()

obj.epi <- read_rds("Step4 tumor cell/Epi.rds") %>% colnames()

obj.cd8 <- read_rds("Step2 Lymphoid cell/CD8/cd8.rds") %>% colnames()

obj.cd4.treg <- read_rds("Step2 Lymphoid cell/CD4/cd4.rds") %>% 
  subset(CD4.subset %>% str_detect('Treg')) %>% colnames()

obj.cd4 <- read_rds("Step2 Lymphoid cell/CD4/cd4.rds") %>% 
  subset(CD4.subset %>% str_detect('Treg') == F) %>% colnames()

obj.dc <- read_rds ('Step3 Myeloid cell/DC/DC.rds') %>% colnames()











obj <- read_rds('1.general.ident.assigned.rds')


### Infer CNV

DefaultAssay(obj) <- 'RNA'
obj[['SCT']] <- NULL

immune <- c(obj.cd4,obj.cd8,obj.cd4.treg,obj.macro,obj.dc,obj.nk) %>% sample(1000)

Epi <- obj.epi %>% sample(1000)



sub_obj <- obj[,c(immune,Epi)]



sub_obj$celltype.general.1 %>% unique()

counts_matrix = GetAssayData(sub_obj,assay='RNA',layer = 'counts')

metadata <- data.frame(type=sub_obj$celltype.general.1)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=metadata,
                                    delim="\t",
                                    gene_order_file="Step4 tumor cell/inferCNV/mm10 gene position.txt",
                                    ref_group_names=c('NK cell','T cell',"DCs","Macrophage")) 

#rm(sub_obj)
rm(counts_matrix)
options(scipen = 100)
infercnv_obj = try(infercnv::run(infercnv_obj,
                                 cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir='Step4 tumor cell/inferCNV', 
                                 num_threads = 12,
                                 window_length=101,
                                 max_centered_threshold=3,
                                 cluster_by_groups=T,
                                 plot_steps=F,
                                 denoise=T,
                                 sd_amplifier=2,
                                 analysis_mode = "subclusters",
                                 HMM = T,
                                 leiden_resolution = 0.01))





#### process infercnv

cnv_objs_rds <- list.files('Step4 tumor cell/inferCNV',
                           pattern = '.rds',
                           full.names = T,recursive = T)
metadata <- data.frame(type=sub_obj$celltype.general.1)
CNV_score_hold <- 1.025
CNV_correlation_hold <- 0.5
CNV.dt <- process_infercnv(cnv_objs_rds) 




CNV.dt$celltype <- metadata[rownames(CNV.dt),]


write_rds(CNV.dt,'CNV.dt.rds')






#### plotting


###1``
CNV.dt <- read_rds('CNV.dt.rds')

p<- ggplot(CNV.dt,aes(x=CNV_correlation,y=CNV_score))+
  geom_point(aes(color=cell),alpha=1,size=0.5)+
  geom_hline(yintercept = CNV_score_hold,linetype='dashed',alpha=0.5)+
  geom_vline(xintercept = CNV_correlation_hold ,linetype='dashed',alpha=0.5)+
  scale_color_hue(labels=c('Observations','References'))+
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  ylim(c(0.9,1.1))

ggsave('cnv_plot1.png',p,width = 5,height = 4)

###2

# p <- ggplot(CNV.dt,aes(x=CNV_correlation,y=CNV_score))+
#   facet_wrap(~ celltype)+
#   geom_point(aes(color=celltype),alpha=1,size=0.5)+
#   geom_hline(yintercept = CNV_score_hold,linetype='dashed',alpha=0.5)+
#   geom_vline(xintercept = CNV_correlation_hold ,linetype='dashed',alpha=0.5)+
#   ylim(c(0.9,1.1))+
#   theme_bw()
# 
# ggsave('cnv_plot_by_celltype.png',p,width = 7,height = 5)



###3

normal_pct.dt <- data.frame()
for(s in unique(CNV.dt$celltype)){
  dt <- CNV.dt %>% filter(celltype == s)
  
  
  
  normal_pct <- which(dt$CNV_score <= CNV_score_hold & (dt$CNV_correlation <= CNV_correlation_hold)) %>% length() /
    nrow(dt) * 100
  
  normal_pct.dt <- rbind(normal_pct.dt,data.frame(nom_pct = round(normal_pct,3) ,celltype=s))
}


normal_pct.dt$nom_pct <-  normal_pct.dt$nom_pct%>% round(1)

ggplot(CNV.dt ,aes(x=CNV_correlation ,y=CNV_score))+
  facet_wrap(~ celltype)+
  geom_point(aes(color=malignancy),alpha=1,size=0.5)+
  geom_hline(yintercept = CNV_score_hold,linetype='dashed',alpha=0.5)+
  geom_vline(xintercept = CNV_correlation_hold ,linetype='dashed',alpha=0.5)+
  ggtitle('')+
  ylim(c(0.9,1.1))+
  geom_text(data = normal_pct.dt ,
            aes(x = 0.1, y = 0.95, label = paste0('Normal pct: ',nom_pct,'%')),size = 3,alpha=0.75)+
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(panel.spacing = unit(0.8, "lines"))+
  scale_color_manual(values = c('#FF966D','#946DFF'))



ggsave('cnv_normal_pct_by_celltype.png',width = 8,height = 5)




