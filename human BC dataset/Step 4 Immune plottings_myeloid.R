library(tidyverse)
library(ggplot2)
library(Seurat)
library(readxl)
library(readr)
library(stringr)
library(scales)
library(splitstackshape)
options(future.globals.maxSize= 32000*1024^2)


mye <- read_rds('Myeloid.rds')
mye <- mye %>% subset(expansion != 'n/a')
DimPlot(mye,label = T)+
DotPlot(mye,features=c('CSF1R','ITGAM','ITGAX','CD14','CD80','XCR1','BATF3','CLEC9A','S100A8','S100A9','CCR7','HLA-DRB1'))


macro <- mye %>% subset(seurat_clusters %in% c(1,5,8,9,17,19,20,21,24,27,28,29) == F)

DimPlot(macro)

dt <- read_tsv('NEDD8_phenotype_label.txt')

macro$N8_status <- lapply(macro$patient_id,function(x) dt$N8_status[which(dt$patient_id == x)]) %>% unlist()

sig <- read_xlsx('macrophage markers.xlsx') %>% as.data.frame()

score <- calculate_sig(macro,sig,h2m_conversion = F,split.by = ',')
macro <- score[[1]]
dt <- score[[2]]


macro$group <- paste(macro$timepoint,macro$N8_status)

## Fig 1 E
avg_heatmap(macro,names(dt),group_by = 'group',wrap_y=26,y.order = names(dt),
            x.order = c('Pre N8_high','Pre N8_low','On N8_high','On N8_low') )+
  theme_light() +
  theme(axis.text.y = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10,angle = 90))+
  labs(fill = "Scaled Score")

### Fig 1 F
### 
### 

signatures <- list(`M1 like` = c('IL1B','TNF','IL6','CXCL9','CXCL10','CXCL11','CCL5',
                                 'NFKBIA','STAT1','IRF1','IRF5','SOCS3'),
                   `M2 like` = c('CD163','MRC1','MSR1','MARCO','IL10','TGFB1','ARG1',
                                 'FN1','VCAN','SPP1'
                   ),
                   `MHC-II` = c('HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74'),
                   `MHC-I` = c('HLA-A','HLA-B','HLA-C','HLA-E','B2M','TAP1','TAP2'))

DotPlot(macro,group.by = 'group',features = signatures)+
  scale_color_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.4,'cm'))
