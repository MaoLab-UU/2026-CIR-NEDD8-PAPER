library(tidyverse)
library(ggplot2)
library(Seurat)
library(readxl)
library(readr)
library(stringr)
library(scales)

######## Tumor bulk DEG analysis
obj <- read_rds('all_cells_all_type_cancer.rds')
Tumor <- obj %>% subset(cellType == 'Cancer_cell')

### Tumor DEG Pre_E vs Pre_NE single cell level 
#Fig 1A 

Idents(Tumor) <- 'expansion'
mks <- FindMarkers(Tumor %>% subset(timepoint == 'Pre'),ident.1 = 'E',ident.2 = 'NE')
write.csv2(mks,'./DEG/single cell/Tumor pre deg.csv')


#mks<-read.csv2('./DEG/single cell/Tumor pre deg.csv')
#rownames(mks) <- mks$X

mks <- mks %>% filter(pct.1 >0.05,pct.2>0.05)
mks$pct.difference <- mks$pct.1/mks$pct.2

highlight.keys <- read_tsv('./DEG/single cell/highlight_key.txt.txt',col_names = F) %>% unlist()

ggplot(mks,aes(x=pct.difference,y=avg_log2FC))+
  geom_point(color='grey',alpha=0.5,size=1)+
  geom_point(data = mks[highlight.keys,],size=1,
             color=ifelse(mks[highlight.keys,]$avg_log2FC <=0 ,'#5dade2','#e25d5d'))+
  
  ggrepel::geom_text_repel(data = mks[highlight.keys,],
                           label = highlight.keys,
                           color ='black',
                           max.overlaps = length(highlight.keys)+1,
                           segment.size = 0.1,
                           segment.linetype = "dashed",
                           nudge_x = runif(n=length(highlight.keys), min=-0.2, max=0.8),
                           nudge_y = runif(n=length(highlight.keys), min=-0.8, max=0.2),
                           show.legend = F,
                           size=2.5,
                           alpha=0.75)+
  xlab('pct of cell express in E / pct of cell express in NE')+
  theme_bw()+
  geom_hline(yintercept = 0,alpha=0.2,linetype='dashed')+
  geom_vline(xintercept = 1,alpha=0.2,linetype='dashed')


## FIG S1A Tumor pseudo bulk RNAseq DEG and GSEA


Tumor$group <- paste(Tumor$timepoint,Tumor$expansion)
Idents(Tumor) <- 'group'

library(DESeq2)

exp <- AggregateExpression(Tumor %>% subset(group %in% c('Pre E','Pre NE')),
                           assays = 'RNA',group.by = 'patient_id',
                           normalization.method=NULL)[['RNA']] %>% as.data.frame()

meta <- unique(subset(Tumor,
                      group %in% c('Pre E','Pre NE'))@meta.data[,c('group','patient_id')])
rownames(meta) <- meta$patient_id
colnames(exp) <- colnames(exp) %>% str_replace('-','_')
meta <- meta[colnames(exp),]

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = meta,
                              design = ~ group)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$group <- factor(dds$group, levels = c("Pre NE","Pre E"))

dds <- DESeq(dds)
res <- results(dds)
result <- res %>% as.data.frame()

write.csv2(result,'./DEG/bulk/DESeq2_tumor.csv')

result['NEDD8',]

result <- result %>% filter(is.na(result$pvalue) == F)
result$label <- NA
result['NEDD8','label'] <- 'NEDD8'

#result <- read.csv2('./DEG/bulk/DESeq2_tumor.csv',row.names = 1)
highlights <- read_tsv('./DEG/bulk/highlight_key.txt.txt',col_names = F) %>% unlist()

png('TumorDEG_-VOLCANO.png',width = 5,height = 4,units = 'in',res = 800)
FindMks_Volcano(result,
                p_adj.hold = 0.05,
                avg_lfc.hold = 0.25,
                top_n_plot = 0,
                gene.highlight = highlights,
                log2fc = 'log2FoldChange',
                p_val = 'pvalue')

dev.off()


### GSEA following DESeq2 
## Fig S1D

result$gene <- rownames(result)
write_tsv(result[order(result$log2FoldChange,decreasing = T),c('gene','log2FoldChange')] ,
          './DEG/bulk/tumor pre E_VS_NE.rnk',col_names = F)

GSEA_batch.rnk(
  GSEA_installed_path = "../R tools/GSEA_4.3.2",
  DEG_path='DEG/bulk',
  species = 'human',
  gene_sets =c(`hallmark gene sets`='h.all.v2025.1.Hs.symbols.gmt'),
  symbol_chip='Human_Gene_Symbol_with_Remapping_MSigDB.v2025.1.Hs.chip',
  out_dir='GSEA',
  GSEA_plots_number=30,
  collapse='Collapse'
)

GSEA_bubble_3('GSEA')

### TUMOR nedd8 expression by group 
### Fig 1B

exp <- get_average_expression(obj=Tumor,'NEDD8',assay = 'SCT',
                              group.by = c('timepoint','expansion'))

exp$group.by <- factor(exp$group.by,levels = c('Pre_NE','Pre_E','On_NE','On_E'))

png('NEDD8 expression based on N8_LOW_HIGH in Pre samples.png',width = 4.5,height = 4,units = 'in',res=800)

ggplot(exp,aes(x=group.by,y=value,colour = group.by))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = group.by))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45,hjust=-0.1),
        axis.text.x.bottom = element_blank(),
        axis.ticks.length.x = unit(0,'cm'))+
  scale_shape_manual(values = c(16, 17, 16, 17))+
  facet_wrap(~gene,ncol = 5,scales = 'free_y')+
  geom_signif(comparisons = list(c('Pre_NE','Pre_E'),c('On_NE','On_E')),test = 'wilcox.test')
  ylab('Normalized expression')

dev.off()


### paired tumor NEDD8 expression before after treatment
### #Fig S1C

exp <- get_average_expression(Tumor,'NEDD8',assay = 'SCT',
                              group.by = c('timepoint','expansion'),
                              remove_outliers_by_group = F)
pj <- position_jitter(0)

exp$group.by <- factor(exp$group.by,levels = c('Pre_NE','On_NE','Pre_E','On_E'))

png('NEDD8 exp in tumor paired expression boxplot.png',width = 6,height = 4.5,units = 'in',res = 800)

ggplot(exp,aes(x=group.by,y=value,colour = group.by,group=patient_id))+
  geom_boxplot(aes(group=NULL),outlier.shape = NA)+
  geom_point(aes(shape = group.by),position = pj)+
  geom_line(position = pj)+
  theme_bw()+
  geom_signif(comparison = list(c("Pre_E", "On_E"),c('Pre_NE','On_NE')), 
              test = "wilcox.test", test.args = list(paired = TRUE))+
  theme(axis.text.x = element_text(angle = -45,hjust=-0.1),
        axis.text.x.bottom = element_blank(),
        axis.ticks.length.x = unit(0,'cm'))+
  scale_color_manual(values = c('#b20f0f','#b20f0f','#0F99B2','#0F99B2'))+
  scale_shape_manual(values = c(16, 17, 16, 17))+
  facet_wrap(~gene,ncol = 5,scales = 'free_y')+
  ylab('Normalized expression')

dev.off()

####### define n8_low n8_high

dt <- AverageExpression(Tumor %>% subset(timepoint == 'Pre'),features = 'NEDD8',
                        group.by = 'patient_id',assays = 'SCT')  %>% as.data.frame() %>% 
  t() %>% as.data.frame()

colnames(dt) <- 'NEDD8'

dt$N8_status <- ifelse(dt$NEDD8 > mean(dt$NEDD8),'N8_high','N8_low' )

dt$patient_id <- rownames(dt)%>% str_remove_all('SCT.') %>% str_replace_all('\\.','_')

write_tsv(dt,'NEDD8_phenotype_label.txt')

Tumor$N8_status <- lapply(Tumor$patient_id,function(x) dt$N8_status[which(dt$patient_id == x)]) %>% unlist()


####### Tumor DEG + GSEA by NEDD8 expression
####### Fig S1D
Tumor$group <- paste(Tumor$timepoint,Tumor$N8_status)
Idents(Tumor) <- 'group'

library(DESeq2)

exp <- AggregateExpression(Tumor %>% subset(group %in% c('Pre N8_low','Pre N8_high')),
                           assays = 'RNA',group.by = 'patient_id',
                           normalization.method=NULL)[['RNA']] %>% as.data.frame()

meta <- unique(subset(Tumor,
                      group %in% c('Pre N8_low','Pre N8_high'))@meta.data[,c('group','patient_id')])
rownames(meta) <- meta$patient_id
colnames(exp) <- colnames(exp) %>% str_replace('-','_')
meta <- meta[colnames(exp),]

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = meta,
                              design = ~ group)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$group <- factor(dds$group, levels = c("Pre N8_high","Pre N8_low"))

dds <- DESeq(dds)
res <- results(dds)
result <- res %>% as.data.frame()

write.csv2(result,'./DEG/bulk/DESeq2_tumor.csv')

result['NEDD8',]

result <- result %>% filter(is.na(result$pvalue) == F)
result$label <- NA
result['NEDD8','label'] <- 'NEDD8'

#result <- read.csv2('./DEG/bulk/DESeq2_tumor.csv',row.names = 1)


FindMks_Volcano(result,
                p_adj.hold = 0.05,
                avg_lfc.hold = 0.25,
                top_n_plot = 0,
                gene.highlight = c(),
                log2fc = 'log2FoldChange',
                p_val = 'pvalue')


### GSEA following DESeq2 
## Fig S1D

result$gene <- rownames(result)
write_tsv(result[order(result$log2FoldChange,decreasing = T),c('gene','log2FoldChange')] ,
          './DEG/bulk/tumor pre E_VS_NE.rnk',col_names = F)

GSEA_batch.rnk(
  GSEA_installed_path = "../R tools/GSEA_4.3.2",
  DEG_path='DEG/bulk',
  species = 'human',
  gene_sets =c(`hallmark gene sets`='h.all.v2025.1.Hs.symbols.gmt'),
  symbol_chip='Human_Gene_Symbol_with_Remapping_MSigDB.v2025.1.Hs.chip',
  out_dir='GSEA',
  GSEA_plots_number=30,
  collapse='Collapse'
)

GSEA_bubble_3('GSEA')


## Fig S1B NEDD8 level in Pre tumor vs number of expanded T cell clone 
## 
## 
## 
T.expand <- read_csv2("TCR/T-cell_expansion.csv")
exp <- get_average_expression(Tumor,'NEDD8',assay = 'SCT',group.by=c('BC_type','expansion','timepoint'),remove_outliers_by_group.by = F) %>% filter(timepoint == 'Pre')

exp$patient_id <- exp$patient_id %>% str_replace_all('-','_')
plot.dt <- merge(exp,T.expand,by='patient_id')


ggplot(plot.dt,aes(x=value,y=`Number of expanded T cell clone`,color=expansion))+
  geom_point(aes(shape=BC_type)) +
  geom_smooth(data = plot.dt %>% filter(expansion=='E'), aes(),method = 'lm',se=F,linetype='dashed')+
  geom_smooth(data = plot.dt %>% filter(expansion=='NE'), method = 'lm',se=F,linetype='dashed')+
  #stat_poly_line() +
  #stat_poly_eq(use_label(c("eq", "R2"))) +
  ylab('Avg NEDD8 expression in tumor')+
  xlab('NEDD8 expression in PRE tumor')+
  stat_cor(method="pearson")








