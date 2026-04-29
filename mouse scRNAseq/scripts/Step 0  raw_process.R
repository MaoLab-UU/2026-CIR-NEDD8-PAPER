library(tidyverse)
library(readxl)
library(Seurat)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 32000*1024^2)
lapply(list.files("E:/bioinfo/R tools",full.names = T),source) %>% invisible()

### create raw seurat object
addr <- './raw'

obj <- Read10X(addr)

###### obly keep annotated cells

decode <-  fromJSON('cells_per_tag.json') %>% melt()
obj <- obj[,which(colnames(obj) %in% decode[,1])]

### rename deprecated genes
deprecated.genes.index <- grepl('DEPRECATED_',rownames(obj)) %>% which()

new.names <- rownames(obj)[deprecated.genes.index] %>% 
  str_remove('DEPRECATED_') %>%
  ### renamed deprecated genes
  gconvert(organism="mmusculus",target="ENTREZGENE",filter_na = F) %>%
  dplyr::select('target') %>% unlist() %>% as.vector()

rownames(obj)[deprecated.genes.index] <- new.names 

# ### convert gene alias to official symbol
new.names <- alias2SymbolTable(rownames(obj), species = "Mm")
use.names <- ifelse(is.na(new.names),rownames(obj),new.names)
rownames(obj) <- use.names

### remove dupliacted genes , if duplicated, take the maximum
no.redundant.obj <- obj[which(!(rownames(obj) %in%
                             use.names[duplicated(use.names)])),]

redundant.obj <- obj[which((rownames(obj) %in%
                               use.names[duplicated(use.names)])),] 


extra_dt <- data.frame()
for (i in unique(rownames(redundant.obj))) {
  j <- redundant.obj[which(rownames(redundant.obj) == i),]
  j <- apply(j, 2, FUN = max) %>% as.data.frame() %>%t()
  rownames(j) <- i
  
  if(nrow(extra_dt)==0){
    extra_dt <-j
  }else{
    extra_dt <- rbind(extra_dt,j)
  }}

extra_dt <- Matrix::Matrix(data=extra_dt %>% as.matrix(),sparse = T)
no.redundant.obj <- rbind(no.redundant.obj,extra_dt)

obj <- no.redundant.obj 

###
rm(no.redundant.obj)
rm(redundant.obj)

#### create meta data
meta <- read_xlsx('../metadata.xlsx')

meta.data <- data.frame(sample_id = decode$L1[order(decode$value,colnames(obj))])
meta.data$genotype <- lapply(meta.data$sample_id,function(x)meta$Genotype[which(meta$Sample==x)])%>% unlist()
meta.data$treatment <- lapply(meta.data$sample_id,function(x)meta$Treatment[which(meta$Sample==x)])%>% unlist()
meta.data$group <- paste(meta.data$genotype,meta.data$treatment,sep = '_')

### creat seurat object
obj <- CreateSeuratObject(counts = obj,meta.data = meta.data)

obj <- obj %>% subset(sample_id != 'BC011')
### add cell cycle data AND QC
obj <- QC(obj,species = 'mm')

obj <-  NormalizeData(obj)

#DotPlot(obj,c('Itgam','Ptprc','Cd14'))

obj <- cell_cycle_score(obj,species = 'mm')
#obj@assays$RNA@layers$data <- NULL


obj <- QC_plot(obj,filtering = T)

#write_rds(obj,'all_cells_raw.rds')












