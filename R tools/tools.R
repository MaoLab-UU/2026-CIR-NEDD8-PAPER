library(Seurat)
library(readr)
library(readxl)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(AUCell)
library(GSEABase)
library(lattice)
library(gprofiler2)
library(ggnewscale)
library(CellChat)
library(ggrepel)
#library(ComplexHeatmap)
#library(ggpubr)
library(ggrepel)
#library(glmGamPoi)
library(harmony)

######################
pct_plot <- function(obj,
                     celltype = "celltype",
                     sample_id = 'sample_id',
                      group = 'group',
                     significance = 'p', ## or 'p.adj',
                     show.p.value = F,
                     compare=NA){
  
  
  
  c <- obj@meta.data[,c(celltype,group,sample_id)]  
  colnames(c) <- c('celltype','group','sample_id')
  
  
  c <- c %>% group_by(group,sample_id,celltype) %>% summarise (n = n()) %>%
    mutate(Freq = n / sum(n) * 100)
  
  
  c$celltype <- factor(c$celltype,levels = unique(c$celltype))
  c$group <- factor(c$group,levels = levels(obj@meta.data[[group]]))
  
  if(length(compare) == 0){
    stat.test <- c %>%
      group_by(celltype) %>%
      wilcox_test(Freq ~ group) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
    
  }else{
    stat.test <- c %>%
      group_by(celltype) %>%
      wilcox_test(Freq ~ group,comparisons = compare) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")%>%
      add_significance("p")
  }

  
  
  # Create a box plot
  bxp <- ggboxplot(
    c, x = 'celltype', y = "Freq", 
    color = "group"
  )
  
  # Add p-values onto the box plots
  stat.test <- stat.test %>%
    add_xy_position(x = 'celltype', dodge = 0.8) 
  
  signif <- ifelse(show.p.value,'',".signif")

  
  bxp +  
    geom_point(aes(x=celltype,y=Freq,colour=group),
               position = position_jitterdodge(jitter.width=0.2) ) + 
    stat_pvalue_manual(
      stat.test,  label = paste0(significance,signif),size = 3.5
    ) # # Add 10% spaces between the p-value labels and the plot border
  #   + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  
}

####### Downsampling before run DEG or avgerage expression comparison

Prep_DEG_downsample <- function(obj=epi,
                                mode='median',
                                group.id = 'group',
                                sample_id = 'sample_id',
                                exclude.outliers = T,
                                cell.id='cell',
                                remove.cutoff.pct = 0.2){
  
  tab <- data.frame(obj@meta.data[,c(sample_id,cell.id,group.id)]) 
  colnames(tab) <- c('sample_id','cell.id','group.id')
  
  med <- tab %>% dplyr::count(sample_id, group.id) %>% group_by(group.id) %>% summarise(median=median(n))

  tab_downsmapled <- lapply(unique(med$group.id),function(x) {
    med.num <- med$median[med$group.id==x]
    
    tab.sub <- tab %>% filter(group.id == x) %>% group_by(sample_id) %>% 
      slice_sample(n = as.integer(med.num) ,replace = F)
    
    if(exclude.outliers){
      freq <- tab.sub %>% dplyr::count(sample_id)
      
      outliers <-   freq$sample_id[which(freq$n < as.integer(med.num*remove.cutoff.pct))]
      
      tab.sub <- tab.sub %>% filter(sample_id %in% outliers == F)
    }}
    ) 
  
  tab_downsmapled <- do.call('rbind',tab_downsmapled)
  #tab_downsmapled%>% count(sample_id, group.id)
  
  ### any sample less than 20% of the median will be excluded from DEG

  obj$cell.id <- obj@meta.data[,cell.id]

  return(obj %>% subset(cell.id %in% tab_downsmapled$cell.id))
}




avg_heatmap <-function(obj=obj,
                       modules=names(signatures.mm),
                       group_by='sample_id',
                       cols=c('#5dade2','#EEEEEF','#e25d5d'),
                       label = F,
                       x.order = NULL,
                       y.order = NULL,
                       return.dt = F,
                       wrap_y = NULL){
  
  means<- obj@meta.data[,modules]
  
  #means$group <- obj.cd8$group
  #means$group <- paste0(obj.cd8$sample_id,obj.cd8$group)
  means$group <- obj[[group_by]]
  
  
  
  means <- means  %>% arrange(group)
  
  plot_dt <- data.frame()
  for ( g in unique(means$group)%>%unlist()){
    avg_score <- means[which(means$group == g),modules] %>% colMeans() %>% as.data.frame()
    colnames(avg_score) <- g
    if(length(plot_dt)==0){
      plot_dt <-avg_score
    }else{
      plot_dt <-cbind(plot_dt,avg_score)
    }
  }
  if(return.dt == F){
    plot_dt <- t(apply(plot_dt , 1, rescale, to=c(-1, 1))) %>% as.data.frame()
  }
  
  plot_dt$sig <- rownames(plot_dt)
  
  
  plot_dt <- plot_dt %>% melt(id.vars = 'sig')
  
  plot_dt$sig <- factor(plot_dt$sig,levels = modules)
  
  if(x.order %>% length() != 0){
    plot_dt$variable <- factor(plot_dt$variable,levels = x.order)
  }
  
  if(y.order %>% length() != 0){
    plot_dt$y.order <- factor(plot_dt$sig,levels = y.order) %>% order()

  

    if(length(wrap_y)!=0){
      plot_dt$sig <- str_replace_all(plot_dt$sig,"_",' ')
      plot_dt$sig <- str_wrap(plot_dt$sig,width = wrap_y,
                              indent = 2,whitespace_only = T)
      
      plot_dt$sig <- factor(plot_dt$sig, levels = plot_dt$sig[plot_dt$y.order] %>% unique())
      
    }else{
      plot_dt$sig <- factor(plot_dt$sig,levels = y.order)
    }}

  
   
  
  if(label == T){
    p <- ggplot(plot_dt,aes(x=variable,y=sig,fill=value))+
      geom_tile()+
      geom_text(aes(label = round(value, 3)),size = 3) +
      scale_fill_gradientn(colors = cols)+
      labs(fill='scaled_avg_score')
  }else{
    p <- ggplot(plot_dt,aes(x=variable,y=sig,fill=value))+
      geom_tile()+
      scale_fill_gradientn(colors = cols)+
      labs(fill='scaled_avg_score')
  }
  
 if(return.dt){
   return(plot_dt)
 }else{
   return(p)
 }
  
}


#######################

calculate_sig <- function(obj,signature,h2m_conversion=T,split.by=NULL,ctrl=100){
  signature <- as.list(signature)
  
  signature <- sapply(signature,function(x) as.vector(na.omit(x) ))


  
  if(length(split.by) != 0){
    signature <- sapply(signature,function(x) strsplit(x,split.by)%>%unlist())
  }
 
    
  if(h2m_conversion == T){
    signature.mm <- sapply(signature,function(x) gene.hs.to.mm(x)$gene.m %>% na.omit %>% as.vector)
  }else{
    signature.mm <- signature
  }
  
  # intersect genes and remove length == 0 signatures
  signature.mm <- sapply(signature.mm,function(x) intersect(x,Features(obj)))
  signature.mm <- signature.mm[lapply(signature.mm,length)>0]
  print(signature.mm)
  try(obj@meta.data[,names(signature)]<-NULL,silent = T)
  
  
  obj <- AddModuleScore(obj,features = signature.mm,ctrl = ctrl)
  
  colnames(obj@meta.data)[(length(colnames(obj@meta.data))-length(names(signature.mm))+1)
                          :length(colnames(obj@meta.data))] <- names(signature.mm) %>% str_remove_all('\\d+$')
  
  return(list(obj,signature.mm))
  
}











#### set the pattern to ^mt- for mouse ,default is human

QC <- function(obj,species='mm'){
  obj$log10FeaturePerlog10UMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  if(species == 'mm'){
    print('Mitochondial genes detected:')
    print( grep(pattern = "^mt-", x = rownames(obj@assays[["RNA"]]), value = TRUE))
    
    print('Hemogloblin genes detected:')
    print( grep(pattern = "^Hb[^(p)]", x = rownames(obj@assays[["RNA"]]), value = TRUE))
    
    print('Heatshock protein detected:')
    print( grep(pattern = "^Hsp", x = rownames(obj@assays[["RNA"]]), value = TRUE))
    
    obj <- obj %>% 
      PercentageFeatureSet( "^mt-",col.name = 'mitoRatio')  %>% 
      PercentageFeatureSet( "^Hb[^(p)]",col.name = 'HbRatio')  %>% 
      PercentageFeatureSet( "^Hsp",col.name = 'HspRatio')  
    
  }
  if(species == 'hs'){
    print('Mitochondial genes detected:')
    print( grep(pattern = "^MT-", x = rownames(obj@assays[["RNA"]]), value = TRUE))
    
    print('Hemogloblin genes detected:')
    print( grep(pattern = "^HB[^(P)]", x = rownames(obj@assays[["RNA"]]), value = TRUE))
    
    print('Heatshock protein detected:')
    print( grep(pattern = "^HSP", x = rownames(obj@assays[["RNA"]]), value = TRUE))
    
    obj <- obj %>% 
      PercentageFeatureSet( "^MT-",col.name = 'mitoRatio')  %>% 
      PercentageFeatureSet( "^HB[^(P)]",col.name = 'HbRatio')  %>% 
      PercentageFeatureSet( "^HSP",col.name = 'HspRatio')  
    
  }
 
  return(obj)
}

QC_plot <- function(obj,filtering = F,name = 'All'){


  if(ncol(obj) > 2000){
    plot_obj <- obj[,sample(obj %>% colnames(),2000)]
  }else{
    plot_obj <- obj
  }
  
    nCount_high_hold <- quantile(obj$nCount_RNA, 0.95)
    nCount_low_hold <- quantile(obj$nCount_RNA, 0.05)
    
    nFeature_RNA_high_hold <- quantile(obj$nFeature_RNA, 0.95)
    nFeature_RNA_low_hold <- quantile(obj$nFeature_RNA, 0.05)
    
    log10GenesPerUMI_low_hold <- quantile(obj$log10FeaturePerlog10UMI,0.05)

    mitoRatio_high_hold <- quantile(obj$mitoRatio,0.95) %>% min(15) %>% max(10)
    mitoRatio_low_hold <- quantile(obj$mitoRatio,0.05)
  
    Hsp_high_hold <- quantile(obj$HspRatio,0.95)
    Hb_high_hold <- quantile(obj$HbRatio,0.95)
  
    meta <- plot_obj@meta.data[,c('mitoRatio','nFeature_RNA','nCount_RNA','log10FeaturePerlog10UMI','HspRatio','HbRatio')] %>% melt()

    plot <- ggplot(meta,aes(x=variable,y=value))+ theme_bw()+
     
    facet_wrap(~variable, scales="free") + geom_violin(fill='#609fd6', width = 0.6)+ geom_point(position = position_jitter(width=0.2),size=0.1) + 
      
     geom_hline(data=filter(meta,variable=='mitoRatio'),aes(yintercept=mitoRatio_high_hold),linetype='dashed',col='#d66060') +
      geom_hline(data=filter(meta,variable=='mitoRatio'),aes(yintercept=mitoRatio_low_hold),linetype='dashed',col='#d66060') +
     geom_text(data=filter(meta,variable=='mitoRatio'), aes( x=0.8, y=mitoRatio_high_hold * 1.1, 
                                                              label=paste0( '<=',mitoRatio_high_hold %>% round(3),'& | >=',mitoRatio_low_hold%>% round(3),'|',filter(meta,variable=='mitoRatio') %>% 
                                                                            filter(value <= mitoRatio_high_hold & value >= mitoRatio_low_hold )%>%
                                                                            nrow() / filter(meta,variable=='mitoRatio') %>% nrow() * 100,
                                                                            '%')),size=3)  +
      
      
     geom_hline(data=filter(meta,variable=='nFeature_RNA'),aes(yintercept=nFeature_RNA_low_hold),linetype='dashed',col='#d66060') +
     geom_hline(data=filter(meta,variable=='nFeature_RNA'),aes(yintercept=nFeature_RNA_high_hold),linetype='dashed',col='#d66060') +
     geom_text(data=filter(meta,variable=='nFeature_RNA'), aes( x=0.80, y=nFeature_RNA_high_hold * 1.1, 
                                                              label=paste0( '<=',nFeature_RNA_high_hold%>% round(3),'& | >=',nFeature_RNA_low_hold%>% round(3),'|',filter(meta,variable=='nFeature_RNA') %>% 
                                                                              filter(value >= nFeature_RNA_low_hold & value <= nFeature_RNA_high_hold)%>%
                                                                              nrow() / filter(meta,variable=='nFeature_RNA') %>% nrow() * 100,
                                                                            '%')),size=3)  +
      
      
     geom_hline(data=filter(meta,variable=='nCount_RNA'),aes(yintercept=nCount_low_hold),linetype='dashed',col='#d66060') +
      geom_hline(data=filter(meta,variable=='nCount_RNA'),aes(yintercept=nCount_high_hold),linetype='dashed',col='#d66060') +
      geom_text(data=filter(meta,variable=='nCount_RNA'), aes( x=0.8, y=nCount_high_hold * 1.1, 
                                                              label=paste0( '>=',nCount_low_hold%>% round(3),'& | >=',nCount_high_hold%>% round(3),'|',filter(meta,variable=='nCount_RNA') %>% 
                                                                              filter(value >= nCount_low_hold & value <= nCount_high_hold)%>%
                                                                              nrow() / filter(meta,variable=='nCount_RNA') %>% nrow() * 100,
                                                                            '%')),size=3)   +
    
      
     geom_hline(data=filter(meta,variable=='log10FeaturePerlog10UMI'),aes(yintercept=log10GenesPerUMI_low_hold),linetype='dashed',col='#d66060')+
      geom_text(data=filter(meta,variable=='log10FeaturePerlog10UMI'), aes( x=0.6, y=log10GenesPerUMI_low_hold * 0.9, 
                                                              label=paste0( '>=',log10GenesPerUMI_low_hold%>% round(3),'|',filter(meta,variable=='log10FeaturePerlog10UMI') %>% 
                                                                              filter(value >= log10GenesPerUMI_low_hold )%>%
                                                                              nrow() / filter(meta,variable=='log10FeaturePerlog10UMI') %>% nrow() * 100,
                                                                            '%')),size=3)  +
      
      geom_hline(data=filter(meta,variable=='HspRatio'),aes(yintercept=Hsp_high_hold),linetype='dashed',col='#d66060')+
      geom_text(data=filter(meta,variable=='HspRatio'), aes( x=0.6, y=Hsp_high_hold * 1.1, 
                                                                            label=paste0( '<=',Hsp_high_hold%>% round(3),'|',filter(meta,variable=='HspRatio') %>% 
                                                                                            filter(value <= Hsp_high_hold )%>%
                                                                                            nrow() / filter(meta,variable=='HspRatio') %>% nrow() * 100,
                                                                                          '%')),size=3)  +
      
      geom_hline(data=filter(meta,variable=='HbRatio'),aes(yintercept=Hb_high_hold),linetype='dashed',col='#d66060')+
      geom_text(data=filter(meta,variable=='HbRatio'), aes( x=0.6, y=Hb_high_hold * 1.1, 
                                                                            label=paste0( '<=',Hb_high_hold%>% round(3),'|',filter(meta,variable=='HbRatio') %>% 
                                                                                            filter(value <= Hb_high_hold )%>%
                                                                                            nrow() / filter(meta,variable=='HbRatio') %>% nrow() * 100,
                                                                                          '%')),size=3) 
      
      dt <- obj@meta.data[,c('nFeature_RNA','nCount_RNA')]
      p2 <- ggplot()+ theme_classic()+
      geom_point(data=dt,aes(x=nCount_RNA,y=nFeature_RNA),col='black',alpha=0.3)+
      geom_hline(aes(yintercept=nFeature_RNA_low_hold),linetype='dashed',col='#d66060')+
      geom_hline(aes(yintercept=nFeature_RNA_high_hold),linetype='dashed',col='#d66060')+
      geom_vline(aes(xintercept=nCount_low_hold),linetype='dashed',col='#d66060')+
      geom_vline(aes(xintercept=nCount_high_hold),linetype='dashed',col='#d66060')+
      geom_text( aes( x=nCount_low_hold*1.05,y=nFeature_RNA_low_hold*1.05,label=paste0( 'in area',filter(dt,nCount_RNA <= nCount_high_hold & nCount_RNA >= nCount_low_hold) %>% 
                                                                             filter(nFeature_RNA <=nFeature_RNA_high_hold & nCount_RNA >= nFeature_RNA_low_hold )%>% unlist() %>%
                                                                             length() / dt  %>% unlist() %>%
                                                                             length() * 100,
                                                                             '%')),size=3,family='Arial')
  options(warn=-1)   
  pdf(paste0(name,' ','QC_plot.pdf'),width = 21,height = 7)
  print(plot+p2)
  dev.off()
      
  if(filtering == T){
    obj <- cell_filter(obj,
                          nCount_high_hold,
                          nCount_low_hold,
                          
                          nFeature_RNA_high_hold,
                          nFeature_RNA_low_hold,
                          
                          log10GenesPerUMI_low_hold,
                          
                          mitoRatio_high_hold,
                          mitoRatio_low_hold,
                          
                          Hsp_high_hold,
                          Hb_high_hold)
    return(obj)
  }else{
    return(obj) 
  }    
      
  
  
}

cell_filter <- function(obj,
                        nCount_high_hold,
                        nCount_low_hold,
                        
                        nFeature_RNA_high_hold,
                        nFeature_RNA_low_hold,
                        
                        log10GenesPerUMI_low_hold,
                        
                        mitoRatio_high_hold,
                        mitoRatio_low_hold,
                        
                        Hsp_high_hold,
                        Hb_high_hold){
  filtered_seurat <- subset(x = obj, 
                            subset= ( nCount_RNA >= nCount_low_hold ) & 
                              ( nCount_RNA >= nCount_low_hold ) &
                              ( nFeature_RNA >= nFeature_RNA_low_hold) & 
                              ( nFeature_RNA <= nFeature_RNA_high_hold) & 
                              (log10FeaturePerlog10UMI > log10GenesPerUMI_low_hold) & 
                              (mitoRatio <= mitoRatio_high_hold) &
                              (mitoRatio >= mitoRatio_low_hold) &
                              (Hb_high_hold <= Hb_high_hold) &
                              (Hsp_high_hold <= Hsp_high_hold)
                            )
  return(filtered_seurat)
}



find_optimal_pcs <- function(obj,reduction='pca'){
  
  #obj <- obj %>% RunPCA(nfeatures.print = 1)
  # determine optimal number of PC
  stdv <- obj[[reduction]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  return(min_pc)
}


Estimate_multiplet_rate <- function(number_of_cell_recovered){
  multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                    'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                    'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))

  multiplet_rate <- multiplet_rates_10x %>% 
    dplyr::filter(Recovered_cells < number_of_cell_recovered) %>% 
    dplyr::slice(which.max(Recovered_cells)) %>% 
    dplyr::select(Multiplet_rate) %>% as.numeric(as.character())
  return(multiplet_rate)
  }


Identify_doublets <- function(obj,resolution = 0.3,sct = F){
  
  if(!sct){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
    obj <- ScaleData(obj,features=Features(obj))
    
  }
  obj <- obj %>% RunPCA(nfeatures.print = 10)
  print('Make sure object is preprocessed (done normalization / PCA)') 
  
  min_pc <- find_optimal_pcs(obj)
  # calculate optimal pk pn
  sweep_list <- paramSweep(obj, PCs = 1:min_pc, sct = sct)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  
  obj <- RunUMAP(obj, dims = 1:min_pc)
  obj <- FindNeighbors(object = obj, dims = 1:min_pc)              
  obj <- FindClusters(object = obj, resolution = resolution)
  
  annotations <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(Estimate_multiplet_rate(nrow(obj@meta.data)) * nrow(obj@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  obj <- doubletFinder(seu = obj, 
                       PCs = 1:min_pc, 
                       pK = optimal.pk,
                       nExp = nExp.poi.adj,
                       sct = sct)
  metadata <- obj@meta.data
  colnames(metadata)[colnames(metadata) %>% length] <- "doublet_finder"
  obj@meta.data <- metadata 
  
  p <- DimPlot(obj,reduction = 'umap',group.by = 'doublet_finder')+
    ggtitle('Doublet ratio % : ',length(grep('Doublet',obj@meta.data$doublet_finder))/nrow(obj@meta.data) * 100 )
  
  png(paste0(unique(obj$sample_id),'_doublets.png'),width=6,height = 6,units = 'in',res = 800)
  print(p)
  dev.off()
  # subset and save
  #obj <- subset(obj, doublet_finder == "Singlet")
  return(metadata)
}

standard_process <- function(obj,normalization = F,var.to.regress=c()){
  if(normalization == T){

    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj,features=Features(obj),vars.to.regress = var.to.regress)
  }

  obj <- RunPCA(obj, features = VariableFeatures(obj = obj))
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 2)
  obj <- RunUMAP(obj, dims = 1:20,
                    seed.use = 123)
  return(obj)
}

#obj=matrix.list[[1]]

RunSCT <- function(obj,
                   vars.to.regress=c(),
                   join = F,
                   split = c(),
                   variable.features.n = 3000,
                   clip.range = c(),
                   res = c(0.8,2)){
  
  DefaultAssay(obj) <- 'RNA'
  obj@reduction <- list()
  if('SCT' %in% Assays(obj)){
    obj[['SCT']] <- NULL
  }
  
  if(join == T){
    obj <- JoinLayers(obj)
  }
  if(length(split) > 0){
    obj <- split(obj,f=obj@meta.data[[split]])
  }
  
  if(length(clip.range) >0 ){
    obj <- SCTransform(obj,
                       vars.to.regress = vars.to.regress,
                       #clip.range = c(-5,5),
                       variable.features.n = variable.features.n,
                       verbose = TRUE,
                       vst.flavor = 'v2',clip.range = clip.range)
  }else{
    obj <- SCTransform(obj,
                       vars.to.regress = vars.to.regress,
                       #clip.range = c(-5,5),
                       variable.features.n = variable.features.n,
                       verbose = TRUE,
                       vst.flavor = 'v2')
  }
  
  obj <- obj %>%
    RunPCA(assay = "SCT", npcs = 50)
  
  pcs <- find_optimal_pcs(obj) + 15
  
  obj <- PrepSCTFindMarkers(obj)
  
  obj <- obj %>% 
    FindNeighbors(reduction = "pca", dims = 1:pcs,nn.method = "rann") %>%
    FindClusters( resolution = res) %>%
    RunUMAP(reduction = "pca", dims = 1:pcs, reduction.name = "umap.pca",seed.use = 123)
  
  return(obj)
}  




RunIntegration <- function(obj, harmony = T,pcs=30){
  if(harmony == T){
    
    obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                           orig.reduction = "pca", new.reduction = 'harmony',
                           assay = "SCT", verbose = FALSE)
    
    
    
    obj <- obj %>% 
      FindNeighbors(reduction = "harmony", dims = 1:pcs,nn.method = "rann") %>%
      FindClusters( resolution = c(0.8,2), cluster.name = c("harmony_clusters_0.8","harmony_clusters_2") ) %>%
      RunUMAP(reduction = "harmony", dims = 1:pcs, reduction.name = "umap.harmony",seed.use = 123)
    
  }else{
    
    
    obj <- IntegrateLayers(object = obj, method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = 'integrated.cca',
                           assay = "SCT", verbose = FALSE)  %>%
      
      #pcs <- find_optimal_pcs(obj,reduction = 'integrated.cca') + 15- 
      pcs <- npcs
      
      obj <- obj %>% 
        FindNeighbors(dims = 1:pcs, reduction = "integrated.cca",nn.method = "rann") %>%
        FindClusters(resolution = c(0.8,2), cluster.name = c("cca_clusters_0.8","cca_clusters_2")) %>%
        RunUMAP(dims = 1:pcs,seed.use = 123,reduction='integrated.cca')
      
  }
  return(obj)
}















cell_cycle_score <- function(obj,species = 'hs'){

  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  
  if(species == 'mm'){
    s.genes = gorth(s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    g2m.genes = gorth(g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
  }
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes)
  obj$cc.difference <- obj$S.Score - obj$G2M.Score
  
  return(obj)
}


DEG2RNK <- function(DEG,p_hold , log2fc_hold ,name){
  
  
  up <- DEG[which((DEG$avg_log2FC > 0)  & (abs(DEG$avg_log2FC) > log2fc_hold) & (DEG$p_val_adj < p_hold)) ,]
  
  down <- DEG[which ((DEG$avg_log2FC < 0) & (abs(DEG$avg_log2FC) > log2fc_hold) & (DEG$p_val_adj < p_hold)) ,]
  
  
  up <- up[order(up$avg_log2FC,decreasing = T),]
  
  down <- down[order(down$avg_log2FC,decreasing = F),]
  
  DEG <- rbind(up,down)
  
  rnk <- DEG[order(DEG$avg_log2FC,decreasing = T),]
  
  write_tsv(data.frame(row=rownames(rnk),log2fc=rnk[,'avg_log2FC']) ,paste0(name,'.rnk'),col_names = F)
  
  
  
}

calculate_AUCell <- function(obj=malignant_cell,geneset){
  
  cells_AUC <- sapply(names(geneset), function(x) {
      exprMatrix <- GetAssayData(obj) 
      exprMatrix <- as(exprMatrix, "dgCMatrix")
      
      geneSets <- GeneSet(geneset[[x]]%>%unlist(), setName=x)
      
      AUCell_run(exprMatrix, geneSets)
    })
  return(cells_AUC)
}



# deprecated in Seurat v5

# SCT_harmony <- function(obj_list=obj,
#                         splitted_by=c('patient_id','treatment_phase'),
#                         vars_to_regress='percent.mt',
#                         do.harmony=F,
#                         do.integration=F,
#                         variable.features.n = 2000){
#   
#   transfer <- function(obj_list){
#     if(length(obj_list) > 1 ){
#       merged_obj <- merge(obj_list[[1]],
#                           obj_list[2:length(obj_list)],
#                           merge.data = TRUE)
#     }else(
#       merged_obj <- obj_list
#     ) 
#     return(merged_obj)
#   }
# 
#   
#   if(do.harmony == T){
#   merged_obj <- transfer(obj_list) %>%
#     SCTransform(vars.to.regress = c(vars_to_regress), vst.flavor = "v2",variable.features.n = variable.features.n) %>%
#     RunPCA(assay = "SCT", npcs = 50)
#     harmonized_obj <- merged_obj %>%  
#       RunHarmony(merged_obj, 
#                  group.by.vars = splitted_by, 
#                  reduction = "pca", assay.use = "SCT", reduction.save = "harmony")%>%
#       FindNeighbors(reduction = "harmony") %>%
#       FindClusters(resolution = c(0.8,2)) %>%
#       RunUMAP(reduction = "harmony", assay = "SCT", dims = 1:20,seed.use = 123) 
#       
#     
# 
#     return(harmonized_obj)
#   }else{
#     if(do.integration == T) {
#       
#       merged_obj <- lapply(merged_obj, function(x) SCTransform(x,vars.to.regress = c(vars_to_regress), vst.flavor = "v2"))
#       features <- SelectIntegrationFeatures(obj.list = merged_obj, nfeatures = 3000)
#       merged_obj <- PrepSCTIntegration(obj.list = merged_obj, anchor.features = features)
#       merged_obj.anchors <- FindIntegrationAnchors(obj.list = merged_obj, normalization.method = "SCT",
#                                                    anchor.features = features)
#       merged_obj <- IntegrateData(anchorset = merged_obj.anchors, normalization.method = "SCT")
#     }
#     
#     merged_obj <- merged_obj%>%
#     RunPCA(assay = "SCT", npcs = 50)%>%
#     FindNeighbors( dims = 1:20, reduction = "pca")%>%
#     FindClusters(resolution = c(0.5,1,2)) %>%
#     RunUMAP(dims = 1:20,seed.use = 123,reduction='pca')
#      
#      
#     return(merged_obj)
#   }
#   
# }


remove_outliers <- function(dataframe,value='value',split.by='gene',group.by='group.by'){
  
  if(length(group.by) != 1){
    dataframe$`group.by` <- lapply(1:nrow(dataframe), function(x) paste0(dataframe[x,group.by],collapse = '_')) %>% unlist()
  }
  
  return_dt <- data.frame(matrix(ncol = ncol(dataframe), nrow = 0))
  colnames(return_dt) <- colnames(dataframe)
  
  dataframe <- dataframe %>% as.data.frame()
  
  for (gene in unique(dataframe[[split.by]])){
    
    subdt <- dataframe[which(dataframe[[split.by]] == gene),]
    outlier_value <- sapply(unique(subdt[[group.by]]),function(x)  boxplot.stats( subdt[[value]][which(subdt[[group.by]]==x)])$out)
    
    for (group in names(outlier_value) ){
      
      if(length(outlier_value[[group]] != 0)){
        return_dt <- rbind(return_dt,subdt[which(subdt[[group.by]] == group & subdt$value != outlier_value[[group]]),])
      }else{
        return_dt <- rbind(return_dt,subdt[which(subdt[[group.by]] == group),])
      }
      #index <- append(index,which(dataframe[[split.by]] %in% gene & dataframe[[group.by]] %in% group & dataframe[[value]] %in% outlier_value[[group]]))
      }}
  return(return_dt)}





## calculate average expression of gene(s) for each sample *
#(grouby_by = timepoint + patient id + xxx... sample conditions), 
# outlier_individuals are removed for each condition 
#which is basically grouping without considering patient_id
#Users can reduce the outlier factors, thus pooling one or several condition together 
#(like pooling the BC_type together so BC_type is not used in group.by)

get_average_expression <- function(obj=tnbc.tumor,
                                   features,
                                   assay,
                                   group.by,
                                   individual_labels='patient_id',
                                   reshape_for_ggplot = T,
                                   remove_outliers_by_group.by = T){
  
  features <- features[which(features %in% Features(obj,assay = assay))]
  
  avg_exp <- AverageExpression(obj,features = features,group.by=append(individual_labels,group.by),assays = assay)[[assay]] %>% as.data.frame()
  
  num <- 4+length(group.by)
  
  if(reshape_for_ggplot){
    if(length(features) != 1){
      avg_exp$gene <- rownames(avg_exp)
    }else{
      avg_exp$gene <- features
    }
    
    
    avg_exp <- avg_exp %>% melt(id.vars = 'gene' )
    avg_exp <- concat.split(avg_exp, "variable", "_")
    avg_exp$`group.by` <- lapply(1:nrow(avg_exp), function(x) paste0(avg_exp[x,5:num],collapse = '_')) %>% unlist()
    
    colnames(avg_exp)[4:num] <- append(individual_labels,group.by)
    
  }
  if(remove_outliers_by_group.by){
    
    avg_exp <- remove_outliers(avg_exp,split.by = 'gene' ,group.by='group.by')
  }
  
  return(avg_exp)
  
}


#### Users can rbind single gene get_average_expression and multiple genes get_average_expression
# dataframes and run the correlation plot function to see the expression correlation between that single gene and 
# each gene in the second multiple gene set.

# In development
# average_expression_correlation <- function(average_expression_obj,
#                                            method='t.test'){
# }


Run_cellchat <- function(SeuratObj = sub_obj,
                         assay = 'SCT',
                         manual_dt = NULL,
                         manual_meta = NULL,
                         sample_col = 'patient_id',
                         label_col = 'cellType',
                         ccDB = CellChatDB.human,
                         search = 'Cell-Cell Contact',
                         smooth = PPI.human,
                         spatial.coord=NULL){
  
  if(is.null(manual_dt) ){
  dt <- SeuratObj[[assay]]$data # normalized data matrix
  }else{
  dt <- manual_dt
  }
  
  if(is.null(manual_meta) ){
    
    if(length(sample_col)!=0 ){
      meta <- data.frame(labels = SeuratObj[[label_col]], row.names = colnames(SeuratObj),samples = SeuratObj[[sample_col]]) # create a dataframe of the cell labels
      colnames(meta) <- c('labels','samples')
    }else{
      meta <- data.frame(labels = SeuratObj[[label_col]], row.names = colnames(SeuratObj)) # create a dataframe of the cell labels
      colnames(meta) <- c('labels')
    }
    
  }else{
    meta <- manual_meta
  }
  
  
  
  if(length(spatial.coord) == 0){
    cc <- createCellChat(object = dt, meta = meta, group.by = 'labels')
  }
  
  if(length(spatial.coord) != 0){
    
    colnames(spatial.coord)[c(1,2)] <- c("imagerow", "imagecol")
    spatial.coord <- spatial.coord[,c(1,2)]
    
    #spatial.coord <- spatial.coord*10
    
    conversion.factor = 0.18 
    d = computeCellDistance(spatial.coord)
 
    spot.size = mean(d)*conversion.factor 
    print(spot.size)
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2 * 1.5)
    
    cc <- createCellChat(object = dt, meta = meta, group.by = "labels",
                          datatype = "spatial", coordinates = spatial.coord, 
                         spatial.factors = spatial.factors)
  }
  
  
  
  ccDB.use <- subsetDB(ccDB,search = search)
  cc@DB <- ccDB.use
  
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc,thresh.pc = 0.05)
  cc <- identifyOverExpressedInteractions(cc,variable.both = F)
  cc <- projectData(cc, adj = smooth)
  
  if(length(spatial.coord) != 0){
    cc <- computeCommunProb(cc, type = "triMean",raw.use = FALSE,
                            contact.range = spot.size*1.5,distance.use = TRUE,scale.distance = 36,
                            interaction.range =  spot.size*3.5)
  }else{
    cc <- computeCommunProb(cc, type = "triMean",raw.use = FALSE)
  }
  

  cc <- filterCommunication(cc, min.cells = 10)

  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  print('idents level :')
  print(as.character(cc@idents %>% levels())) 
  
  return(cc)
}



# density plot template:

# {
#   plot_dts <- split(cor_dt,cor_dt$gene.y)  
#   
#   plot.list <- lapply(plot_dts, function(x) ggplot(x,aes(x=value.x,y=value.y,group = expansion.x))+
#                         geom_point(aes(shape = BC_type.x,color=expansion.x) )+
#                         geom_density_2d_filled(aes(alpha=..level..,fill=expansion.x),bins=4)+
#                         scale_alpha_manual(values = c(0,0.05,0.15,0.3))+
#                         theme_bw()+
#                         theme(legend.position="none")+
#                         xlab('NEDD8')+
#                         ylab('Genes')
#   )
#   
#   ggarrange(plotlist=plot.list , widths = c(1,1), common.legend = T,labels=names(plot.list),
#             hjust = -0.5,
#             vjust = -0.5,font.label = list(size = 10, color = "black", face = "bold", family = NULL))
#   
#   
# }



# perason correlation plot template:

# {
#   ggplot(cor_dt,aes(x=value.x,y=value.y))+
#     facet_wrap(~gene.y,scales = 'free')+
#     geom_point(data = cor_dt, aes(color=expansion.x,shape=BC_type.x))+
#     stat_cor()+
#     geom_smooth(method='lm', se=FALSE,alpha=0.2)+
#     ylab('Genes')+
#     xlab('Nedd8')
# }


# box_plotter <- function(plot_dt){
#   
# }

legend_plot.list <- function(plot.list,extension,name,height.factor=3.5,width.factor=3.5){
  
  width = length(plot.list)^0.5 %>% ceiling
  
  p <- ggarrange(plotlist=plot.list , widths = c(1,1), common.legend = T,labels=names(plot.list),legend="right",
                 hjust = 0,
                 vjust = -0.5,font.label = list(size = 10, color = "black", face = "bold", family = NULL) ) +
    theme(plot.margin = margin(1,1,1,1, "cm")) 
  
  pdf(paste0(name,'_',extension,'_plot.pdf'),width = width * width.factor ,height = ((length(plot.list)/width) %>% ceiling) * height.factor )
  print(p)
  dev.off()
}


####################################### find_DEG_bewteen_groups (in specific clusters)

find_DEG_between_groups <- function(obj_reduced,
                                    subset_name='All_cluster',
                                    subset_cluster=c(),
                                    control_group=1,
                                    variable_group=2,
                                    logfc.threshold =DEG.logfc.threshold,
                                    min.pct=DEG.min.pct,
                                    save_folder='DEG'){
  obj_reduced@active.assay = 'RNA'  
  if(length(subset_cluster)!=0 ){
    if(subset_cluster[1] > 0){
      print(paste0('Looking for DEG in clusters : ',paste0(subset_cluster,collapse = ',') ))
      obj_reduced <- subset(reduced_data,idents = subset_cluster,invert = FALSE)
    }
    if( subset_cluster[1] < 0){
      clusters= unique(levels(reduced_data$seurat_clusters))
      print(paste0('Looking for DEG in clusters : ',paste0(clusters[clusters %in%abs(subset_cluster)==F],collapse = ',')))
      obj_reduced <- subset(reduced_data,idents = unique(reduced_data$seurat_clusters)[unique(reduced_data$seurat_clusters)%in%abs(subset_cluster)==F],invert = FALSE)
    }}else{
      print('Looking for DEG across all clusters ')
    }
  
  dir <- dir_create(save_folder,subset_name,'')
  if(is.numeric(control_group[1])==T){
    control_group=unique(obj_reduced$group)[control_group]
  }
  
  if(is.numeric(variable_group[1])==T){
    variable_group=unique(obj_reduced$group)[variable_group]
  }
  
  print(paste0('set : ', control_group,' as control') )
  print(paste0('set : ', variable_group,' as variable') )
  
  DEG <- FindMarkers(obj = obj_reduced ,ident.1=variable_group,ident.2=control_group,group.by='group',logfc.threshold =logfc.threshold ,min.pct=min.pct)
  
  write.csv(DEG,paste0(dir,'/',subset_name,'_group_',paste0(variable_group,collapse='_'),'_vs_group_',paste0(control_group,collapse='_'),'.csv'))
  #return(DEG)
}


#############################

GSEA_bubble <- function(GSEA_folder='GSEA',
                        height_factor=1,
                        width_factor=1,
                        GSEA_fdr_hold=0.5,
                        fdr_top=20){
  files <- list.files(GSEA_folder,pattern='gsea_report_for',recursive = T,full.names = T)
  files <- data.frame(file_name=grep('.tsv',files,value = T))
  files$ident <- lapply(files$file_name,function(x) str_sub(x,-17,-1)) %>% unlist() %>% invisible()
  save_dirs <- files$file_name %>% dirname() %>% dirname()
  
  batch <- function(iden=files$ident[3],fl = files$file_name[3]){
    print(paste0('Launch GSEA plot for ' ,dirname(files$file_name[files$ident==iden]) %>% basename() %>% unique()) )
    data <- lapply(files$file_name[files$ident==iden],function(x) read.table(x,sep="\t",quot="",header = T)%>%invisible())
    
    plot <- function(dt=data[[1]]){
      if(nrow(dt)==0){
        return('Ship empty data')
      }
      if(dt$NES[1]<0){
        group='Suppressed'
      }else{
        group='Activated'
      }
      
      dt <-  dt[dt$`FDR.q.val`<GSEA_fdr_hold,]
      dt <- dt[order(dt$`FDR.q.val`),]
      if(nrow(dt) > fdr_top){
        dt <- dt[1:20,]
      }
      
      dt$NAME <- gsub('_',' ',dt$NAME)
      dt$NAME <- lapply(dt$NAME,function(x) str_extract(x," .*")%>% tolower()) %>% unlist() 
      
      dt$NAME <- str_wrap(dt$NAME, width = 40,  indent = 2,whitespace_only = T)
      dt$NES <- abs(dt$NES)
      dt$NAME <- factor(dt$NAME,levels = dt$NAME[order(dt$NES)])
      
      ggplot(dt,aes(x=NES,y= NAME,size=SIZE,color=`FDR.q.val`))+
        geom_point()+
        theme_bw()+
        scale_color_continuous(high='#42C0EF',low='#F37E00')+
        ylab(NULL)+
        xlab('Absolute NES')+
        ggtitle(paste0(basename(fl),'_GSEA_enrichment_',group))+
        theme(axis.text.y=element_text(size=10))
      
      ggsave(paste0(strsplit(dirname(fl),'.Gsea')[[1]][1],'_',group,'.pdf'),height = max(7,0.6*nrow(dt))*height_factor,width = 8*width_factor)
    }
    lapply(1:length(data),function(x) plot(data[[x]]) )
    
  }
  lapply(1:nrow(files),function(x) try(batch(files$ident[x],files$file_name[x])))
  return('GSEA plot finished ....')
}




volcano <- function(DEG_path='DEG',
                    p_hold = 0.05,
                    log2_fc_hold = 0.2,
                    tops=10,
                    highlight_top_by='avg_log2FC',
                    highlight_by_keys=F,
                    width=7,
                    height=7)
                    {
  if(highlight_top_by %in% c('p_val_adj','avg_log2FC')==F){
    return('Plot failed, avaliable highlight_top_by options : p_val_adj,avg_log2FC')
  }
  
  folders <- list.dirs(DEG_path,recursive = F)
  
  
  
  
  batch <- function(path=folders[2]){
    csv_file <- list.files(path,recursive=T,pattern = '.csv',full.names = T)
    
    highlight_keys <- list.files(path,recursive=T,pattern = 'highlight_key.txt',full.names = T)
    
    cat(paste0('Plotting: ',basename(csv_file),'\n'))
    
    data <- read.csv2(csv_file,row.names = 1)
    #data <- data[data$p_val_adj<p_hold,]
    
    data$significance <- 'Stable'
    data$significance[data$avg_log2FC>log2_fc_hold & data$p_val_adj<p_hold] <- 'Up'
    data$significance[data$avg_log2FC< -log2_fc_hold & data$p_val_adj<p_hold] <- 'Down'
    data$gene <- rownames(data)
    

    Colors <- c(Down='#42C0EF',Stable='grey',Up='#F37E00')
    
    
    if(length(highlight_keys) != 0 & highlight_by_keys==T){
      highlight_keys <- read.table(highlight_keys,col.names = F) %>% unlist()
      
      data_highlight <- data[highlight_keys,] %>% na.omit()
      tops <- length(highlight_by_keys)
      if(nrow(data_highlight)==0){
        return('All of highlight_key genes were not found in DEG')
      }
      highlight_by = 'keys'
      
    }else{

      if(length(highlight_keys) == 0 & highlight_by_keys ==T){
        cat('"highlight_key.txt" file not found, use tops instead\n')
      }
      
      if(highlight_top_by=='p_val_adj'){
        data_highlight <- rbind(top_n(data[data$significance=='Up',],-tops,
                                      data[data$significance=='Up',][,highlight_top_by]),
                                
                                top_n(data[data$significance=='Down',],-tops,
                                      data[data$significance=='Down',][,highlight_top_by])) 
      }else{
        data_highlight <- rbind(top_n(data[data$significance=='Up',],tops,
                                      data[data$significance=='Up',][,highlight_top_by]),
                                
                                top_n(data[data$significance=='Down',],-tops,
                                      data[data$significance=='Down',][,highlight_top_by])) 
      }
      highlight_by = highlight_top_by
    }
    
    data$significance <- factor(data$significance,levels = c('Up','Stable','Down'))
    plot<-ggplot(data,aes(y=-log10(p_val_adj),x=avg_log2FC,color=significance))+
      theme_classic()+
      
      geom_point(size=1,alpha=0.15)+
      scale_color_manual(values=Colors)+
      #geom_vline(xintercept =log2_fc_hold,linetype='dotted',alpha=0.5 )+
     # geom_vline(xintercept =-log2_fc_hold,linetype='dotted',alpha=0.5)+
      #geom_hline(yintercept =-log10(p_hold ),linetype='dotted',alpha=0.5 )+
      geom_point(data=data_highlight,aes(x=avg_log2FC,y=-log10(p_val_adj)),colour='#f3e529',size=0.4,fill='white',alpha=0.7)+
      geom_text_repel(data=data_highlight,max.overlaps = Inf,
                      aes(x=avg_log2FC,y=-log10(p_val_adj),color=significance,label = gene),
                     show.legend = FALSE,size=3,
                      nudge_x = ifelse(data_highlight$significance == 'Up',unit(0.4,'pt'),unit(-0.4,'pt')),
                     segment.size = 0.1)+
      geom_vline(xintercept = c(-log2_fc_hold,log2_fc_hold),linetype='dashed',alpha=0.5 )+
      geom_hline(yintercept = -log10(p_hold),linetype='dashed',alpha=0.5)+
      theme_light()
    
    
    ggsave(paste0(strsplit(csv_file,'.csv')[[1]][1],'_',highlight_by,'_volcano.png') ,plot,width = width,height=height)
    
    return('Volcano plot Complete\n')}
  
  for (i in folders){cat(paste0(batch(i),'\n'))}
}


process_infercnv <- function(cnv_addr = cnv_objs_rds[[1]]){
  cnv_obj <- read_rds(cnv_addr)
  
  observations <- cnv_obj@expr.data[,cnv_obj@observation_grouped_cell_indices%>%unlist()] %>% t() %>% 
    as.data.frame() 
  
  observations$cell <- 'observations'
  
  references <- cnv_obj@expr.data[,cnv_obj@reference_grouped_cell_indices%>%unlist()] %>% t() %>% 
    as.data.frame() 
  
  references$cell <- 'references'
  
  
  all <- rbind(observations,references)
  
  all$totalCNV <- rowMeans(all[,1:(dim(all)[2]-1)]^2)

  # then generate correlation
  avg <- colMeans(all[,1:(dim(all)[2]-2)]) %>% as.data.frame()
  
  rownames(avg) == rownames(all[,1:(dim(all)[2]-2)] %>% t())
  
  CNV_cor <- cor(avg,all[,1:(dim(all)[2]-2)] %>% t(),use = "pairwise.complete.obs")%>%as.data.frame()
  CNV_cor[2,] <- all[colnames(CNV_cor),'totalCNV']
  
  CNV_cor <- CNV_cor %>% t() %>% data.frame()
  
  colnames(CNV_cor) <- c('CNV_correlation','CNV_score')
  
  CNV_cor$cell <- all[rownames(CNV_cor),'cell']
  
  CNV_cor$malignancy <- lapply(1:nrow(CNV_cor),function(x) ifelse(( 
                                                                      CNV_cor$CNV_score[x] <= CNV_score_hold &
                                                                      CNV_cor$CNV_correlation[x] <= CNV_correlation_hold )
                                                                     ,'Non-Malignant','Malignant'))%>%unlist()

  CNV_cor$sample <- str_extract(string =  cnv_addr ,pattern =  'BC[:digit:]*')
  
  return(CNV_cor)
  
}




FindMks_Volcano <- function(DEG,
                            p_adj.hold = 0.01,
                            avg_lfc.hold = 0.5,
                            top_n_plot = 10,
                            gene.highlight = '',
                            log2fc = 'avg_log2FC',
                            p_val = 'p_val_adj',
                            show.tops=T){
  
  DEG$avg_log2FC <- DEG[[log2fc]]
  DEG$p_val_adj <- DEG[[p_val]]

  
  DEG$gene <- rownames(DEG)
  DEG$group <- ifelse(DEG$p_val_adj < p_adj.hold & DEG$avg_log2FC > avg_lfc.hold,'Up',
                      ifelse(DEG$p_val_adj < p_adj.hold & DEG$avg_log2FC < -avg_lfc.hold,'Down',
                             'Stable'))
  DEG$group <- factor(DEG$group,levels = c('Up','Stable','Down'))
  DEG$label <- NA
  if(show.tops){
    n = top_n_plot
    
    DEG$p_val_adj_z <- rescale(-log10(DEG$p_val_adj),to=c(0,1))
    DEG$avg_log2FC_z <- rescale(DEG$avg_log2FC,to=c(-1,1)) 
    DEG$dist <- DEG$p_val_adj_z**2 + DEG$avg_log2FC_z**2
    
    highlights <- rbind(DEG %>% filter(group=='Up') %>% top_n(n,dist),
                        DEG %>% filter(group=='Down') %>% top_n(n,dist)) %>% rownames()

    DEG[highlights,'label']  <- DEG[highlights,'gene'] 
  }
 
  
  if(length(gene.highlight) >0){
    DEG[gene.highlight,'label'] <- gene.highlight
  }

  
  volcano <- ggplot(DEG,aes(x=avg_log2FC,y=-log10(p_val_adj),colour = group))+
    geom_point()+
    geom_hline(yintercept = -log10(p_adj.hold),linetype ='dashed',alpha=0.2)+
    geom_vline(xintercept = c(-avg_lfc.hold,avg_lfc.hold),linetype ='dashed',alpha=0.2)+
    geom_text_repel(max.overlaps = Inf,
                    aes(label = label),show.legend = FALSE,size=2.5,colour='#2a2a2a',min.segment.length = 0)+
    scale_color_manual(values = c('#f56969','grey','#69bbf5'))+
    theme_bw()
  return(volcano)
}


                               

DEG_pipeline <- function(obj,
                         GSEA_installed_path = "E:/bioinfo/GSEA_4.3.2",
                         comparisons=c('Jak2 KO_IgG2a','Jak2 Nedd8 KO_IgG2a'),
                         DEG_out_dir = 'tumor analysis/DEG',
                         GSEA_out_dir = 'tumor analysis/GSEA',
                         species = 'mouse' ### mouse or human
                         ){
  
  mks <- FindMarkers(obj,ident.1 = comparisons[1],ident.2 = comparisons[2])
  mks <- filter(mks, pct.2 > 0.05 )
  
  #### if comparison name contain any non ASCII words convert them
  #comparisons <- comparisons   %>% str_replace_all('α','a')
  comparisons <- comparisons %>% iconv("UTF-8", "ASCII", sub = "")
  
  
  DEG_out_dir.new <- paste0(DEG_out_dir,'/',paste0(comparisons,collapse = ' vs '))
  if(! dir.exists(DEG_out_dir.new)){
    dir.create(DEG_out_dir.new)
  }
  
  ########## set minium p_val to 10 e-300
  mks$p_val_adj <- mks$p_val_adj + 10^(-300)
  
  csv_name <- paste0(DEG_out_dir.new,'/',paste0(comparisons,collapse = ' vs '),'.csv')
  write.csv2(mks,csv_name)
  
  volcano(DEG_out_dir,p_hold = 0.01,log2_fc_hold = 0.4)
  
  mks_sig <- filter(mks, abs(avg_log2FC) > 0.3 & p_val_adj < 0.01)
  DEG2RNK(mks_sig,csv_name %>% str_remove('.csv'),log2fc_hold = 0.3 ,p_hold = 0.01)
  
  if(species == 'mouse'){
    GSEA_batch.rnk(
      GSEA_installed_path = GSEA_installed_path,
      DEG_path = DEG_out_dir.new,
      species = 'mouse',
      gene_sets = c(`hallmark gene sets`='mh.all.v2025.1.Mm.symbols.gmt',
                    `GOBP gene sets` = 'm5.go.bp.v2025.1.Mm.symbols.gmt',
                    `Reactome gene sets` = 'm2.cp.reactome.v2025.1.Mm.symbols.gmt'),
      symbol_chip = 'Mouse_Gene_Symbol_Remapping_MSigDB.v2025.1.Mm.chip',
      out_dir = GSEA_out_dir,
      GSEA_plots_number = 30,
      collapse = 'Collapse') 
  }else{
    GSEA_batch.rnk(
      GSEA_installed_path = GSEA_installed_path,
      DEG_path=DEG_out_dir.new,
      species = 'human',
      gene_sets =c(`hallmark gene sets`='h.all.v2025.1.Hs.symbols.gmt',
                   `Reactome` = 'c2.cp.reactome.v2025.1.Hs.symbols.gmt',
                   KEGG = 'c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt'),
      symbol_chip='Human_Gene_Symbol_with_Remapping_MSigDB.v2025.1.Hs.chip',
      out_dir=GSEA_out_dir,
      GSEA_plots_number=30,
      collapse='Collapse'
    )}
GSEA_bubble_3(GSEA_out_dir)}





GSEA_bubble_2 <- function(GSEA_folder='./Tumor cell/GSEA',
                          p.value = c('NOM.p.val','FDR.q.val')[1],
                        height_factor=1,
                        width_factor=1,
                        GSEA_fdr_hold=0.1,
                        topn=20,
                        pos_group = 1,
                        min.fdr.display=10**-10){ ## if want to plot the comparison reversely specify pos_group = 2
  
  
  
  files <- list.files(GSEA_folder,pattern='gsea_report_for',recursive = T,full.names = T)
  files <- data.frame(file_name=grep('.tsv',files,value = T))
  files$ident <- lapply(files$file_name,function(x) str_sub(x,-17,-1)) %>% unlist() %>% invisible()
  save_dirs <- files$file_name %>% dirname() %>% dirname()
  
  batch <- function(iden=files$ident[3],fl = files$file_name[3]){
    print(paste0('Launch GSEA plot for ' ,dirname(files$file_name[files$ident==iden]) %>% basename() %>% unique()) )
    data <- lapply(files$file_name[files$ident==iden],function(x) read.table(x,sep="\t",quot="",header = T)%>%invisible())
    
    data[[1]]$p_val <- data[[1]][,p.value]
    data[[2]]$p_val <- data[[2]][,p.value]
    
    data <- lapply(data, function(x){
      x <- x %>% filter(p_val < GSEA_fdr_hold) 
      x <- x[order(x$p_val),]
      
    } )
    data <- rbind(data[[1]] %>% top_n(topn,-data[[1]]$p_val),data[[2]] %>% top_n(topn,-data[[2]]$p_val)) 
    
    if(pos_group == 2){
      data$NES <- -data$NES
    }
    
    plot <- function(dt=data){
      if(nrow(dt)==0){
        return('Skip empty data')
      }
      
      dt$NAME <- gsub('_',' ',dt$NAME)
      dt$NAME <- lapply(dt$NAME,function(x) str_extract(x," .*")%>% tolower()) %>% unlist() 
      
      dt$NAME <- str_wrap(dt$NAME, width = 40,  indent = 2,whitespace_only = T)
      #dt$NES <- abs(dt$NES)
      dt$group <- ifelse(dt$NES > 0 ,'Up','Down')
      dt$NAME <- factor(dt$NAME,levels = dt$NAME[order(dt$NES,decreasing = F)])
      dt$`Set size` <- dt$SIZE
      
      dt$p_val <- ifelse(dt$p_val == 0, dt$p_val + min(min.fdr.display,dt$p_val[dt$p_val!=0]),
                               dt$p_val)
      
      p<-ggplot()+
        geom_point(data = dt[dt$group == 'Down',],aes(x=abs(NES)*1.2,y= NAME,size=`Set size`,color=-log10(p_val) ))+
        geom_col(data = dt[dt$group == 'Down',],aes(x=abs(NES),y= NAME,fill=-log10(p_val)),width = 0.45)+
        
        scale_color_gradientn(colors = c('#c8d7f8','#3270fa'),name =  paste0('-log10.',p.value,'.Down'),
                              limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))),
                              breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))))+
        scale_fill_gradientn(colors = c('#c8d7f8','#3270fa'),name =  paste0('-log10.',p.value,'.Down'),
                             limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))),
                             breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))))+
        
        new_scale_fill()+
        new_scale_color()+
        
        geom_point(data = dt[dt$group == 'Up',],aes(x=abs(NES)*1.2,y= NAME,size=`Set size`,color=-log10(p_val)))+
        geom_col(data = dt[dt$group == 'Up',],aes(x=abs(NES),y= NAME,fill=-log10(p_val)),width = 0.45)+
        scale_color_gradientn(colors = c('#efb8b8','#f7373a'),name =  paste0('-log10.',p.value,'.Up'),
                              limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))),
                              breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))))+
        scale_fill_gradientn(colors = c('#efb8b8','#f7373a'),name = paste0('-log10.',p.value,'.Up'),
                             limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))),
                             breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))))+

        theme_bw()+
        xlim(-0.05,max(abs(dt$NES)) * 1.4)+
        geom_hline(yintercept = nrow(dt[dt$group == 'Down',]) + 0.5,linetype='dashed',)+
        ylab(NULL)+
        xlab('Absolute NES')+
        ggtitle(basename(strsplit(dirname(fl),'.Gsea')[[1]][1]))+
        theme(axis.text.y=element_text(size=10),
              legend.key.size = unit(0.5, "cm"),
              legend.key.height =  unit(0.5, "cm"),
              legend.key.width =  unit(0.5, "cm"))
      
      png(paste0(strsplit(dirname(fl),'.Gsea')[[1]][1],'.png'),height = max(5.6,0.3*nrow(dt))*height_factor,width = 5.5*width_factor,units = 'in',res=800)
      
      print(p)
      
      dev.off()
      
      

      #ggsave(paste0(strsplit(dirname(fl),'.Gsea')[[1]][1],'.png'),height = max(7,0.3*nrow(dt))*height_factor,width = 4*width_factor)
    }
    
    plot(data)
    
  }
  lapply(1:nrow(files),function(x) try(batch(files$ident[x],files$file_name[x])))
  return('GSEA plot finished ....')
}
#############################

########################################

### similar as GSEA_bubble_2 but no bars included

GSEA_bubble_3 <- function(GSEA_folder='./Tumor cell/GSEA',
                          p.value = c('NOM.p.val','FDR.q.val')[1],
                          height_factor=1,
                          width_factor=1,
                          GSEA_fdr_hold=0.1,
                          topn=10,
                          pos_group = 1,
                          min.fdr.display='auto'){ ## if want to plot the comparison reversely specify pos_group = 2
  
  
  
  files <- list.files(GSEA_folder,pattern='gsea_report_for',recursive = T,full.names = T)
  files <- data.frame(file_name=grep('.tsv',files,value = T))
  files$ident <- lapply(files$file_name,function(x) str_sub(x,-17,-1)) %>% unlist() %>% invisible()
  save_dirs <- files$file_name %>% dirname() %>% dirname()
  
  batch <- function(iden=files$ident[3]){
    print(paste0('Launch GSEA plot for ' ,dirname(files$file_name[files$ident==iden]) %>% basename() %>% unique()) )
    data <- lapply(files$file_name[files$ident==iden],function(x) read.table(x,sep="\t",quot="",header = T)%>%invisible())
    fl <- files$file_name[files$ident==iden][1]
    data[[1]]$p_val <- data[[1]][,p.value]
    data[[2]]$p_val <- data[[2]][,p.value]
    
    data <- lapply(data, function(x){
      x <- x %>% filter(p_val < GSEA_fdr_hold) 
      x <- x[order(x$p_val),]
      
    } )

    dt <- rbind( top_n(data[[1]],topn,-p_val), top_n(data[[2]],topn,-p_val))  
    
    if(pos_group == 2){
      dt$NES <- -dt$NES
    }
    
    plot <- function(dt){
      if(nrow(dt)==0){
        return('Skip empty data')
      }
      
      dt$NAME <- gsub('_',' ',dt$NAME)
      dt$NAME <- lapply(dt$NAME,function(x) str_extract(x," .*")%>% tolower()) %>% unlist() 
      
      dt$NAME <- str_wrap(dt$NAME, width = 40,  indent = 2,whitespace_only = T)
      #dt$NES <- abs(dt$NES)
      dt$group <- ifelse(dt$NES > 0 ,'Up','Down')
      dt$NAME <- factor(dt$NAME,levels = dt$NAME[order(dt$NES,decreasing = F)])
      dt$`Set size` <- dt$SIZE
      
      if(min.fdr.display == 'auto'){
        Up <- dt$p_val[which(dt$group == 'Up')]
        if(length(Up) > 0){
          if(length (which(Up != 0)) == 0){
            dt$p_val[which(dt$group == 'Up')] <- 10 ** -10
          }else{
            dt$p_val[which(dt$group == 'Up')] <- ifelse(Up  == 0, 
                                                        Up  + 
                                                          min(Up[Up !=0]/10) %>% 
                                                          max(10**-10),
                                                        Up )
          }}
        Down <- dt$p_val[which(dt$group == 'Down')]
        if(length(Down) > 0){
          if(length (which(Down != 0)) == 0){
            dt$p_val[which(dt$group == 'Down')] <- 10 ** -10
          }else{
            dt$p_val[which(dt$group == 'Down')] <- ifelse(Down  == 0, 
                                                          Down  + 
                                                          min(Down[Down !=0]/10) %>% 
                                                          max(10**-10),
                                                          Down )
          }}
      }else{
        dt$p_val <- ifelse(dt$p_val == 0, dt$p_val + min(min.fdr.display,dt$p_val[dt$p_val!=0]),
                           dt$p_val)
      }
      
      
      p<-ggplot()+
        geom_point(data = dt[dt$group == 'Down',],aes(x=abs(NES),y= NAME,size=`Set size`,color=-log10(p_val) ))+
       # geom_col(data = dt[dt$group == 'Down',],aes(x=abs(NES),y= NAME,fill=-log10(p_val)),width = 0.45)+
        
        scale_color_gradientn(colors = c('#c8d7f8','#3270fa'),name =  paste0('-log10.',p.value,'.Down'),
                              limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))),
                              breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))))+
        scale_fill_gradientn(colors = c('#c8d7f8','#3270fa'),name =  paste0('-log10.',p.value,'.Down'),
                             limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))),
                             breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Down','p_val']))))+
        
        new_scale_fill()+
        new_scale_color()+
        
        geom_point(data = dt[dt$group == 'Up',],aes(x=abs(NES),y= NAME,size=`Set size`,color=-log10(p_val)))+
        #geom_col(data = dt[dt$group == 'Up',],aes(x=abs(NES),y= NAME,fill=-log10(p_val)),width = 0.45)+
        scale_color_gradientn(colors = c('#efb8b8','#f7373a'),name =  paste0('-log10.',p.value,'.Up'),
                              limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))),
                              breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))))+
        scale_fill_gradientn(colors = c('#efb8b8','#f7373a'),name = paste0('-log10.',p.value,'.Up'),
                             limits =  c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))),
                             breaks = c(-log10(GSEA_fdr_hold),-log10(min(dt[dt$group == 'Up','p_val']))))+
        
        theme_bw()+
        xlim(min(abs(dt$NES))*0.9,max(abs(dt$NES)) * 1.1)+
        geom_hline(yintercept = nrow(dt[dt$group == 'Down',]) + 0.5,linetype='dashed',)+
        ylab(NULL)+
        xlab('Absolute NES')+
        ggtitle(basename(strsplit(dirname(fl),'.Gsea')[[1]][1]))+
        theme(axis.text.y=element_text(size=10),
              legend.key.size = unit(0.5, "cm"),
              legend.key.height =  unit(0.5, "cm"),
              legend.key.width =  unit(0.5, "cm"))
      
      png(paste0(strsplit(dirname(fl),'.Gsea')[[1]][1],'.png'),height = max(5.6,0.3*nrow(dt))*height_factor,width = 5.5*width_factor,units = 'in',res=800)
      
      print(p)
      
      dev.off()
      
      
      
      #ggsave(paste0(strsplit(dirname(fl),'.Gsea')[[1]][1],'.png'),height = max(7,0.3*nrow(dt))*height_factor,width = 4*width_factor)
    }
    
    plot(dt)
    
  }
  lapply(unique(files$ident),function(x) try(batch(x)))
  return('GSEA plot finished ....')
}



# in development
# enrichR_bubble <- function(enrich_result){
#   enrich_result$Gene_ratio <- lapply(enrich_result$Overlap,
#                                      function(x) (as.numeric(strsplit(x,'/')[[1]][1])) / 
#                                        as.numeric(strsplit(x,'/')[[1]][2])) %>% unlist()
#   
#   ggplot(enrich_result,aes(x=Gene_ratio,y=-log10(Adjusted.P.value)))+
#     geom_point()
#   
# }

############################### Spatial toolkit

## RUN RCTD
## Calculate co-Localization
## Calculate infiltration 


















