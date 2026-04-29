library(tidyverse)
library(readxl)
library(Seurat)
library(spacexr)
#library(SingleR)
library(scales)
options(timeout = 60000)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 48000*1024^2)
options(max.print = 10)

lapply(c(list.files("../R tools",full.names = T,pattern = '.R'),
         list.files("../R tools",full.names = T,pattern = '.r')),source) %>% invisible()

obj1 <- readRDS('seuratObject_241127_Yumeng.1.RDS')
obj2 <- readRDS('seuratObject_241127_Yumeng.2.RDS')
obj <- merge(obj1,obj2,merge.data = T,merge.dr = T)

obj <- obj%>% subset(qcCellsFlagged == F &
                       qcFlagsFOV == 'Pass' &
                       qcFlagsCellPropNeg == 'Pass' &
                       qcFlagsCellComplex == 'Pass' & 
                       qcFlagsCellArea == 'Pass' &
                       qcFlagsCellCounts == 'Pass')

obj$group <- 'None'

obj$group[which(obj$fov %in% c(1:66) & obj$slide_ID_numeric == 1)] <- 'N8KO_aPD1'
obj$group[which(obj$fov %in% c(67:161) & obj$slide_ID_numeric == 1)] <- 'WT_IgG'

obj$group[which(obj$fov %in% c(1:67) & obj$slide_ID_numeric == 2)] <- 'N8KO_IgG'
obj$group[which(obj$fov %in% c(68:169) & obj$slide_ID_numeric == 2)] <- 'WT_aPD1'

write_rds(obj,"seuratObject_241127_Yumeng.all.filtered.RDS")


#### assign general cell identity 
##### RCTD
## reference
### both the scRNAseq and Spatial should be un-transformed raw counts

obj <- read_rds("seuratObject_241127_Yumeng.all.filtered.RDS")

### using annotated mouse scRNAseq as reference
ref.obj <- read_rds('../irineos_mice_scRNA_final_version/cleaned.cells.rds')

sub_meta <- ref.obj@meta.data[,c('celltype','nCount_SCT')] %>% as.data.frame()
colnames(sub_meta)[1] <- 'celltype'

sub_meta$cellname <- sub_meta%>%rownames()
sub_index <-  sub_meta  %>%
  group_by(celltype) %>% slice_sample(n=300,replace=F)

sub.ref.obj <- ref.obj[,sub_index$cellname]

ct <- sub.ref.obj$celltype %>% str_replace('/','_') %>%as.factor()

DefaultAssay(sub.ref.obj) <- 'RNA'
sub.ref.obj$celltype <- ct
sub.ref.obj[['SCT']] <- NULL

ref.counts <- GetAssayData(sub.ref.obj,assay = 'RNA',layer = 'counts')
names(ct) <- colnames(sub.ref.obj)

ref.nUMI <- sub.ref.obj$nCount_RNA
names(ref.nUMI) <- colnames(sub.ref.obj)

ref <- Reference(ref.counts,ct ,ref.nUMI)

## query
query.counts <- GetAssayData(obj,assay = 'RNA',layer = 'counts')
query.coord <- GetTissueCoordinates(obj,image = 'X241127_Yumeng.1') %>% 
  rbind(GetTissueCoordinates(obj,image = 'X241127_Yumeng.2'))

rownames(query.coord) <- query.coord$cell
query.coord$cell <- NULL

query <- SpatialRNA(query.coord,query.counts,colSums(query.counts))

## run RCTD

RCTD <- create.RCTD(query,ref,max_cores = 12)

rm(list = ls()[ls() != 'RCTD'])
gc()
library(spacexr)

RCTD <- run.RCTD(RCTD,doublet_mode = 'doublet')

saveRDS(RCTD,'RCTD.rds')

################ assign cell idents

obj <- read_rds("seuratObject_241127_Yumeng.all.filtered.RDS")
RCTD <- read_rds('RCTD.rds')

re.df <- RCTD@results$results_df
re <- re.df$first_type
names(re) <- rownames(re.df)

levels(re) <- c(levels(re),'Unassigned')

obj$celltype_rctd <- re
obj$celltype_rctd[is.na(obj$celltype_rctd)] <- 'Unassigned'


#obj <- obj[,(is.na(obj$celltype_rctd) == F)]
obj <- UpdateSeuratObject(obj)

write_rds(obj,'obj_RCTD_assigned_cell_type.rds')

########
obj <- read_rds('obj_RCTD_assigned_cell_type.rds')

DotPlot(obj ,c('Ptprc',
               'Ncr1','Nkg7','Klrb1c','Ncam1',
               'Krt8','Krt19','Epcam',
               'Il3ra','Tcf7','Plac8',
               'Pecam1','Cav1','Eng',
               'Col1a1','Col1a2','Col3a1','Col6a1','Col6a2','Fap',
               'Cd80','Cd68','Cd14','Itgax','Itgam','Csf1r',
               'C1qa','C1qb',
               'Prf1','Gzmb','Gzma',
               'Cd8a','Cd8b1','Cd3d','Cd3e',
               'Cd4','Foxp3','Il2ra'),group.by = 'celltype_rctd')+
  scale_colour_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'))

############# plotting 
############# 
obj$celltype_rctd %>% unique()

p1<- ImageDimPlot(obj %>% subset((celltype_rctd %in% c('Macrophage','CD8 T cell','NK cell' ))
                            & (group == 'WT_aPD1')),group.by = 'celltype_rctd',
             ,
             cols = c('#EC6363','#63A5EC','#6BCC5E'))+
  ggtitle('WT_aPD1')
  
p2<-  ImageDimPlot(obj %>% subset((celltype_rctd %in% c('Macrophage','CD8 T cell','NK cell' ))
                              & (group == 'N8KO_aPD1')),group.by = 'celltype_rctd',
               ,coord.fixed = T,
               cols = c('#EC6363','#63A5EC','#6BCC5E'))+
    ggtitle('N8KO_aPD1')
  
p3<-  ImageDimPlot(obj %>% subset((celltype_rctd %in% c('Macrophage','CD8 T cell','NK cell' ))
                              & (group == 'WT_IgG')),group.by = 'celltype_rctd',
               ,coord.fixed = T,
               cols = c('#EC6363','#63A5EC','#6BCC5E'))+
    ggtitle('WT_IgG')
  
p4<-  ImageDimPlot(obj %>% subset((celltype_rctd %in% c('Macrophage','CD8 T cell','NK cell' ))
                              & (group == 'N8KO_IgG')),group.by = 'celltype_rctd',
               ,coord.fixed = T,
               cols = c('#EC6363','#63A5EC','#6BCC5E'))+
    ggtitle('N8KO_IgG')

p <- p1+p2+p3+p4


ggsave('all.svg',p,width=20,height=20)



## cell pct plot Fig S7A
plot.dt <- obj@meta.data[,c('group','celltype_rctd')] %>% filter(celltype_rctd != 'Unassigned') %>% table() %>% prop.table(margin =1) %>% as.data.frame() %>% filter(celltype_rctd != 'Unassigned')

ggplot(plot.dt,aes(x=group,y=Freq,fill = celltype_rctd))+
  geom_col()+
  geom_text(aes(label = ifelse(Freq > 0.02,paste0(round(Freq,2)*100, ' %'),'')),position = position_stack(vjust = 0.5))


###cell-cell co-localization

obj <- readRDS('obj_RCTD_assigned_cell_type.rds')
obj <- obj %>% subset(celltype_rctd != 'Unassigned')

obj$celltype_rctd <- as.character(obj$celltype_rctd)

coord <- GetTissueCoordinates(obj,image = 'X241127_Yumeng.1') %>% 
  rbind(GetTissueCoordinates(obj,image = 'X241127_Yumeng.2'))

coord <- cbind(coord,obj@meta.data[coord$cell,c('celltype_rctd','group','slide_ID_numeric')])

rownames(coord) <- coord$cell

radius <- 0.1
background.size = 2000
test.size = 500


#####
#write_rds(coord,'all.coord.rds')
#rm(obj)

coord <- read_rds('all.coord.rds')

#sample <- coord$group[[1]]

co_local_score <- function(coord,  radius, sample='WT_aPD1'){
  ## subet coord of one sample
  coord.sub <- coord %>% filter(group == sample)
  
  ### select sample number of cell for each cell type as test 
  test_cells <- coord[,c('cell','group','celltype_rctd')] %>% 
    filter(celltype_rctd %in% c('Mast cell') == F ) %>% 
    filter(group == sample) %>% 
    group_by(celltype_rctd) %>% 
    slice_sample(n=test.size,replace = F)
  
  ### random select cells as background
  background_cells <- coord[,c('cell','group')] %>%
    filter(group == sample) %>% 
    slice_sample(n=background.size,replace = F)
  
  ### calculate numbers of each type of coefficient cells
  adjacent_cell_freq <- function(coord.matrix, cells,radius){
    freq.table <- NULL
    
    for(n in cells){  
      cell_a <- coord.matrix %>% filter(cell == n)
      cell_b <- coord.matrix %>% filter(cell != n)
      ## pre_filter
      adjacent_cells <- cell_b[(abs(cell_b$x - cell_a$x) <= radius &
                                  abs(cell_b$y - cell_a$y) <= radius) ,]
      ## circular region
      adjacent_cells <- adjacent_cells[((adjacent_cells$x - cell_a$x)^2 + 
                                          (adjacent_cells$y - cell_a$y)^2 <= radius^2),]
      
      stat <- merge(cell_a,data.frame(table(adjacent_cells$celltype_rctd)),all=T) %>% as.data.frame()
      
      freq.table <- rbind(freq.table,stat)
      
    } 
    return(freq.table)
  }
  ## Mean background as expectation
  background.freq <- adjacent_cell_freq(coord.sub, background_cells$cell,radius)
  background.freq.mean <- background.freq %>% group_by(Var1) %>% mutate(expected.Freq=mean(Freq)) %>%
    dplyr::select(Var1,expected.Freq) %>% unique()
  
  ## test combine with background
  test.freq <- adjacent_cell_freq(coord.sub, test_cells$cell,radius)
  test.freq <- merge(test.freq,background.freq.mean,by='Var1')
  
  ## calculate coefficient for cell-cell paris
  #  coeff <- (observation of adjacent cell number  / expectation of adjacent cell number) -1
  coefficient <- test.freq %>% group_by(cell) %>% mutate(coefficient = (Freq/expected.Freq)-1 ) %>% 
    ## then mean coeff
    
    group_by(celltype_rctd,Var1) %>% mutate(mean.coefficient = mean(coefficient)) %>%
    dplyr::select(Var1,celltype_rctd,mean.coefficient) %>% unique()
  return(coefficient %>% as.data.frame())
  
}

scores <- lapply(coord$group %>% unique(), function(x) co_local_score(coord, radius, x))

names(scores) <- coord$group %>% unique()
#write_rds(scores,'co-score.rds')

scores <- read_rds('co-score.rds')

for(n in names(scores)){
  scores[[n]]$group <- n
}


scores <- do.call('rbind',scores)
dt <- scores %>% filter(Var1 %in% c('Unassigned','Mast cell',
                                    'NK cell','Endothelial cell')==F &
                          celltype_rctd %in% c('Unassigned','Mast cell',
                                               'NK cell','Endothelial cell')==F )

###
dt <- dt %>% group_by(Var1,celltype_rctd) %>% mutate(scaled.coeff = rescale(mean.coefficient,to = c(-1,1) ))

ggplot( dt,
        aes(x=Var1,y=celltype_rctd,fill=mean.coefficient))+
  facet_wrap(~group,scales = 'free')+
  geom_tile(width = 0.95, height = 0.9)+
  theme_bw()+
  scale_fill_gradientn(colors = c('#2760DA','#EEEEEF','#E0191B'),breaks = c(-1,0,1),limits =c(-1.5,1.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=12))+
  xlab('Adjacent cell')+
  ylab('Central cell')

### differential co-local score
dt$pair <- paste0(dt$Var1,dt$celltype_rctd)
dt.list <- dt %>% split.data.frame(dt$group)
dt.diff <- dt %>% distinct(Var1,celltype_rctd,pair)
dt.diff <- dt.diff[order(dt.diff$pair),]

dt.diff$diff.coeff <- dt.list[['N8KO_aPD1']][order(dt.list[['N8KO_aPD1']]$pair),"mean.coefficient"] - 
  dt.list[['WT_aPD1']][order(dt.list[['WT_aPD1']]$pair),"mean.coefficient"]

## check order
which(dt.diff$pair != dt.list[['N8KO_aPD1']][order(dt.list[['N8KO_aPD1']]$pair),'pair'])
which(dt.diff$pair != dt.list[['WT_aPD1']][order(dt.list[['WT_aPD1']]$pair),'pair'])


## Fig 4D

png('co-local-score.group_scaled.png',width=7,height = 6,units = 'in',res=600)

ggplot( dt,
        aes(x=Var1,y=celltype_rctd,fill=scaled.coeff))+
  facet_wrap(~group,scales = 'free')+
  geom_tile(width = 0.95, height = 0.9)+
  theme_bw()+
  scale_fill_gradientn(colors = c('#5dade2','#EEEEEF','#e25d5d'),
                       breaks = c(-1,1),labels =c('Min','Max'),name = 'Coeff.scaled' ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=12))+
  xlab('Adjacent cell')+
  ylab('Central cell')

dev.off()

######################## fluorescence plottings

library(ggspatial)

fluo_plot <- function(coord,cells.plot,pt.size=0.5,cols,alpha=1){
  
  lapply(unique(coord$group),function(z){
    
    coord.sub <- coord %>% subset(group == z & celltype_rctd == cells.plot)
    
    ylim <- coord.sub$y %>% range()
    height <- ylim[2]-ylim[1]
    xlim <- coord.sub$x %>% range()
    width <- xlim[2]- xlim[1]
    
    colors <- as.list(cols)
    names(colors) <- cells.plot
    
    p<-ggplot(coord.sub,aes(x=x,y=y))+
      geom_point(aes(color = celltype_rctd),size=pt.size,alpha=alpha)+
      coord_fixed()+
      theme(legend.position = 'none',
            plot.margin = margin(0,0,0,0,'in'),
            axis.title = NULL,
            axis.text = NULL)+
      scale_color_manual(values = colors)+
      
      geom_segment(
        aes(x = xlim[1]+0.09, xend = xlim[1]+0.61, 
            y = ylim[1]+0.1, yend = ylim[1]+0.1),
        color = '#ffffff',
        linewidth = 0.5
      ) +
      geom_segment(
        aes(x = xlim[1]+0.1, xend = xlim[1]+0.1, 
            y = ylim[1]+0.1, yend = ylim[1]+0.2),  
        color = '#ffffff',
        linewidth = 0.5
      ) +
      geom_segment(
        aes(x = xlim[1]+0.6, xend = xlim[1]+0.6, 
            y = ylim[1]+0.1, yend = ylim[1]+0.2), 
        color = '#ffffff',
        linewidth = 0.5
      ) +
      
      
      annotate(
        "text",
        x = xlim[1]+0.1, y = ylim[1]+0.4,
        label = "0",
        color = '#ffffff',
        size = 4
      ) +
      annotate(
        "text",
        x = xlim[1]+0.6, y = ylim[1]+0.4,
        label = "0.5",
        color = '#ffffff',
        size = 4
      ) +
      annotate(
        "text",
        x = xlim[1]+0.35, y = ylim[1]-0.1,
        label = "Pixels",
        color = '#ffffff',
        size = 4
      ) +
      
      theme_dark()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            plot.margin = margin(0,0,0,0),
            panel.background = element_rect(fill = "#000000"),
            legend.title = element_text(size=10),
            legend.text = element_text(size=10)
      )+
      guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
    
    svg(paste0(paste0(cells.plot,collapse = '_'),' legend.svg',sep = ' '),width = 4,height = 4)
    print(get_legend(p)%>%as_ggplot())
    dev.off()
    
    svg(paste0(z,' -- ',paste0(cells.plot,collapse = '_'),'.svg',sep = ' '),width = width/1.2,height = height/1.2)
    print(p+theme(legend.position = 'none'))
    dev.off()  
    
    
  })
}

# Fig 4A.B

fluo_plot(coord,cells.plot= c('CD8 T cell','CD4 T cell'),pt.size=0.5,cols=c('#FE9900','#55A871'),alpha = 0.8)

fluo_plot(coord,c('Tumor cell','Fibroblast'),pt.size=0.005,cols=c('#789BC7','#E98091'),alpha = 0.5)

fluo_plot(coord,c('CD8 T cell','Macrophage'),pt.size=0.05,cols=c('#FE9900','#62BBCA'),alpha = 0.5)


#################################################################################### infiltration analysis
# library(ggalt)  
library(area)
library(StereoMorph)

coord <- read_rds('all.coord.rds')

calculate_hull_dis <- function(coord,sample){

  smallest_dist <- function(point.x,point.y,hull_lines){
    
    dist <- apply(hull_lines, 1, function(x) {
      
      distancePointToLine(c(point.x,point.y),c(x[1],x[2]),c(x[3],x[4]))
    })
    
    return(min(dist))
  }
  
  hulls <- function(coord,gp,sample=sample){
    coord <- coord %>% filter(group == gp)
    
    if(sample){
      coord <- coord %>% slice_sample(prop=sample)
    }
    
    
    hull_indices <- chull(coord$x, coord$y)
    hull_points <- coord[hull_indices, ]
    #print('done')
    
    png(paste0(gp,'convex.png'),width = 4,height = 4,units = 'in',res=600 ) 
    plot(coord$x, coord$y, xlab = "X", ylab = "Y", main = paste0(gp," Convex Hull Encircle"))
    lines(hull_points$x, hull_points$y, col = "red", lwd = 2) # Draw convex hull
    points(c(1,hull_points$x), c(1,hull_points$y), col = "blue", pch = 16) # Highlight hull points
    dev.off()
    
    #### to hull distance
    
    len <- nrow(hull_points)
    hull_lines <- data.frame(point.1.x = hull_points$x[1:len],
                             point.1.y = hull_points$y[1:len],
                             point.2.x = hull_points$x[c(2:len,1)],
                             point.2.y = hull_points$y[c(2:len,1)])
    
    
    
    coord$min_dis2edge <- apply(coord,1 ,function(x)
    {smallest_dist(as.numeric(x[1]),as.numeric(x[2]),hull_lines)} )
    
    ### 

    avg_radius <- (polygon_area(hull_points[,c('x','y')]) / pi) ^ 0.5
    
    coord$`min_dis2edge/radius` <- coord$min_dis2edge/avg_radius
    
    
    coord$tumor_radius <- avg_radius
    return(coord)
  }
  
  return(lapply(unique(coord$group),function(x) hulls(coord,x,sample)))
  
}

dis <- do.call('rbind',calculate_hull_dis(coord,sample = F))

#write.csv(dis,'distance.csv')

dis <- read.csv('distance.csv')

library(ggridges)

png('infiltration.png',width=8,height = 5,units = 'in',res=600)

dt <- dis %>% group_by(group) %>% slice_sample(n=40000) %>% ungroup() %>% 
  filter(celltype_rctd %in% c('Mast cell','NK cell','Endothelial cell','Tumor cell')==F)


dt <- dt %>% group_by(group,celltype_rctd) %>% mutate(`mean.min_dis2edge.radius` = mean(`min_dis2edge.radius`))

### Fig 4C
ggplot(dt ,aes(x=`min_dis2edge.radius`, y=group,fill = celltype_rctd ))+
  facet_wrap(~celltype_rctd )+
  #geom_density_ridges2(colour = 'grey',alpha=0.65)+
  theme_bw()+
  geom_density_ridges2( quantile_lines=TRUE, 
                        quantile_fun=function(`min_dis2edge.radius`,...)mean(`min_dis2edge.radius`),
                        alpha=0.5,linetype='dashed',color=alpha('black',0.3),scale=1.3,
                        rel_min_height = 0.001
  )+
  geom_text(aes(x= 0.7,label=round(mean.min_dis2edge.radius,2)))




#### ############### Macrohage niche definition

coord.sub <- coord %>% subset(celltype_rctd == 'Macrophage')


### SNN clustering
library(dbscan)


coord.sub$cluster <- NA

for(g in coord.sub$group %>% unique()){
  
  cl <- dbscan(coord.sub%>% filter(group==g) %>% dplyr::select(x,y), eps = 0.03, minPts = 3)
  
  #### only keep cluster with cell number greater than x
  cl$cluster[cl$cluster %in% (table(cl$cluster)[table(cl$cluster) < 25] %>% names())] <- '0'
  
  coord.sub$cluster[which(coord.sub$group==g)] <- cl$cluster %>% as.character()
  
  
}
rm(cl)


coord.sub$cell_idents <- paste(coord.sub$group,coord.sub$cluster)


#write_rds(coord.sub,'macro.coord.niche.assigned.rds')


coord.sub <- read_rds('macro.coord.niche.assigned.rds')

ggplot(coord.sub %>% filter(cluster!=0),aes(x=x,y=y,colour = cell_idents))+
  
  facet_wrap(~group,scales = 'free')+
  #geom_point(data =  coord %>% subset(group=='N8KO_aPD1'),color='white')+
  geom_point(data =  coord.sub %>% filter(cluster==0),color='grey',size=0.5)+
  geom_point(size=0.5)+
  scale_color_discrete()+
  theme(legend.position = 'none')


ggplot(coord.sub %>% filter(cluster!=0 & group =='N8KO_aPD1'),aes(x=x,y=y,colour = cell_idents))+
  
  facet_wrap(~group,scales = 'free')+
  #geom_point(data =  coord %>% subset(group=='N8KO_aPD1'),color='white')+
  geom_point(data =  coord.sub %>% filter(cluster==0& group =='N8KO_aPD1'),color='grey',size=0.5)+
  geom_point(size=0.5)+
  scale_color_discrete()+
  theme(legend.position = 'none')+
  
  
  ggplot(coord.sub %>% filter(cluster!=0 & group =='WT_aPD1'),aes(x=x,y=y,colour = cell_idents))+
  
  facet_wrap(~group,scales = 'free')+
  #geom_point(data =  coord %>% subset(group=='N8KO_aPD1'),color='white')+
  geom_point(data =  coord.sub %>% filter(cluster==0& group =='WT_aPD1'),color='grey',size=0.5)+
  geom_point(size=0.5)+
  scale_color_discrete()+
  theme(legend.position = 'none')

rm(obj)


### add niche region number to obj

obj$niche <- NA

obj$niche[coord.sub$cell] <- coord.sub$cell_idents


ImageDimPlot(obj %>% subset(celltype_rctd == 'Macrophage') ,
             fov='X241127_Yumeng.1',group.by = 'niche',na.value ='grey')+
  theme(legend.position = 'none')


#### add niche info to neibor cells

cells <- coord.sub$cell[coord.sub$cluster!=0]

neighbor.rad <- 0.05

coord.sub.ex <- coord %>% subset(celltype_rctd != 'Macrophage')

coord.sub.ex$cluster <- NA

coord.neighbor <- data.frame()

for(cell in  cells){
  
  pt <- coord.sub[cell,c('x','y','cluster','group')]
  
  cell.neighbor <- coord.sub.ex[ ((coord.sub.ex$x - pt$x)^2 < neighbor.rad^2) &
                                   ((coord.sub.ex$y - pt$y)^2 < neighbor.rad^2) & 
                                   coord.sub.ex$group == pt$group,]
  
  if(nrow(cell.neighbor) > 0 ){
    cell.neighbor$host <- cell
    cell.neighbor$cluster <- pt$cluster
    if(nrow(coord.neighbor)==0){
      coord.neighbor <- cell.neighbor
    }else{
      coord.neighbor <- rbind(coord.neighbor,cell.neighbor)
    }
  }
}


write_rds(coord.neighbor,'coord.neighbor.all.rds')

coord.neighbor$cell_idents <- paste(coord.neighbor$group,coord.neighbor$cluster)

## unique the cells for each niche

coord.neighbor.new <- data.frame()

for(niche in coord.neighbor$cell_idents %>% unique()){
  index <- which(coord.neighbor$cell_idents == niche)
  coord.niche <- coord.neighbor[index,]
  
  coord.neighbor.new <- rbind(coord.neighbor.new,
                              distinct(coord.niche,cell,.keep_all = T))
  
}

write_rds(coord.neighbor.new,'Macro_neighbors.rds')



##################
#coord.neighbor<- readRDS('Macro_neighbors.rds')


test <- coord.neighbor.new[unique(coord.neighbor.new$cell),]
ggplot(test %>% filter(cluster!=0 & group=='WT_aPD1'),aes(x=x,y=y,colour = cell_idents))+
  #facet_wrap(~cluster)+
  #geom_point(data =  coord %>% subset(group=='WT_aPD1'),color='white')+
  geom_point(data =  coord.sub %>% filter(group=='WT_aPD1' & cluster==0),color='grey')+
  geom_point(size=0.1)+
  scale_color_discrete()+
  theme(legend.position = 'none')


test <- coord.neighbor.new[unique(coord.neighbor.new$cell),]
ggplot(test %>% filter(cluster!=0 & group=='N8KO_aPD1'),aes(x=x,y=y,colour = cell_idents))+
  #facet_wrap(~cluster)+
  #geom_point(data =  coord %>% subset(group=='WT_aPD1'),color='white')+
  geom_point(data =  coord.sub %>% filter(group=='N8KO_aPD1' & cluster==0),color='grey')+
  geom_point(size=0.1)+
  scale_color_discrete()+
  theme(legend.position = 'none')


## extract tumor subsets based on niche

obj <- obj %>% SCTransform()

obj.list <- list()

for(nich in coord.neighbor.new$cell_idents %>% unique()){
  sample.index <- c(coord.neighbor.new$host[coord.neighbor.new$cell_idents==nich] %>% unique(),
                    coord.neighbor.new$cell[coord.neighbor.new$cell_idents==nich])
  
  
  obj.list <- append(obj.list,obj@assays$SCT$data[,sample.index])
}

names(obj.list) <- coord.neighbor.new$cell_idents %>% unique()

rm(obj)

write_rds(obj.list,'macro_niches.obj.list.rds')

### run cell chat

obj.list <- read_rds('macro_niches.obj.list.rds')
coord <- read_rds('all.coord.rds')


meta <- data.frame(label=NULL)


cc.list <- list()

obj.list <- obj.list[c(grep('N8KO_aPD1',names(obj.list),value = T),grep('WT_aPD1',names(obj.list),value = T)) ]


### in general calculation

for(niche in names(obj.list)){
  colnames(obj.list[[niche]])  <- paste0(niche,'+',obj.list[[niche]] %>% colnames)
}

obj.list.integrated <- do.call('cbind',obj.list)


labels <- coord[obj.list.integrated %>% colnames() %>% str_remove('.+\\+'),'celltype_rctd']
samples <- obj.list.integrated %>% colnames() %>% str_remove('\\+.+')
groups <- obj.list.integrated %>% colnames() %>% str_remove('[:space:][:digit:]+\\+.+')

meta <- data.frame(row.names = obj.list.integrated %>% colnames(),
                   labels= labels,
                   samples = samples,
                   groups = groups)


future::plan("multisession", workers = 8)

for( group in meta$groups %>%unique()){
  
  index <- meta %>% filter(groups==group) %>% rownames()
  
  cc <- Run_cellchat(manual_dt = obj.list.integrated[,index] ,
                     manual_meta = meta[index,],
                     ccDB = CellChatDB.mouse,
                     search = c('Cell-Cell Contact','Secreted Signaling'),
                     smooth = PPI.mouse)
  
  cc.list <- append(cc.list,cc)
  
  
}
names(cc.list) <- meta$groups %>%unique()

write_rds(cc.list,'Macro_spatial_cc_group_integrated.rds')


############################

### per niche level calculation

for(niche in names(obj.list)){
  
  labels <- coord[obj.list[[niche]] %>% colnames(),'celltype_rctd']
  
  colnames(obj.list[[niche]])  <- paste0(niche,obj.list[[niche]] %>% colnames)
  
  meta <- data.frame(row.names = obj.list[[niche]] %>% colnames(),
                     labels= labels,
                     samples = rep(niche,ncol(obj.list[[niche]])))
  
  
  future::plan("multisession", workers = 8)
  
  cc <- Run_cellchat(manual_dt = obj.list[[niche]] ,
                     manual_meta = meta,
                     ccDB = CellChatDB.mouse,
                     search = c('Cell-Cell Contact','Secreted Signaling'),
                     smooth = PPI.mouse)
  
  
  cc.list <- append(cc.list,cc)
}

names(cc.list) <- names(obj.list)
write_rds(cc.list,'spatial-cell-chat.Macro_niche.rds')


######## analysis

cc.list <- read_rds('spatial-cell-chat.Macro_niche.rds')
cc.list <- sapply(cc.list,function(x) subsetCellChat(x,idents.use = c('CD8 T cell',
                                                                      'CD4 T cell',
                                                                      'CD4 Treg cell',
                                                                      'Macrophage',
                                                                      'DCs')))

dt <- lapply(names(cc.list), function(x){
  sub = cc.list[[x]]@net$prob %>% melt() %>% dcast(Var3 ~ Var1+Var2)
  rownames(sub) <- sub$Var3
  sub <- sub[,-1]
  colnames(sub) <- paste0(colnames(sub),' + ',x) 
  sub$rn <- rownames(sub)
  return(sub)
}  )

names(dt) <- names(cc.list)


cc_matrix <- Reduce(function(x, y) merge(x, y, by='rn',all=T), dt)

cc_matrix[is.na(cc_matrix)] <- 0

rownames(cc_matrix) <- cc_matrix$rn
cc_matrix <- cc_matrix[,-1]


meta <- data.frame(row.names = colnames(cc_matrix),
                   niche = colnames(cc_matrix) %>% str_remove('.+\\+ '),
                   pairs = colnames(cc_matrix) %>% str_remove(' \\+.+'),
                   group = colnames(cc_matrix) %>% str_remove('.+\\+ ') 
                   %>% str_remove(' [:digit:]+') 
)

niche.cc.obj <- CreateSeuratObject(counts = cc_matrix %>%as.matrix() ,
                                   data = cc_matrix %>%as.matrix(),
                                   meta.data = meta,min.features = 5,min.cells = 2)

write_rds(niche.cc.obj,'niche.cc.obj.rds')

niche.cc.obj <- ScaleData(niche.cc.obj)
niche.cc.obj <- RunPCA(niche.cc.obj, npcs = 50, features = rownames(niche.cc.obj))
ElbowPlot(niche.cc.obj)
niche.cc.obj <- FindNeighbors(niche.cc.obj, reduction = "pca", dims = 1:15)
niche.cc.obj <- FindClusters(niche.cc.obj, resolution = 0.3)

niche.cc.obj <- RunUMAP(niche.cc.obj, dims = 1:15,seed.use = 123)

p1 <-DimPlot(niche.cc.obj,group.by = 'group')
p2 <-DimPlot(niche.cc.obj,group.by = 'pairs')+ theme(legend.position = 'none')
p3 <-DimPlot(niche.cc.obj,label = T)
p4 <- DimPlot(niche.cc.obj,group.by = 'niche') + theme(legend.position = 'none')

p1+p2+p3+p4

features <-  Features(niche.cc.obj)[intersect(grep('H2-',Features(niche.cc.obj),invert = T),
                                              grep('MIF',Features(niche.cc.obj),invert = T))]


test.s <- sapply(niche.cc.obj$pairs %>% unique(), function(x) {
  
  dt <- niche.cc.obj %>% subset(pairs == x)
  Idents(dt) <- 'group'
  mks <- FindAllMarkers(dt,features = features ,
                        min.pct = 0.1)
  #if(nrow(mks)>0){
  #mks <- mks %>% filter(p_val_adj <0.01 & pct.1>0.1 & avg_log2FC>0)
  return(mks)
  #}
})

k<-test.s[[2]]
l<-do.call('rbind',test.s)

write_csv2(l,'cc.csv')

plot.cc <- mks %>% group_by(cluster) %>% top_n(1,avg_log2FC)
plot.matrix <- cc_matrix[plot.cc$gene %>% str_replace_all('-','_'),]


plot.matrix <- apply(plot.matrix,1,function(x) rescale(x,to=c(0,1))) %>% t() %>% melt()

plot.matrix$pair <- plot.matrix$Var2 %>% str_remove(' \\+.+')
plot.matrix$group <- plot.matrix$Var2 %>% str_remove('.+\\+ ') %>% str_remove(' [:digit:]+')
plot.matrix$niche <- plot.matrix$Var2 %>% str_remove('.+\\+ ')

ggplot(plot.matrix,aes(x=niche,y=Var1,fill = value))+geom_tile()



niche.cc.obj$ident <- paste(niche.cc.obj$group,niche.cc.obj$pairs)
Idents(niche.cc.obj) = 'ident'
mks <- FindAllMarkers(niche.cc.obj,features = Features(niche.cc.obj)[grep('H2-',Features(niche.cc.obj),invert = T)],
                      min.pct = 0.1)

mks <- mks %>% filter(p_val_adj <0.01&pct.1>0.1&avg_log2FC>0)

write_csv2(mks,'cc.mks.csv')

table(mks$cluster)

plot.mks <- rbind( mks %>% group_by(cluster) %>% top_n(5,avg_log2FC),
                   mks %>% group_by(cluster) %>% top_n(5,-avg_log2FC)) %>% unique()

write_csv2(mks,'cc.mks.uni.csv')


unique(plot.mks$gene)

## cluster niche pct
pct <- niche.cc.obj@meta.data[,c('niche','seurat_clusters')] %>% table() %>% prop.table(margin = 2) %>% melt()

ggplot(pct,aes(x=seurat_clusters,y=value,fill=ifelse(value>0.05, niche,NA)))+
  geom_bar(stat = 'identity')+ theme(legend.position = 'none')

#########

k<-niche.cc.obj@meta.data %>%filter(seurat_clusters==8) 
k<-niche.cc.obj@meta.data %>%filter(seurat_clusters==5) 

options(max.print = 30)
c<- k[,c('niche')] %>% table()

########################################

result <- read_rds('Macro_spatial_cc_group_integrated.rds')
result <- sapply(result,function(x) subsetCellChat(x,idents.use = c('CD8 T cell',
                                                                    'CD4 T cell',
                                                                    'CD4 Treg cell',
                                                                    'Macrophage',
                                                                    'Tumor cell',
                                                                    'Fibtoblast',
                                                                    'DCs',
                                                                    'NK cell')))
result <- sapply(result, netAnalysis_computeCentrality)

cellchat <- mergeCellChat(result,add.names =names(result))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


gg2 <- netVisual_heatmap(cellchat,
                         comparison = c('WT_aPD1','N8KO_aPD1'),
                         measure = c('weight'),
                         slot.name = 'net',
                         title.name = 'WT αPD1_vs_N8KO αPD1_num_of_diff_interactions')

#> Do heatmap based on a merged object
p <- (gg1 + gg2) + plot_layout(ncol = 2)

ggsave('heatmap_diff_interactions.png',p)

num.link <- sapply(result, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(result)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(result[[i]], 
                                               title = names(result)[i], 
                                               weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg) & xlim(0,3) & ylim(0,3) & 
  geom_vline(xintercept = 1.5,linetype='dashed',alpha=0.4) &
  geom_hline(yintercept = 1.5,linetype='dashed',alpha=0.4)


netAnalysis_signalingChanges_scatter(cellchat, 
                                     idents.use = "Macrophage", 
                                     comparison = c(1,2))

netAnalysis_signalingChanges_scatter(cellchat, 
                                     idents.use = "CD4 T cell", 
                                     comparison = c(1,2))

netAnalysis_signalingChanges_scatter(cellchat, 
                                     idents.use = "CD8 T cell", 
                                     comparison = c(1,2))

## niche spatial cell chat

netVisual_bubble(cellchat, sources.use = 'CD4 T cell', 
                 targets.use = c('Macrophage'),  
                 comparison = c(1,2),
                 max.dataset = 1,
                 thresh = 1,
                 angle.x = 45,
                 remove.isolate = T,return.data = T)


netVisual_bubble(cellchat, 
                 sources.use = 2, 
                 targets.use = c(1,3),  
                 comparison = c(1,2), 
                 max.dataset = 1, 
                 title.name = paste0("Increased signaling in ",names(result)[1]), 
                 angle.x = 45, 
                 remove.isolate = T)


