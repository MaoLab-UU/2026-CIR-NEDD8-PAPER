

################# human analysis template
# GSEA_batch.rnk(GSEA_installed_path = "E:/bioinfo/GSEA_4.3.2",
#                DEG_path="final plots/ranks_diff_2000 result",
#                species = 'human',
#                gene_sets =c(`hallmark gene sets`='h.all.v2025.1.Hs.symbols.gmt',
#                             `Reactome` = 'c2.cp.reactome.v2025.1.Hs.symbols.gmt',
#                             KEGG = 'c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt'),
#                symbol_chip='Human_Gene_Symbol_with_Remapping_MSigDB.v2025.1.Hs.chip',
#                out_dir='final plots',
#                GSEA_plots_number=30,
#                collapse='Collapse')




############## mouse analysis template
# GSEA_batch.rnk(
#   GSEA_installed_path = "E:/bioinfo/GSEA_4.3.2",
#   DEG_path = "final plots/ranks_diff_2000 result",
#   species = 'mouse',
#   gene_sets = c(`hallmark gene sets`='mh.all.v2025.1.Mm.symbols.gmt',
#                 `GOBP gene sets` = 'm5.go.bp.v2025.1.Mm.symbols.gmt',
#                 `Reactome gene sets` = 'm2.cp.reactome.v2025.1.Mm.symbols.gmt'),
#   symbol_chip = 'Mouse_Gene_Symbol_Remapping_MSigDB.v2025.1.Mm.chip',
#   out_dir = 'final plots',
#   GSEA_plots_number = 30,
#   collapse = 'Collapse'
# )






















## reset system PATH
##  usethis::edit_r_profile()


working_dir <- 'cd _cd_'


command.rnk <- ('gsea-cli.bat GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/msigdb/_species_/gene_sets/_gene_set_ -collapse _collapse_ -mode Abs_max_of_probes -norm meandiv -nperm 1000 -rnd_seed timestamp -rnk _rnk_file_ -scoring_scheme weighted -rpt_label _the_label_ -chip ftp.broadinstitute.org://pub/gsea/msigdb/_species_/annotations/_chip_ -create_svgs false -include_only_symbols true -make_sets true -plot_top_x _numberplot_ -set_max 500 -set_min 15 -zip_report false -out _out_dir_')


command.counts <- ('gsea-cli.bat GSEA -res _count_file_ -cls _comparison_ -gmx ftp.broadinstitute.org://pub/gsea/msigdb/_species_/gene_sets/_gene_set_ -collapse _collapse_ -norm meandiv -nperm 1000 -permute phenotype -rnd_seed timestamp -rnd_type no_balance -scoring_scheme weighted -rpt_label _the_label_ -metric Signal2Noise -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/msigdb/_species_/annotations/_chip_ -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x _numberplot_ -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out _out_dir_' )




############################################   GSEA batch pipeline # reauires GSEA command line package

path_syntax_change <- function(path_in_R){
  path_in_windows <- gsub('/','\\\\',path_in_R,fixed = TRUE)
  path_in_windows <- paste0('\\"',path_in_windows,'\\"',sep='')
  return(path_in_windows)
}








GSEA_batch.rnk <- function(
    GSEA_installed_path = "E:/bioinfo/GSEA_4.3.2",
    DEG_path='DEG',
    species = 'mouse',
    gene_sets =c(`hallmark gene sets`='mh.all.v2023.1.Mm.symbols.gmt',
                 `positional gene sets`='m1.all.v2023.1.Mm.symbols.gmt',
                 `curated gene sets`='m2.all.v2023.1.Mm.symbols.gmt',
                 `regulatory target gene sets`='m3.all.v2023.1.Mm.symbols.gmt',
                 `ontology gene sets`='m5.all.v2023.1.Mm.symbols.gmt',
                 `cell type signature gene sets`='m8.all.v2023.1.Mm.symbols.gmt'),
    symbol_chip='Mouse_Gene_Symbol_Remapping_MSigDB.v2023.1.Mm.chip',
    out_dir='GSEA',
    GSEA_plots_number=30,
    collapse='Remap_Only'
){
  work_dir <- getwd()
  
  GSEA_installed_path.new <-  path_syntax_change(GSEA_installed_path)
  
  rnks <- list.files(DEG_path,recursive=T,pattern = '.rnk',full.names = T,include.dirs=T)
  
  if(grep(work_dir,rnks[[1]]) %>% length() == 0){
    rnks <- paste0(work_dir,'/',rnks)
  }
  if(grep(work_dir,out_dir) %>% length() == 0){
    out_dir <- paste0(work_dir,'/',out_dir)
  }
  
  
  
  
  batch <- function(rnk,gene_set,out_dir=out_dir,name){
    
    rnk_names <- strsplit(rnk,'\\/')[[1]]
    
    comparison_folder <- paste0(out_dir,'/',rnk_names[length(rnk_names)] )
    
    if(file.exists(comparison_folder) != T){
      dir.create(comparison_folder)
    }
    
    out_dir.new <- path_syntax_change(comparison_folder)
    rnk.new <- path_syntax_change(rnk)
   
    command <- command.rnk
    
    #Replace Values
    command <- gsub("_rnk_file_", rnk.new, command)
    command <- gsub("_chip_", symbol_chip, command)
    command <- gsub("_gene_set_", gene_set, command)
    
    if(out_dir == 'auto'){
      command <- gsub("_out_dir_",dirname(rnk) %>% path_syntax_change, command) 
    }else{
      command <- gsub("_out_dir_", out_dir.new, command)
    }
    
    
    
    
    command <- gsub("_the_label_", paste0('GSEA_result_',gsub(' ','_',name)), command,fixed = T)
    
    working_dir  <- gsub("_cd_",GSEA_installed_path.new, working_dir )
    command <- gsub("_numberplot_",GSEA_plots_number, command)
    command <- gsub("_collapse_",collapse, command)
    command <- gsub("_species_",species, command)
    print(paste0('Performing enrichment analysis of ---  ',basename(rnk),'  --- with ---  ' ,name,' : ',gene_set,' .....' ))
    Sys.setenv(JAVA_HOME=paste0(GSEA_installed_path,'/jdk-11'))
    system("cmd.exe",input=c(working_dir,command),show.output.on.console = F,wait = T)
    Sys.sleep(1)
  }
  for(i in 1:length(rnks)){
    for(name in names(gene_sets) ){
      batch(rnks[i],gene_sets[name],out_dir,name)}}
 
}









GSEA_batch.counts <- function(
    GSEA_installed_path = "E:/bioinfo/GSEA_4.3.2",
    counts_path ='counts',
    species = 'mouse',
    gene_sets = c(`hallmark gene sets`='mh.all.v2024.1.Mm.symbols.gmt',
                 `positional gene sets`='m1.all.v2024.1.Mm.symbols.gmt',
                 `curated gene sets`='m2.all.v2024.1.Mm.symbols.gmt',
                 `regulatory target gene sets`='m3.all.v2024.1.Mm.symbols.gmt',
                 `ontology gene sets`='m5.go.v2024.1.Mm.symbols.gmt',
                 `cell type signature gene sets`='m8.all.v2024.1.Mm.symbols.gmt',
                 `wiki pathways` = 'm2.cp.reactome.v2024.1.Mm.symbols.gmt',
                 `reactome gene sets` = 'm2.cp.wikipathways.v2024.1.Mm.symbols.gmt',
                 `GOBP gene sets`='m5.go.bp.v2024.1.Mm.symbols.gmt'),
    symbol_chip = 'Mouse_Gene_Symbol_Remapping_MSigDB.v2024.1.Mm.chip',
    out_dir = "E:/bioinfo/Irineos_mice_scRNA/old threshold, more cell/GSEA/results/aggregated/save",
    GSEA_plots_number = 30,
    collapse ='Remap_Only',
    cls_path = "E:/bioinfo/Irineos_mice_scRNA/old threshold, more cell/GSEA/cls"
){
  
  count_files <- list.files(counts_path,recursive=T,pattern = '.txt',full.names = T,include.dirs=T)
  cls_files <- list.files(cls_path,recursive=T,pattern = '.cls',full.names = T,include.dirs=T)
  
  GSEA_installed_path.new <-  path_syntax_change(GSEA_installed_path)
  out_dir.new <- path_syntax_change(out_dir)
  
  
  batch <- function(count=count_files[1],cls=cls_files[1],gene_set,out_dir=out_dir,name){
    
    comparison <- basename(cls)
    dt_name <- basename(count)
    
    assay.label <- paste( gsub(' ','_',name),str_remove(dt_name, '.txt') ,str_remove(comparison, '.cls'), sep = '_' )
    assay.label <- gsub(' ','_',assay.label)
    
    
    cls.extension <- str_remove(comparison, '.cls') %>% str_replace('_vs_','_versus_')
    cls <- paste0(cls,'#',cls.extension)
    #print(cls)
    
    cls <- path_syntax_change(cls)
    count <- path_syntax_change(count)
        
    command <- command.counts

    #Replace Values
    command <- gsub("_count_file_", count, command)
    command <- gsub("_comparison_", cls, command)
    command <- gsub("_chip_", symbol_chip, command)
    command <- gsub("_gene_set_", gene_set, command)
    
    if(out_dir == 'auto'){
      command <- gsub("_out_dir_",dirname(rnk) %>% path_syntax_change, command) 
    }else{
      command <- gsub("_out_dir_", out_dir.new, command)
    }
    
    
    command <- gsub("_the_label_",assay.label , command,fixed = T)
    working_dir  <- gsub("_cd_",GSEA_installed_path.new, working_dir )
    command <- gsub("_numberplot_",GSEA_plots_number, command)
    command <- gsub("_collapse_",collapse, command)
    command <- gsub("_species_",species, command)
    
    print(paste0('Performing enrichment analysis of ---  ',str_remove(dt_name, '.txt'),' + ',str_remove(comparison, '.cls'),'  --- with ---  ' ,name,' : ',gene_set,' .....' ))
    Sys.setenv(JAVA_HOME=paste0(GSEA_installed_path,'/jdk-11'))
    system("cmd.exe",input=c(working_dir,command),show.output.on.console = F,wait = T)
    Sys.sleep(1)
  }
  for(i in 1:length(count_files)){
    for(n in 1:length(cls_files)){
      for(name in names(gene_sets) ){
        batch(count_files[i],cls_files[n],gene_sets[name],out_dir,name)
      }
    }
  }
}


# command <- 'cd _path_ 
# cd E:\\\\'
#   
# command <- gsub('_path_',out_dir, command)
# 
# system("cmd.exe",input=command,show.output.on.console = T)




