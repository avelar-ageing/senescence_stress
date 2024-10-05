#Script to download and analyse temporal DEGs
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')

cs_signatures_studies=read.csv(paste0(save_dir_csv,'1_cellage.csv'))
cs_signatures_studies=cs_signatures_studies[cs_signatures_studies$source!='Tao et al',]
#####

ensembl100=useMart(host='https://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)

human_pc_entrez=getBM(attributes=c('external_gene_name','entrezgene_id','ensembl_gene_id'),
                      filters = 'biotype',
                      values = c('protein_coding'),
                      mart = ensembl100)
#####
cs_signatures_studies$dir[cs_signatures_studies$dir=='Induces']='Induce CS'
cs_signatures_studies$dir[cs_signatures_studies$dir=='Inhibits']='Inhibit CS'
cs_signatures_studies$dir[cs_signatures_studies$dir=='Overexpressed']='Up in RS'
cs_signatures_studies$dir[cs_signatures_studies$dir=='Underexpressed']='Down in RS'

#Temporal DEG download
#time series
search_1_time=build_search(column='sra',
                           search_terms=c("ERP021140"))

recount_filtered_quiescence_fixed=read.csv('/Volumes/GoogleDrive/My Drive/PhD_copy/To Publish/CS Clusters/CS_studies/recount_final.csv')
# Keratinocyte
time_analysis(cell_type='Keratinocyte',
              recount_pheno=recount_filtered_quiescence_fixed,
              ensembl_dictionary=human_pc,
              grab_all=FALSE,
              enrich_degs_test=FALSE,
              entrez_dictionary=human_pc_entrez,
              exclude_search=NULL,
              ensembl_gene_col='ensembl_gene_id',
              treatment_time_col='time_after_treatment',
              time_levels=c('none', '4_days', '10_days', '20_days'),
              pval_cutoff=0.05,
              log2fc_cutoff=log2(1.5),
              independentFiltering=TRUE,
              db=cellage,
              db_col='Senescence',
              exclude_db=NULL,
              search_term_study=search_1_time,
              sep='/',
              results_dir=save_dir,
              fix_batch=TRUE,
              batch_1=c('ERR1805235','ERR1805236','ERR1805238','ERR1805230','ERR1805231','ERR1805224',
                        'ERR1805223','ERR1805222','ERR1805239','ERR1805240','ERR1805241','ERR1805229'),
              run_simulations = TRUE,
              simulation_n = 10000,
              pca_plot_group=c('cell_state','time_after_treatment'),
              vst_n=500)

# Melanocyte
time_analysis(cell_type='Melanocyte',
              recount_pheno=recount_filtered_quiescence_fixed,
              ensembl_dictionary=human_pc,
              grab_all=FALSE,
              enrich_degs_test=FALSE,
              entrez_dictionary=human_pc_entrez,
              exclude_search=NULL,
              ensembl_gene_col='ensembl_gene_id',
              treatment_time_col='time_after_treatment',
              time_levels=c('none', '4_days', '10_days', '20_days'),
              pval_cutoff=0.05,
              log2fc_cutoff=log2(1.5),
              independentFiltering=TRUE,
              db=cellage,
              db_col='Senescence',
              exclude_db=NULL,
              search_term_study=search_1_time,
              sep='/',
              results_dir=save_dir,
              fix_batch=FALSE,
              run_simulations = TRUE,
              simulation_n = 10000,
              pca_plot_group=c('cell_state','time_after_treatment'),
              vst_n=500)

# Fibroblast
time_analysis(cell_type='Fibroblast',
              recount_pheno=recount_filtered_quiescence_fixed,
              ensembl_dictionary=human_pc,
              grab_all=FALSE,
              enrich_degs_test=FALSE,
              entrez_dictionary=human_pc_entrez,
              exclude_search=NULL,
              ensembl_gene_col='ensembl_gene_id',
              treatment_time_col='time_after_treatment',
              time_levels=c('none', '4_days', '10_days', '20_days'),
              pval_cutoff=0.05,
              log2fc_cutoff=log2(1.5),
              independentFiltering=TRUE,
              db=cellage,
              db_col='Senescence',
              exclude_db=NULL,
              search_term_study=search_1_time,
              sep='/',
              results_dir=save_dir,
              fix_batch=FALSE,
              run_simulations = TRUE,
              simulation_n = 10000,
              pca_plot_group=c('cell_state','time_after_treatment'),
              vst_n=500)

#####
#All DEG analyses
# time_degs=compile_time_files(dir_use='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/',
#                    deg_file='degs.csv')
# save_csv(time_degs,file_name = 'all_time_degs',
#          path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')
time_degs=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/all_time_degs.csv')
time_degs=factor_column_and_modify(df = time_degs,
                                   column = 'group_1',
                                   old_list = rev(c('4_days','10_days','20_days')),
                                   keyword = NULL)

time_degs=factor_column_and_modify(df = time_degs,
                                   column = 'cell_type',
                                   old_list = c('fibroblast','keratinocyte',
                                                'melanocyte'),
                                   keyword = NULL)
#####
#simulations
sims=compile_time_files(dir_use='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/',
                   deg_file='simulation_overlap.csv')
save_csv(sims,file_name = 'sims_by_type',path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/')
sims_cumsum=compile_time_files(dir_use='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/',
                               deg_file='reverse_cumsum_recount.csv')
save_csv(sims_cumsum,file_name = 'reverse_cumsum_recount',path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/')

#####
time_degs$ingroup=paste0(time_degs$cell_type,'_',
                         time_degs$group_1_dir)

time_deg_overlap_sim=simulate_overlaps(relevant_degs = time_degs,
                  facet_col='accession',
                  ingroup_col='ingroup',
                  outgroup_col='group_2_dir',
                  outgroup_relevant=c('none_up',
                                      'none_down'),
                  gene_col='ensembl',
                  seed=1,
                  simulation_n=10000)

save_csv(time_deg_overlap_sim,file_name = 'all_sim_overlap',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

max(time_deg_overlap_sim$Freq[time_deg_overlap_sim$outgroup=='none_up'&
                                time_deg_overlap_sim$Var1==9])
max(time_deg_overlap_sim$Var1[time_deg_overlap_sim$outgroup=='none_up'&
                                time_deg_overlap_sim$Freq!=0])
max(time_deg_overlap_sim$Var1[time_deg_overlap_sim$outgroup=='none_down'&
                                time_deg_overlap_sim$Freq!=0])

max(time_deg_overlap_sim$Freq[time_deg_overlap_sim$outgroup=='none_down'&
                                time_deg_overlap_sim$Var1==9])

time_degs_sig=time_degs[time_degs$sig=='y',]

table(time_degs_sig%>%dplyr::select(group_1,direction_1,cell_type))
#####
#Plot DEG overlaps with each other
time_degs$direction_1=gsub(time_degs$direction_1,pattern='up',
                           replacement='↑')
time_degs$direction_1=gsub(time_degs$direction_1,pattern='down',
                           replacement='↓')
time_degs$direction_2=gsub(time_degs$direction_2,pattern='up',
                           replacement='↑')
time_degs$direction_2=gsub(time_degs$direction_2,pattern='down',
                           replacement='↓')
time_deg_overlaps = compare_degs_between_groups(
  all_degs = time_degs,
  group_comparison_col_1 = "group_1",
  group_comparison_col_2 = "group_2",
  outgroup_col = "outgroup_col_default",
  group_col = "cell_type",
  gene_col = 'ensembl',
  accession_col = 'dir_accession',
  return_overlap = TRUE,
  outgroup = "none", 
  direction_col_1 = 'direction_1',
  direction_col_2 = 'direction_2',
  xlab = NULL,
  ylab = NULL,
  temp_order = c('4 Days','10 Days','20 Days'),
  graph_scale = 2
)

save_csv(time_deg_overlaps$df,file_name = 'temporal_comparison',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

save_csv(time_deg_overlaps$df,file_name = 'temp_temporal_comparison',path = '/Users/ravelarvargas/Downloads')

save_p(plot = time_deg_overlaps$p,
       file_name = 'temporal_comparison',
       save_dir ='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 13,
       p_height = 11)

#####
#just interesting overlaps
all_degs = time_degs
group_comparison_col_1 = "group_1"
group_comparison_col_2 = "group_2"
outgroup_col = "outgroup_col_default"
group_col = "cell_type"
gene_col = 'ensembl'
accession_col = 'dir_accession'
return_overlap = TRUE
outgroup = "none" 
direction_col_1 = 'direction_1'
direction_col_2 = 'direction_2'
xlab = NULL
ylab = NULL
temp_order = c('4 Days','10 Days','20 Days')
graph_scale = 2

all_comparison=unique(all_degs[[group_col]])
all_degs[[group_comparison_col_1]]=as.character(all_degs[[group_comparison_col_1]])
all_degs[[group_comparison_col_2]]=as.character(all_degs[[group_comparison_col_2]])
all_degs[[group_col]]=as.character(all_degs[[group_col]])
comparison_df=vector_to_df(all_comparison)
all_overlaps=do.call('rbind',apply(comparison_df,MARGIN=1,function(per_row){
  group_1=per_row[['x']]
  group_2=per_row[['y']]
  
  degs_1=all_degs[all_degs[[group_col]]==group_1,]
  degs_1_relevant=degs_1[degs_1[[group_comparison_col_2]]==outgroup&degs_1[['sig']]=='y',]
  degs_2=all_degs[all_degs[[group_col]]==group_2,]
  degs_2_relevant=degs_2[degs_2[[group_comparison_col_2]]==outgroup&degs_2[['sig']]=='y',]
  
  background_temp=unique(c(degs_1[[gene_col]],
                           degs_2[[gene_col]]))
  
  accession_1=unique(degs_1_relevant[[accession_col]])
  accession_2=unique(degs_2_relevant[[accession_col]])
  overlap_temp=do.call('rbind',lapply(accession_1,function(per_accession_1){
    degs_1_relevant_temp=degs_1_relevant[degs_1_relevant[[accession_col]]==per_accession_1,]
    degs_1_relevant_temp_gene=degs_1_relevant_temp[[gene_col]]
    
    n_1=length(degs_1_relevant_temp_gene)
    n_1_name=degs_1_relevant_temp[[group_col]][1]
    n_1_dir_ingroup=degs_1_relevant_temp[[direction_col_1]][1]
    n_1_dir_time=degs_1_relevant_temp[[group_comparison_col_1]][1]
    do.call('rbind',lapply(accession_2,function(per_accession_2){
      degs_2_relevant_temp=degs_2_relevant[degs_2_relevant[[accession_col]]==per_accession_2,]
      degs_2_relevant_temp_gene=degs_2_relevant_temp[[gene_col]]
      overlap=testGeneOverlap(newGeneOverlap(listA=degs_1_relevant_temp_gene,
                                             listB=degs_2_relevant_temp_gene,
                                             genome.size=length(background_temp)))
      con_tbl=getContbl(overlap)
      actual=con_tbl[4]
      expected=chisq.test(con_tbl)$exp[4]
      n_diff=actual-expected
      odds=fisher.test(con_tbl)$est
      pval=fisher.test(con_tbl)$p.val
      overlap=paste0(getIntersection(overlap),collapse='/')
      
      n_2=length(degs_2_relevant_temp_gene)
      
      n_2_name=degs_2_relevant_temp[[group_col]][1]
      n_2_dir_ingroup=degs_2_relevant_temp[[direction_col_1]][1]
      n_2_dir_time=degs_2_relevant_temp[[group_comparison_col_1]][1]
      
      if(return_overlap){
        full_return=matrix(c(actual,expected,n_diff,odds,pval,
                             n_1_name,n_1_dir_ingroup,n_1_dir_time,n_1,
                             n_2_name,n_2_dir_ingroup,n_2_dir_time,n_2,
                             length(background_temp),overlap),ncol=15)
        colnames(full_return)=c('actual','expected','diff','odds','pval',
                                'group_1','dir_ingroup_1','time_1','n_1',
                                'group_2','dir_ingroup_2','time_2','n_2',
                                'background_n','overlap')
      }else{
        full_return=matrix(c(actual,expected,n_diff,odds,pval,
                             n_1_name,n_1_dir_ingroup,n_1_dir_time,n_1,
                             n_2_name,n_2_dir_ingroup,n_2_dir_time,n_2,
                             length(background_temp)),ncol=14)
        colnames(full_return)=c('actual','expected','diff','odds','pval',
                                'group_1','dir_ingroup_1','time_1','n_1',
                                'group_2','dir_ingroup_2','time_2','n_2',
                                'background_n')
      }
      
      return(full_return)
    }))
  }))
  overlap_temp=data.frame(overlap_temp)
  return(overlap_temp)
}))
all_overlaps[['actual']]=as.numeric(as.character(all_overlaps[['actual']]))
all_overlaps[['expected']]=as.numeric(as.character(all_overlaps[['expected']]))
all_overlaps[['diff']]=as.numeric(as.character(all_overlaps[['diff']]))
all_overlaps[['pval']]=as.numeric(as.character(all_overlaps[['pval']]))
all_overlaps[['odds']]=as.numeric(as.character(all_overlaps[['odds']]))
all_overlaps[['n_1']]=as.numeric(as.character(all_overlaps[['n_1']]))
all_overlaps[['n_2']]=as.numeric(as.character(all_overlaps[['n_2']]))
all_overlaps[['background_n']]=as.numeric(as.character(all_overlaps[['background_n']]))

all_overlaps$p.adj=p.adjust(all_overlaps[['pval']])
all_overlaps$log2odds=log2(all_overlaps[['odds']])
all_overlaps=all_overlaps%>%rstatix::add_significance('p.adj')

all_overlaps[['accession']]=paste0(all_overlaps[['group_1']],' ',
                                   all_overlaps[['dir_ingroup_1']],
                                   '\n',all_overlaps[['group_2']],' ',
                                   all_overlaps[['dir_ingroup_2']]
)
all_overlaps[['x']]=paste0(all_overlaps[['group_1']],' ',gsub(all_overlaps[['time_1']],pattern='_',replacement=' '))
all_overlaps[['y']]=paste0(all_overlaps[['group_2']],' ',gsub(all_overlaps[['time_2']],pattern='_',replacement=' '))

#Fix some legend parameters
graph_max=max(ceiling(all_overlaps$log2odds/2)[!is.infinite(ceiling(all_overlaps$log2odds))])*2
graph_min=min(floor(all_overlaps$log2odds)[!is.infinite(ceiling(all_overlaps$log2odds))])
if(graph_min>0){
  graph_min=0
}
li <- c(graph_min, graph_max)
la <- c(seq(graph_min,graph_max,graph_scale))
br <- c(seq(graph_min,graph_max,graph_scale))

if(!is.null(temp_order)){
  x_relevant=unique(all_overlaps[['x']])
  x_order=do.call(c,lapply(temp_order,function(per_order){
    x_relevant[grep(x_relevant,pattern=per_order)]
  }))
  all_overlaps[['x']]=factor(all_overlaps[['x']],
                             levels=x_order)
  
  y_relevant=unique(all_overlaps[['y']])
  y_order=do.call(c,lapply(temp_order,function(per_order){
    y_relevant[grep(y_relevant,pattern=per_order)]
  }))
  all_overlaps[['y']]=factor(all_overlaps[['y']],
                             levels=rev(y_order))
}

text_size=12

accession_all=unique(all_overlaps$accession)

temp_p_overlaps=all_overlaps[grepl(all_overlaps$accession,pattern='Melanocyte'),]
temp_p_overlaps[temp_p_overlaps$accession=='Keratinocyte ↓\nMelanocyte ↑',][['accession']]='Melanocyte ↑\nKeratinocyte ↓'
temp_p_overlaps[temp_p_overlaps$accession=='Keratinocyte ↑\nMelanocyte ↑',][['accession']]='Melanocyte ↑\nKeratinocyte ↑'
temp_p_overlaps[temp_p_overlaps$accession=='Keratinocyte ↓\nMelanocyte ↓',][['accession']]='Melanocyte ↓\nKeratinocyte ↓'
temp_p_overlaps[temp_p_overlaps$accession=='Keratinocyte ↑\nMelanocyte ↓',][['accession']]='Melanocyte ↓\nKeratinocyte ↑'
# temp_p_overlaps$accession=factor(temp_p_overlaps$accession,
#                                  levels=c('Melanocyte ↑\nFibroblast ↓',
#                                           'Melanocyte ↓\nFibroblast ↑',
#                                           'Melanocyte ↑\nKeratinocyte ↓',
#                                           'Melanocyte ↓\nKeratinocyte ↑'))


accession_use=do.call('rbind',apply(temp_p_overlaps,MARGIN = 1,function(get_row_column){
  accession_x=ifelse(grepl(pattern='Melanocyte',get_row_column[['x']]),TRUE,FALSE)
  if(accession_x){
    temp_df=data.frame(matrix(c(get_row_column[['x']],
                                get_row_column[['y']]),ncol=2))
  }else{
    temp_df=data.frame(matrix(c(get_row_column[['y']],
                get_row_column[['x']]),ncol=2))
  }
  colnames(temp_df)=c('x_use','y_use')
  return(temp_df)
}))
temp_p_overlaps=cbind(temp_p_overlaps,accession_use)
temp_p_overlaps$x_use=factor(temp_p_overlaps$x_use,
                             levels=c('Melanocyte 4 Days',
                                      'Melanocyte 10 Days',
                                      'Melanocyte 20 Days'))
temp_p_overlaps$y_use=factor(temp_p_overlaps$y_use,
                             levels=rev(c('Fibroblast 4 Days',
                                      'Keratinocyte 4 Days',
                                      'Fibroblast 10 Days',
                                      'Keratinocyte 10 Days',
                                      'Fibroblast 20 Days',
                                      'Keratinocyte 20 Days')))

temp_p_overlaps$p.adj.signif[temp_p_overlaps$p.adj.signif=='ns']=''
p_temp=temp_p_overlaps%>%ggplot(aes(x=x_use,y=y_use,fill=log2odds))+
  geom_tile(colour='black')+
  geom_text(aes(label=p.adj.signif),nudge_y = 0.25, size = text_size / 2)+
  geom_text(aes(label=actual),nudge_y=-0.25, size = text_size / 2)+
  scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = TRUE, 
                                              barwidth = 8,
                                              ticks.colour = 'black',
                                              title.position = 'top',title.hjust=0.5),
                       breaks=br,
                       labels=la,
                       limits=li)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  facet_wrap(~accession,scales='free')+
  theme(axis.text.x=element_text(angle=45,hjust=1, size = text_size),
        axis.text.y = element_text(size = text_size),
        legend.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        strip.text = element_text(size = text_size),
        legend.position='bottom')+
  xlab(xlab)+
  ylab(ylab)
# 
# save_p(p_temp,
#        file_name = 'temporal_comparison_melanocyte',
#        save_dir ='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
#        p_width = 7,
#        p_height = 7)

#simplified
temp_p_overlaps=temp_p_overlaps[,colnames(temp_p_overlaps)!='overlap']

library(stringr)
# Use extract to separate into four columns based on the pattern
df <- temp_p_overlaps %>%
  extract(accession, into = c("CellType1", "Arrow1", "CellType2", "Arrow2"), 
          regex = "([a-zA-Z]+)\\s(↑|↓)\\s([a-zA-Z]+)\\s(↑|↓)") %>%
  mutate(Arrow1 = ifelse(Arrow1 == "↑", "U", "D"),
         Arrow2 = ifelse(Arrow2 == "↑", "U", "D"))
df$Arrow1=paste0(df$time_1,' ',df$Arrow1)
df$Arrow2=paste0(df$time_2,' ',df$Arrow2)

df$fill=ifelse(df$p.adj.signif=='','Not Sig',
               ifelse(df$log2odds<0,'Sig Under-represented','Sig Over-represented'))
df$Arrow1=factor(df$Arrow1,
                 levels=rev(c('4 Days U','10 Days U','20 Days U',
                          '4 Days D', '10 Days D','20 Days D')))

df$Arrow2=factor(df$Arrow2,
                 levels=c('4 Days U','10 Days U','20 Days U',
                              '4 Days D', '10 Days D','20 Days D'))

fill_manual=c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
              "Sig Under-represented" = "#00bfc4")
summary_melanocyte=df%>%ggplot(aes(x=Arrow2,y=Arrow1,fill=fill))+
  geom_tile(colour='black')+
  scale_fill_manual(name=NULL,
                    values = fill_manual)+
  facet_grid(~CellType2)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  geom_hline(yintercept = 3.5,size=1)+
  geom_vline(xintercept = 3.5,size=1)+
  ylab('Melanocyte Temporal DEG Direction')+
  xlab('Temporal DEG Direction')+
  theme(legend.position='bottom',
    strip.text = element_text(size = 12),   # Facet title text size
    axis.text = element_text(size = 10),    # Axis text size
    axis.title = element_text(size = 12),   # Axis title text size
    legend.text = element_text(size = 10),  # Legend text size
    legend.title = element_text(size = 12)  # Legend title text size
    )
save_p(summary_melanocyte,
       file_name = 'temporal_self_summary',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',p_width = 10,
       p_height = 4)

temp_p_overlaps_no_melo=all_overlaps[all_overlaps$group_1!='Melanocyte',]
temp_p_overlaps_no_melo=temp_p_overlaps_no_melo[temp_p_overlaps_no_melo$group_2!='Melanocyte',]

temp_p_overlaps_no_melo=temp_p_overlaps_no_melo[,colnames(temp_p_overlaps_no_melo)!='overlap']

library(stringr)
# Use extract to separate into four columns based on the pattern
df <- temp_p_overlaps_no_melo %>%
  extract(accession, into = c("CellType1", "Arrow1", "CellType2", "Arrow2"), 
          regex = "([a-zA-Z]+)\\s(↑|↓)\\s([a-zA-Z]+)\\s(↑|↓)") %>%
  mutate(Arrow1 = ifelse(Arrow1 == "↑", "U", "D"),
         Arrow2 = ifelse(Arrow2 == "↑", "U", "D"))
df$Arrow1=paste0(df$time_1,' ',df$Arrow1)
df$Arrow2=paste0(df$time_2,' ',df$Arrow2)

df$fill=ifelse(df$p.adj.signif=='','Not Sig',
               ifelse(df$log2odds<0,'Sig Under-represented','Sig Over-represented'))
df$Arrow1=factor(df$Arrow1,
                 levels=rev(c('4 Days U','10 Days U','20 Days U',
                              '4 Days D', '10 Days D','20 Days D')))

df$Arrow2=factor(df$Arrow2,
                 levels=c('4 Days U','10 Days U','20 Days U',
                          '4 Days D', '10 Days D','20 Days D'))

fill_manual=c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
              "Sig Under-represented" = "#00bfc4")
summary_no_melanocyte=df%>%ggplot(aes(x=Arrow2,y=Arrow1,fill=fill))+
  geom_tile(colour='black')+
  scale_fill_manual(name=NULL,
                    values = fill_manual)+
  facet_grid(CellType1~CellType2)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  geom_hline(yintercept = 3.5,size=1)+
  geom_vline(xintercept = 3.5,size=1)+
  ylab('Melanocyte Temporal DEG Direction')+
  xlab('Temporal DEG Direction')+
  theme(legend.position='bottom',
        strip.text = element_text(size = 12),   # Facet title text size
        axis.text = element_text(size = 10),    # Axis text size
        axis.title = element_text(size = 12),   # Axis title text size
        legend.text = element_text(size = 10),  # Legend text size
        legend.title = element_text(size = 12)  # Legend title text size
  )

##
#all
library(stringr)
# Use extract to separate into four columns based on the pattern
df <- all_overlaps %>%
  extract(accession, into = c("CellType1", "Arrow1", "CellType2", "Arrow2"), 
          regex = "([a-zA-Z]+)\\s(↑|↓)\\s([a-zA-Z]+)\\s(↑|↓)") %>%
  mutate(Arrow1 = ifelse(Arrow1 == "↑", "U", "D"),
         Arrow2 = ifelse(Arrow2 == "↑", "U", "D"))
df$Arrow1=paste0(df$time_1,' ',df$Arrow1)
df$Arrow2=paste0(df$time_2,' ',df$Arrow2)

df$fill=ifelse(df$p.adj.signif=='','Not Sig',
               ifelse(df$log2odds<0,'Sig Under-represented','Sig Over-represented'))

fill_manual=c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
              "Sig Under-represented" = "#00bfc4")

df$CellType1=paste0(df$CellType1,'\nTemporal DEGs')
df$CellType2=paste0(df$CellType2,'\nTemporal DEGs')
all_accession=unique(c(df$CellType1,df$CellType2))
df$Arrow2=as.character(df$Arrow2)
df$Arrow1=as.character(df$Arrow1)
for(i in 1:nrow(df)){
  if(grepl(df$CellType2[i],pattern='Melanocyte')){
    df$CellType2[i]=df$CellType1[i]
    df$CellType1[i]="Melanocyte\nTemporal DEGs"
    temp=df$Arrow2[i]
    
    df$Arrow2[i]=df$Arrow1[i]
    df$Arrow1[i]=temp
  }
}

df$Arrow1=factor(df$Arrow1,
                 levels=rev(c('4 Days U','10 Days U','20 Days U',
                              '4 Days D', '10 Days D','20 Days D')))

df$Arrow2=factor(df$Arrow2,
                 levels=c('4 Days U','10 Days U','20 Days U',
                          '4 Days D', '10 Days D','20 Days D'))
all_p=lapply(all_accession,function(exclude_me){
  df_temp=df[df$CellType1!=exclude_me,]
  df_temp=df_temp[df_temp$CellType2!=exclude_me,]
  df_temp$fill[df_temp$p.adj>0.05]='Not Sig'
  df_temp%>%ggplot(aes(x=Arrow2,y=Arrow1,fill=fill))+
    geom_tile(colour='black')+
    scale_fill_manual(name=NULL,
                      values = fill_manual)+
    facet_grid(CellType1~CellType2)+
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    geom_hline(yintercept = 3.5,size=1)+
    geom_vline(xintercept = 3.5,size=1)+
    ylab(NULL)+
    xlab(NULL)+
    theme(legend.position='bottom',
          strip.text = element_text(size = 13),   # Facet title text size
          axis.text = element_text(size = 13),    # Axis text size
          axis.title = element_text(size = 13),   # Axis title text size
          legend.text = element_text(size = 11),  # Legend text size
          legend.title = element_text(size = 13)  # Legend title text size
    )
})

lapply(all_p,function(save_p){
  save_p(save_p,
         file_name = paste0(sample(1:1000)[1]),
         save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
         p_width = 7,p_height = 4)
})

simple_all_temporal=all_p[[1]]+patchwork::plot_spacer()+all_p[[2]]+all_p[[3]]+patchwork::plot_layout(guides = 'collect')&
  theme(legend.position = 'bottom')

save_p(simple_all_temporal,file_name = 'temporal_simple',save_dir ='/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
       p_width = 14,p_height = 5.5)

#####

#Plot DEG overlaps with CellAge
time_degs_sig$cell_type=as.character(time_degs_sig$cell_type)
time_degs_sig$group_1=as.character(time_degs_sig$group_1)
cs_signatures_studies$group_use=paste0(cs_signatures_studies$source,'_',
                                       cs_signatures_studies$dir)
cs_signatures_studies$group_use[cs_signatures_studies$group_use=='Saul et al_']='SenMayo'
time_degs_sig$accession=paste0(time_degs_sig$cell_type,'_',
                               time_degs_sig$group_1_dir)
db_deg_overlaps_time=overlap_function(df_1 = time_degs_sig,
                                      gene_col_1 = 'ensembl',
                                      group_col_1 = 'accession',
                                      carry_col_1 = c('cell_type','direction_1','group_1'),
                                      df_2 = cs_signatures_studies,
                                      group_col_2 = 'group_use',
                                      gene_col_2 = 'gene',
                                      background = time_degs$ensembl,
                                      carry_col_2 = c('dir','source'))

db_deg_overlaps_time$facet=paste0(db_deg_overlaps_time$cell_type,'\n',
                                  db_deg_overlaps_time$direction_1)

db_deg_overlaps_time$group_1=gsub(
  db_deg_overlaps_time$group_1,pattern='_',replacement=' '
)
db_deg_overlaps_time$facet=gsub(db_deg_overlaps_time$facet,
                                pattern='down',replacement='↓')
db_deg_overlaps_time$facet=gsub(db_deg_overlaps_time$facet,
                                pattern='up',replacement='↑')
db_deg_overlaps_time$facet=
  factor(db_deg_overlaps_time$facet,levels=unique(db_deg_overlaps_time$facet)[c(which(grepl(unique(db_deg_overlaps_time$facet),
                                                       pattern='↑')),
                                           which(grepl(unique(db_deg_overlaps_time$facet),
                                                       pattern='↓')))])
db_deg_overlaps_time$group_1=factor(db_deg_overlaps_time$group_1,
                                    levels=c('4 Days','10 Days','20 Days'))

db_deg_overlaps_time$group_use=gsub(db_deg_overlaps_time$group_use,pattern='_',replacement=' ')
db_deg_overlaps_time$group_use=gsub(db_deg_overlaps_time$group_use,
                                    pattern='Signature',replacement='CellAge')
db_deg_overlaps_time$group_use=gsub(db_deg_overlaps_time$group_use,
                                    pattern='Driver Inhibit',replacement='CellAge Inhibits')
db_deg_overlaps_time$group_use=gsub(db_deg_overlaps_time$group_use,
                                    pattern='Driver Induce',replacement='CellAge Induces')
db_deg_overlaps_time$group_use=gsub(db_deg_overlaps_time$group_use,
                                    pattern='down',replacement='Down')
db_deg_overlaps_time$group_use=gsub(db_deg_overlaps_time$group_use,
                                    pattern='up',replacement='Up')

db_deg_overlaps_time=db_deg_overlaps_time[!grepl(db_deg_overlaps_time$group_use,
                                                 pattern='Tao'),]

db_deg_overlaps_time$group_use=
       gsub(db_deg_overlaps_time$group_use,pattern='Induces',replacement='Induce')
db_deg_overlaps_time$group_use=
            gsub(db_deg_overlaps_time$group_use,pattern='Inhibits',replacement='Inhibit')
db_deg_overlaps_time$group_use=factor(db_deg_overlaps_time$group_use,
                                     levels=rev(c('CellAge Induce CS',
                                                  'CellAge Inhibit CS',
                                                  'CellAge Up in RS',
                                              'CellAge Down in RS',
                                              'SenMayo',
                                              'Hernandez-Segura et al Up',
                                              'Hernandez-Segura et al Down',
                                              'Casella et al Up',
                                              'Casella et al Down',
                                              'Cherry et al Up',
                                              'Cherry et al Down')))

compiled_time_db_overlaps=create_overlap_plot(deg_db_overlap = db_deg_overlaps_time,
                    facet_col = 'facet',
                    y = 'group_use',x='group_1',
                    xlab = 'Temporal DEGs',
                    ylab='CS Gene Lists')

save_p(compiled_time_db_overlaps,file_name = 'temporal_vs_cellage',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 9,p_height=9)

db_deg_overlaps_time$direction_1=gsub(db_deg_overlaps_time$direction_1,pattern='↑',
                           replacement='up')
db_deg_overlaps_time$direction_1=gsub(db_deg_overlaps_time$direction_1,pattern='↓',
                           replacement='down')
save_csv(db_deg_overlaps_time,file_name = 'temporal_vs_cellage',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

#simplify summary
db_deg_overlaps_time$summary_y=apply(db_deg_overlaps_time,MARGIN = 1,function(x){
  temp_x=ifelse(grepl(x[['direction_1']],pattern='up'),'Up','Down')
  paste0(x[['group_1']],' ',temp_x)
})

db_deg_overlaps_time$summary_y=gsub(db_deg_overlaps_time$summary_y,pattern='Down',replacement='D')
db_deg_overlaps_time$summary_y=gsub(db_deg_overlaps_time$summary_y,pattern='Up',replacement='U')
db_deg_overlaps_time$summary_y=factor(db_deg_overlaps_time$summary_y,
                                          levels=c("20 Days D",
                                                   "10 Days D",
                                                   "4 Days D",
                                                   "20 Days U",
                                                   "10 Days U",
                                                   "4 Days U"))

db_deg_overlaps_time$group_use=factor(db_deg_overlaps_time$group_use,
                                      levels=c('CellAge Induce CS',
                                               'CellAge Inhibit CS',
                                               'CellAge Up in RS',
                                               'CellAge Down in RS',
                                               'SenMayo',
                                               'Hernandez-Segura et al Up',
                                               'Hernandez-Segura et al Down',
                                               'Casella et al Up',
                                               'Casella et al Down',
                                               'Cherry et al Up',
                                               'Cherry et al Down'))
summary_p_temporal_stress=summarise_overlaps_facet(db=db_deg_overlaps_time,
                                                   x='group_use',
                                                   y='summary_y',
                                                   odds_col='odds',
                                                   pval_col='adj',
                                                   xlab='CS Genelists',
                                                   ylab='Temporal DEGs',
                                                   add_line=TRUE,
                                                   legend_title=NULL,
                                                   label_cq=FALSE,
                                                   facet_col='cell_type',angle = 45)

save_p(summary_p_temporal_stress,file_name = '6b_simplified_temporal_db',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
       p_width = 6,p_height = 4)
#####
time_degs$temporal_accession=paste0(time_degs$cell_type,' ',time_degs$group_1)
metabolism_list=read.csv('/Users/ravelarvargas/Downloads/marian/genelist_pc_arrest_filtered.csv')
metabolism_background=read.csv('/Users/ravelarvargas/Downloads/marian/background_pc_arrest_filtered.csv')
temporal_gsea=rank_gsea_wrapper(genes=time_degs,
                          condition_col='temporal_accession',
                          gene_col='ensembl',
                          pval_col='padj',
                          log2fc_col='log2FoldChange',
                          rank_type = 'both',
                          background=metabolism_background$x,
                          gene_list=metabolism_list,
                          gene_list_gene_col='gene',
                          gene_list_col_subset='class')

p_temp=plot_gsea_enrichment(enrichment_df=temporal_gsea$gsea_enrichment,
                     condition_col='temporal_accession',
                     pathway_col='ID',
                     factor_condition_levels=c(
                       'Fibroblast 4 Days',
                       'Fibroblast 10 Days',
                       'Fibroblast 20 Days',
                       'Melanocyte 4 Days',
                       'Melanocyte 10 Days',
                       'Melanocyte 20 Days',
                       'Keratinocyte 4 Days',
                       'Keratinocyte 10 Days',
                       'Keratinocyte 20 Days'
                     ))

temporal_gsea2=rank_gsea_wrapper(genes=time_degs,
                                condition_col='temporal_accession',
                                gene_col='ensembl',
                                pval_col='padj',
                                log2fc_col='log2FoldChange',
                                rank_type = 'both',
                                background=metabolism_background$x,
                                gene_list=metabolism_list,
                                gene_list_gene_col='gene',
                                gene_list_col_subset='pathway')

p2_temp=plot_gsea_enrichment(enrichment_df=temporal_gsea2$gsea_enrichment,
                     condition_col='temporal_accession',
                     pathway_col='ID',factor_condition_levels=c(
                       'Fibroblast 4 Days',
                       'Fibroblast 10 Days',
                       'Fibroblast 20 Days',
                       'Melanocyte 4 Days',
                       'Melanocyte 10 Days',
                       'Melanocyte 20 Days',
                       'Keratinocyte 4 Days',
                       'Keratinocyte 10 Days',
                       'Keratinocyte 20 Days'
                     ))

temporal_gsea3=rank_gsea_wrapper(genes=time_degs,
                                 condition_col='temporal_accession',
                                 gene_col='ensembl',
                                 pval_col='padj',
                                 log2fc_col='log2FoldChange',
                                 rank_type = 'both',
                                 background=metabolism_background$x,
                                 gene_list=metabolism_list,
                                 gene_list_gene_col='gene',
                                 gene_list_col_subset='class_2')

p3_temp=plot_gsea_enrichment(enrichment_df=temporal_gsea3$gsea_enrichment,
                     condition_col='temporal_accession',
                     pathway_col='ID',factor_condition_levels=c(
                       'Fibroblast 4 Days',
                       'Fibroblast 10 Days',
                       'Fibroblast 20 Days',
                       'Melanocyte 4 Days',
                       'Melanocyte 10 Days',
                       'Melanocyte 20 Days',
                       'Keratinocyte 4 Days',
                       'Keratinocyte 10 Days',
                       'Keratinocyte 20 Days'
                     ))

save_p(plot = p_temp,save_dir = '/Users/ravelarvargas/Downloads',
       file_name = 'test_1',p_width = 10)

save_p(plot = p2_temp,save_dir = '/Users/ravelarvargas/Downloads',
       file_name = 'test_2',p_width = 12)

save_p(plot = p3_temp,save_dir = '/Users/ravelarvargas/Downloads',
       file_name = 'test_3',p_width = 10)

time_degs$accession_full=paste0(time_degs$cell_type,' ',
                                time_degs$group_1,' ',time_degs$direction_1)

time_degs$accession_partial=paste0(time_degs$cell_type,' ',time_degs$group_1)
temp_overlap=overlap_function(df_1 = time_degs[time_degs$sig=='y',],
                 df_2 = metabolism_list,
                 gene_col_1 = 'ensembl',
                 gene_col_2 = 'gene',
                 carry_col_1 = c('accession_partial','direction_1'),
                 group_col_1 = 'accession_full',
                 group_col_2 = 'class_2',
                 background = metabolism_background$x)

temp_overlap$accession_partial=
  factor(temp_overlap$accession_partial,
         levels = rev(c('fibroblast 4_days','fibroblast 10_days',
                    'fibroblast 20_days','melanocyte 4_days',
                    'melanocyte 10_days','melanocyte 20_days',
                    'keratinocyte 4_days','keratinocyte 10_days',
                    'keratinocyte 20_days')))
p_4=create_overlap_plot(deg_db_overlap = temp_overlap,
                    facet_col = 'class_2',
                    y = 'accession_partial',x='direction_1',
                    add_almost_sig = TRUE,remove_nonsig = TRUE)

save_p(plot = p_4,save_dir = '/Users/ravelarvargas/Downloads',
       file_name = 'test_overlap',p_width = 10)

temp_overlap2=overlap_function(df_1 = time_degs[time_degs$sig=='y',],
                              df_2 = metabolism_list,
                              gene_col_1 = 'ensembl',
                              gene_col_2 = 'gene',
                              carry_col_1 = c('accession_partial','direction_1'),
                              group_col_1 = 'accession_full',
                              group_col_2 = 'class',
                              background = metabolism_background$x)

temp_overlap2$accession_partial=
  factor(temp_overlap2$accession_partial,
         levels = rev(c('fibroblast 4_days','fibroblast 10_days',
                        'fibroblast 20_days','melanocyte 4_days',
                        'melanocyte 10_days','melanocyte 20_days',
                        'keratinocyte 4_days','keratinocyte 10_days',
                        'keratinocyte 20_days')))

p_5=create_overlap_plot(deg_db_overlap = temp_overlap2,
                        facet_col = 'class',
                        y = 'accession_partial',x='direction_1',
                        add_almost_sig = TRUE,remove_nonsig = TRUE)

save_p(plot = p_5,save_dir = '/Users/ravelarvargas/Downloads',
       file_name = 'test_overlap2',p_width = 10,p_height=10)

temp_overlap3=overlap_function(df_1 = time_degs[time_degs$sig=='y',],
                               df_2 = metabolism_list,
                               gene_col_1 = 'ensembl',
                               gene_col_2 = 'gene',
                               carry_col_1 = c('accession_partial','direction_1'),
                               group_col_1 = 'accession_full',
                               group_col_2 = 'pathway',
                               background = metabolism_background$x)

temp_overlap3$accession_partial=
  factor(temp_overlap3$accession_partial,
         levels = rev(c('fibroblast 4_days','fibroblast 10_days',
                        'fibroblast 20_days','melanocyte 4_days',
                        'melanocyte 10_days','melanocyte 20_days',
                        'keratinocyte 4_days','keratinocyte 10_days',
                        'keratinocyte 20_days')))

p_6=create_overlap_plot(deg_db_overlap = temp_overlap3,
                        facet_col = 'pathway',
                        y = 'accession_partial',x='direction_1',
                        add_almost_sig = TRUE,remove_nonsig = TRUE)

save_p(plot = p_6,save_dir = '/Users/ravelarvargas/Downloads',
       file_name = 'test_overlap3',p_width = 10,p_height=10)
#####