#Script to download and analyse temporal DEGs
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')
cellage=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/1_cellage.csv')
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

cellage=cellage[cellage$gene_name%in%human_pc$external_gene_name,]
cellage$Senescence[cellage$Senescence=='Induces']='Induce CS'
cellage$Senescence[cellage$Senescence=='Inhibits']='Inhibit CS'
cellage$Senescence[cellage$Senescence=='Overexpressed']='Up in RS'
cellage$Senescence[cellage$Senescence=='Underexpressed']='Down in RS'

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
time_degs=compile_time_files(dir_use='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/',
                   deg_file='degs.csv')
save_csv(time_degs,file_name = 'all_time_degs',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

#GOI
inflammation_genes=c('IL1A','IL1B','IL6','CXCL8',
                     # 'TGFB1',
                     # 'NFKB1',
                     # 'FOS',
                     'TGFA')

cell_cycle_genes=c('CDKN1A','CDKN2A','CCNE1',
                   'CDK4',
                   # 'CDK6',
                   'CDK2',
                   # 'CCND1'
                   'MDM2'
                   # 'CDK1','CCNA2','CCNB1'
                   # 'MKI67'
)

architecture_genes=c('HMGA1','HMGA2',
                     'LMNB1',
                     # 'CBX5',
                     'SUV39H1',
                     'SUV39H2'
                     # 'EZH2'
                     # 'HP1','HIRA',
                     # 'HDAC1','HDAC2'
)

autophagy_lysosome=c(
  # 'LAMP1',
  'LAMP2',
  # 'MTOR',
  # 'ATG7',
  'ULK1',
  'DEPTOR',
  'PTEN'
  # 'LAMP2A',
  # 'CASTOR3'
  )

TFs=c('FOXO3',
      'FOXO1',
      'FOXO4',
      'NFKB1',
      'TP53',
      'ATF4',
      'HIF1A',
      'CEBPB')

apoptosis_genes=c(
  # 'BCL2',
                  # 'BAK1',
                  'BBC3', #PUMA
                  'PMAIP1', #NOXA
                  'BCL2L11',
                  'TNFRSF6B',
                  # 'CASP1',
                  'FAS')

time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = c(inflammation_genes),
                      gene_column = 'ensembl',
                      gene_label = 'Inflammation\nSASP')

time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = c(cell_cycle_genes),
                      gene_column = 'ensembl',
                      gene_label = 'Cell Cycle\nKI67')

time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = architecture_genes,
                      gene_column = 'ensembl',
                      gene_label = 'Chromatin Architecture\nSAHF')

time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = autophagy_lysosome,
                      gene_column = 'ensembl',
                      gene_label = 'Autophagy/Lysosome\nβ-gal staining')

time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = TFs,
                      gene_column = 'ensembl',
                      gene_label = 'TFs')

time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = apoptosis_genes,
                      gene_column = 'ensembl',
                      gene_label = 'Apoptosis')

time_degs=factor_column_and_modify(df = time_degs,
                                   column = 'group_1',
                                   old_list = rev(c('4_days','10_days','20_days')),
                                   keyword = NULL)

time_degs=factor_column_and_modify(df = time_degs,
                                   column = 'cell_type',
                                   old_list = c('fibroblast','keratinocyte',
                                                'melanocyte'),
                                   keyword = NULL)

suppressor_expression=compare_expression(degs=time_degs,
                                         gene_col = 'ensembl',
                                         group = 'group_1',
                                         genes=c(inflammation_genes,
                                                 cell_cycle_genes,
                                                 architecture_genes,
                                                 TFs,
                                                 autophagy_lysosome,
                                                 apoptosis_genes),
                                         facet = 'cell_type',
                                         facet2 = 'group',
                                         flip_axis = FALSE,
                                         scale_temp = 2,fill_cap = 6,
                                         tilt_x = TRUE,
                                         text_size = 13,
                                         legend_bottom = TRUE,
                                         condition_lab='Days Post-Irradiation'
)

save_p(suppressor_expression,file_name = 'temporal_GOI',
       save_dir='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/',
       p_height = 6.5,p_width = 13.2)

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
max(time_deg_overlap_sim$Freq[time_deg_overlap_sim$outgroup=='none_up'&
                                time_deg_overlap_sim$Var1==8])
max(time_deg_overlap_sim$Freq[time_deg_overlap_sim$outgroup=='none_down'&
                                time_deg_overlap_sim$Var1==9])

time_degs_sig=time_degs[time_degs$sig=='y',]

table(time_degs_sig%>%dplyr::select(group_1,direction_1,cell_type))

time_degs_sig$cell_type_time=paste0(time_degs_sig$cell_type,'_',
                                    time_degs_sig$group_1)
time_degs_sig_up=time_degs_sig[time_degs_sig$direction_1=='up',]
time_degs_sig_down=time_degs_sig[time_degs_sig$direction_1=='down',]
time_degs_sig_up_fibroblast=time_degs_sig_up[time_degs_sig_up$cell_type=='fibroblast',]
time_degs_sig_down_fibroblast=time_degs_sig_down[time_degs_sig_down$cell_type=='fibroblast',]
length(table(time_degs_sig_up_fibroblast$ensembl)[table(time_degs_sig_up_fibroblast$ensembl)==3])
length(table(time_degs_sig_down_fibroblast$ensembl)[table(time_degs_sig_down_fibroblast$ensembl)==3])
time_degs_sig_up_melanocyte=time_degs_sig_up[time_degs_sig_up$cell_type=='melanocyte',]
time_degs_sig_down_melanocyte=time_degs_sig_down[time_degs_sig_down$cell_type=='melanocyte',]
length(table(time_degs_sig_up_melanocyte$ensembl)[table(time_degs_sig_up_melanocyte$ensembl)==3])
length(table(time_degs_sig_down_melanocyte$ensembl)[table(time_degs_sig_down_melanocyte$ensembl)==3])
time_degs_sig_up_keratinocyte=time_degs_sig_up[time_degs_sig_up$cell_type=='keratinocyte',]
time_degs_sig_down_keratinocyte=time_degs_sig_down[time_degs_sig_down$cell_type=='keratinocyte',]
length(table(time_degs_sig_up_keratinocyte$ensembl)[table(time_degs_sig_up_keratinocyte$ensembl)==3])
length(table(time_degs_sig_down_keratinocyte$ensembl)[table(time_degs_sig_down_keratinocyte$ensembl)==3])

length(table(time_degs_sig_up$ensembl)[table(time_degs_sig_up$ensembl)==9])
message(paste0(names(table(time_degs_sig_up$ensembl)[table(time_degs_sig_up$ensembl)==9]),
       collapse=', '))
length(table(time_degs_sig_down$ensembl)[table(time_degs_sig_down$ensembl)==9])
message(paste0(names(table(time_degs_sig_down$ensembl)[table(time_degs_sig_down$ensembl)==9]),
               collapse=', '))
#####
#Plot DEG overlaps with each other
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
  temp_order = c('4 days','10 days','20 days'),
  graph_scale = 2
)

save_csv(time_deg_overlaps$df,file_name = 'temporal_comparison',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

save_p(plot = time_deg_overlaps$p,
       file_name = 'temporal_comparison',
       save_dir ='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 13,
       p_height = 11)

#Plot DEG overlaps with CellAge
db_deg_overlaps_time=overlap_function(df_1 = time_degs_sig,
                                      gene_col_1 = 'ensembl',
                                      group_col_1 = c('cell_type','group_1_dir'),
                                      df_2 = cellage,
                                      group_col_2 = 'Senescence',
                                      gene_col_2 = 'gene_name',
                                      background = time_degs$ensembl,
                                      carry_col_1 = c('group_1','direction_1'))

db_deg_overlaps_time$facet=paste0(db_deg_overlaps_time$cell_type,'\n',
                                  db_deg_overlaps_time$direction_1)

db_deg_overlaps_time$group_1=gsub(
  db_deg_overlaps_time$group_1,pattern='_',replacement=' '
)

compiled_time_db_overlaps=time_graph_deg_db(
  degs=db_deg_overlaps_time,
  group_1_deg_col='group_1',
  group_1_dir_col='facet',
  group_2_deg_col='Senescence',
  pval_col='adj',
  odds_col='odds',
  n_col='actual',
  order = c('4 days', '10 days','20 days'),
  facet_col='Senescence',
  facet_order = c('Induce CS','Inhibit CS','Up in RS','Down in RS'),
  nudge_y=45,
  label_ns=FALSE,
  xlab='Days Post-Irradiation',
  ylab='Number of overlapping DEGs',tilt_x = TRUE)

save_p(compiled_time_db_overlaps,file_name = 'temporal_vs_cellage',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',p_width = 8)

save_csv(db_deg_overlaps_time,file_name = 'temporal_vs_cellage',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

#Overlap with stress responses
autophagy_pc=read.csv(paste0(save_dir_csv,'autophagy_genes.csv'))

autophagy_recount3_temporal=overlap_function(df_1=time_degs_sig,
                                    df_2=autophagy_pc,
                                    gene_col_1='ensembl',
                                    gene_col_2='Official.Gene.symbol',
                                    group_col_1=c('group_1','direction_1','cell_type'),
                                    group_col_2=c('func'),
                                    background=time_degs$ensembl)

autophagy_recount3_temporal=factor_column_and_modify(df = autophagy_recount3_temporal,
                                                   column = 'direction_1',
                                                   old_list = c('up','down'),
                                                   keyword = 'in Arrest')

autophagy_recount3_temporal=factor_column_and_modify(df = autophagy_recount3_temporal,
                                                     column = 'group_1',
                                                     old_list = c('4_days',
                                                                  '10_days',
                                                                  '20_days'),
                                                     keyword = '\nPost-Irradiation')

auto_overlap=create_overlap_plot(deg_db_overlap = autophagy_recount3_temporal,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'direction_1',facet_2 = 'cell_type',
                    x = 'group_1',
                    y = 'func',
                    xlab = 'Temporal DEGs',
                    ylab = 'Autophagy Regulators',
                    ggtitle = 'Arrest DEGs vs Autophagy',x_tilt = TRUE)

save_p(auto_overlap,file_name = 'auto_overlap_temporal',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 6.5,p_height = 4)
save_csv(autophagy_recount3_temporal,'auto_overlap_temporal',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

lyso=read.csv(paste0(save_dir_csv,'lyso_genes.csv'))

lyso_recount3_temporal=overlap_function(df_1=time_degs_sig,
                                             df_2=lyso,
                                             gene_col_1='ensembl',
                                             gene_col_2='Official.Gene.symbol',
                                             group_col_1=c('group_1','direction_1',
                                                           'cell_type'),
                                             group_col_2=c('func'),
                                             background=time_degs$ensembl)

lyso_recount3_temporal=factor_column_and_modify(df = lyso_recount3_temporal,
                                                     column = 'direction_1',
                                                     old_list = c('down','up'),
                                                     keyword = 'in Arrest')


lyso_recount3_temporal=factor_column_and_modify(df = lyso_recount3_temporal,
                                                     column = 'group_1',
                                                     old_list = c('4_days',
                                                                  '10_days',
                                                                  '20_days'),
                                                     keyword = '\nPost-Irradiation')

lyso_temporal=create_overlap_plot(deg_db_overlap = lyso_recount3_temporal,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'cell_type',
                    x = 'group_1',
                    y = 'direction_1',
                    xlab = 'Temporal DEGs',
                    ylab = 'Lysosomal Genes',
                    ggtitle = 'Arrest DEGs vs Lysosomal Regulators',x_tilt = TRUE)

save_p(lyso_temporal,file_name = 'lysosome_overlap_temporal',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 6,p_height = 2.5)
save_csv(lyso_recount3_temporal,'lysosome_overlap_temporal',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

inflammation_genes=read.csv(paste0(save_dir_csv,
                                   'go_inflamm.csv'))

inflamm_recount3_temporal=overlap_function(df_1=time_degs_sig,
                                        df_2=inflammation_genes,
                                        gene_col_1='ensembl',
                                        gene_col_2='external_gene_name',
                                        group_col_1=c('group_1','direction_1',
                                                      'cell_type'),
                                        group_col_2=c('go_term'),
                                        background=time_degs$ensembl)

inflamm_recount3_temporal=factor_column_and_modify(df = inflamm_recount3_temporal,
                                       column = 'direction_1',
                                       old_list = c('up','down'),
                                       keyword = 'in Arrest')

inflamm_recount3_temporal=factor_column_and_modify(df = inflamm_recount3_temporal,
                                                   column = 'group_1',
                                                   old_list = c('4_days',
                                                                '10_days',
                                                                '20_days'),
                                                   keyword = '\nPost-Irradiation')

inflammation_recount3_p=create_overlap_plot(deg_db_overlap = inflamm_recount3_temporal,
                                            odds_column = 'odds',
                                            pval_col = 'adj',
                                            facet_col = 'direction_1',
                                            facet_2 = 'cell_type',
                                            x = 'group_1',
                                            y = 'go_term',
                                            xlab = 'Temporal DEGs',
                                            ylab = 'Inflammation GO Term Genes',
                                            ggtitle = 'Arrest DEGs vs Inflammation',
                                            x_tilt = TRUE)

save_p(inflammation_recount3_p,file_name = 'inflamm_overlap_temmporal',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 8,p_height = 4)
save_csv(inflamm_recount3_temporal,'inflamm_overlap_temmporal',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')

all_sasp_fixed=read.csv(paste0(save_dir_csv,'sasp_atlas_genes.csv'))

all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_ras']='Fibroblast RAS'
all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_irradiated']='Fibroblast Irradiated'
all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_atazanivir']='Fibroblast Atazanavir'
all_sasp_fixed$group[all_sasp_fixed$group=='epithelial_irradiated']='Epithelial Irradiated'

sasp_recount3=overlap_function(df_1=time_degs_sig,
                               df_2=all_sasp_fixed,
                               gene_col_1='ensembl',
                               gene_col_2='genes',
                               group_col_1=c('group_1','direction_1','cell_type'),
                               group_col_2=c('group','arrest'),
                               background=time_degs$ensembl)

sasp_recount3$cell_dir=paste0(sasp_recount3$cell_type,'_',
                              sasp_recount3$direction_1)
sasp_recount3=factor_column_and_modify(df = sasp_recount3,
                                       column = 'group_1',
                                       old_list = c('4_days',
                                                    '10_days',
                                                    '20_days'),
                                       keyword = '\nPost-Irradiation')

sasp_recount3=factor_column_and_modify(df = sasp_recount3,
                                       column = 'cell_dir',
                                       old_list = c('fibroblast_up',
                                                    'fibroblast_down',
                                                    'keratinocyte_up',
                                                    'keratinocyte_down',
                                                    'melanocyte_up',
                                                    'melanocyte_down'),
                                       keyword = NULL)
sasp_atlas_recount3_p=create_overlap_plot(deg_db_overlap = sasp_recount3,
                                            odds_column = 'odds',
                                            pval_col = 'adj',
                                            facet_col = 'arrest',
                                            facet_2 = 'cell_dir',
                                            x = 'group_1',
                                            y = 'group',
                                            xlab = 'Temporal DEGs',
                                            ylab = 'SASP Atlas Profiles',
                                            ggtitle = 'Arrest DEGs vs SASP Atlas',
                                            x_tilt = TRUE)

save_p(sasp_atlas_recount3_p,file_name = 'sasp_overlap_temporal',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140',
       p_width = 12,p_height = 3.8)
save_csv(sasp_recount3,'sasp_overlap_temporal',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140')
