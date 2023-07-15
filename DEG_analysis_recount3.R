#DEG analysis
cs_degs=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/cs_degs.csv')
cq_degs=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/cq_degs.csv')
cellage=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/1_cellage.csv')
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')
#####
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)

arrest_degs_merged=rbind(cs_degs,
                         cq_degs)
arrest_degs_merged=arrest_degs_merged[arrest_degs_merged$gene%in%
                                        human_pc$external_gene_name,]

cellage=cellage[cellage$gene_name%in%
                  human_pc$external_gene_name,]
save_csv(cellage,file_name = '1_cellage',path = save_dir_csv)
arrest_degs_merged$arrest=ifelse(grepl(arrest_degs_merged$group_1,pattern = 'CS'),'CS','CQ')
#####
#GOI
inflammation_goi=c('IL1A','IL1B','IL6','CXCL8','TGFB1',
                     # 'NFKB1',
                     # 'FOS',
                     'TGFA')

cell_cycle_genes=c('CDKN1A','CDKN2A','CCNE1',
                   'CDK4','CDK6','CDK2','CCND1',
                   'MDM2'
                   #'CDK1','CCNA2','CCNB1'
                   # 'MKI67'
)

architecture_genes=c('HMGA1','HMGA2',
                     # 'LMNB1',
                     # 'HP1',
                     'SUV39H1',
                     'SUV39H2'
                     # 'EZH2'
                     # 'HIRA','HDAC1','HDAC2'
)

apoptosis_genes=c('BCL2',
                  'BAK1',
                  'BCL2L1',
                  'BBC3', #PUMA
                  'PMAIP1', #NOXA
                  # 'BCL2L11',
                  'TNFRSF6B',
                  # 'TNFRSF10B',
                  # 'CASP1',
                  'FAS'
                  # 'BIRC5'
                  )

autophagy_lysosome=c(
  # 'LAMP1',
                     'LAMP2',
                     # 'MTOR',
                     # 'ATG7',
                     'ULK1',
                     'LAMP2A',
                     'DEPTOR',
                     'CASTOR3'
                     )

TFs=c('FOXO3',
      'FOXO4',
      # 'FOXO1',
      # 'NFKB1',
      'TP53',
      'ATF4',
      # 'HIF1A',
      'CEBPB')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(inflammation_goi),
                               gene_column = 'gene',
                               gene_label = 'Inflammation\nSASP')

arrest_degs_merged=label_genes(arrest_degs_merged,
            new_column = 'group',
            gene_vector = c(apoptosis_genes),
            gene_column = 'gene',
            gene_label = 'Apoptosis')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(cell_cycle_genes),
                               gene_column = 'gene',
                               gene_label = 'Cell Cycle\nKI67')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = architecture_genes,
                               gene_column = 'gene',
                               gene_label = 'Chromatin Architecture\nSAHF')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = autophagy_lysosome,
                               gene_column = 'gene',
                               gene_label = 'Autophagy/Lysosome\nβ-gal staining')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = TFs,
                               gene_column = 'gene',
                               gene_label = 'TFs')

arrest_degs_merged$group_1[arrest_degs_merged$group_1=='Contact_inhibited_CQ']='Contact-inhibited CQ'
arrest_degs_merged$group_1[arrest_degs_merged$group_1=='Serum_starved_CQ']='Serum-starved CQ'
arrest_degs_merged$group_1[arrest_degs_merged$group_1=='Replicative_CS']='RS'
arrest_degs_merged$group_1[arrest_degs_merged$group_1=='Stress_induced_CS']='SIPS'
arrest_degs_merged$group_1[arrest_degs_merged$group_1=='Oncogene_induced_CS']='OIS'
arrest_degs_merged$group_1=factor(arrest_degs_merged$group_1,
                                  levels=c('Contact-inhibited CQ',
                                           'Serum-starved CQ',
                                           'RS',
                                           'SIPS',
                                           'OIS'))
suppressor_expression=compare_expression(degs=arrest_degs_merged,
                                         gene_col = 'gene',
                                         group = 'group_1',
                                         genes=c(inflammation_goi,
                                                 cell_cycle_genes,
                                                 architecture_genes,
                                                 TFs,
                                                 autophagy_lysosome,
                                                 apoptosis_genes),
                                         facet = 'group',
                                         # facet2 = 'arrest',
                                         flip_axis = TRUE,
                                         scale_temp = 2,
                                         fill_cap = 6,
                                         tilt_x = TRUE
)

save_p(suppressor_expression,
       file_name = 'goi_expression',save_dir = save_dir_figure,p_width = 7.5,
       p_height = 4)
#####
#Upset plots
arrest_degs_merged_sig=arrest_degs_merged[arrest_degs_merged$sig=='y',]

arrest_upset_up=upsetR_plot_intersection(
  arrest_degs = arrest_degs_merged_sig,
            dir_col = 'group_1_dir',
            gene_col = 'gene',
            dir = 'up',
            min_intersect_size = 1)

save_p(plot = arrest_upset_up,
       file_name = 'recount_upset_up',
       save_dir = save_dir_figure_si,p_width=18)

arrest_upset_down=upsetR_plot_intersection(
  arrest_degs = arrest_degs_merged_sig,
  dir_col = 'group_1_dir',
  gene_col = 'gene',
  dir = 'down',
  min_intersect_size = 1)

save_p(plot = arrest_upset_down,
       file_name = 'recount_upset_down',
       save_dir = save_dir_figure_si,p_width=18)

arrest_degs_merged_sim=arrest_degs_merged%>%dplyr::select(
  gene,group_1_dir,group_2_dir,accession,sig
)

sim_overlaps=simulate_overlaps(relevant_degs=arrest_degs_merged,
                  facet_col='accession',gene_col='gene',
                  simulation_n=10000)

save_csv(sim_overlaps,
         file_name = 'recount_deg_sim',
         path = save_dir_csv)

#simulation cumsum
sim_overlaps$outgroup[sim_overlaps$outgroup=='Proliferating_up']='Down vs Proliferating'
sim_overlaps$outgroup[sim_overlaps$outgroup=='Proliferating_down']='Up vs Proliferating'
sim_cumsum=plot_top_simulation_cumsum(sim_overlaps,plot_cutoff = FALSE)
save_csv(sim_cumsum$cumsum,file_name = 'reverse_cumsum_recount',
       path = save_dir_csv)

save_p(sim_cumsum$p,file_name = 'reverse_cumsum_recount',
       save_dir  = save_dir_figure_si)

#####
# Overlap with CellAge

arrest_degs_merged$accession=
  arrest_degs_merged$group_1
arrest_degs_merged$accession=gsub(arrest_degs_merged$accession,pattern='_',
                                  replacement=' ')
arrest_degs_merged$dir_accession=ifelse(grepl(arrest_degs_merged$dir_accession,
                                              pattern='Proliferating_down'),'Up in Arrest',
                                        'Down in Arrest')
  # factor(arrest_degs$accession,
  #        levels=c('Contact_inhibited_CQ vs\nProliferating',
  #                 'Serum_starved_CQ vs\nProliferating',
  #                 'Replicative_CS vs\nProliferating',
  #                 'Stress_induced_CS vs\nProliferating',
  #                 'Oncogene_induced_CS vs\nProliferating'))

arrest_degs_merged_sig$group_1=
  as.character(arrest_degs_merged_sig$group_1)
deg_cellage_overlap=overlap_function(df_1=arrest_degs_merged_sig,
                                     df_2=cellage,
                                     gene_col_1='gene',
                                     gene_col_2='gene_name',
                                     group_col_1=c('group_1','direction_1'),
                                     group_col_2=c('Senescence'),
                                     background=arrest_degs_merged$gene)
save_csv(data = deg_cellage_overlap,
         file_name = 'recount_cellage_overlap',
         path = save_dir_csv)

deg_cellage_overlap$direction_1=gsub(deg_cellage_overlap$direction_1,
                                     pattern='up',replacement='Up in\nArrest')
deg_cellage_overlap$direction_1=gsub(deg_cellage_overlap$direction_1,
                                     pattern='down',replacement='Down in\nArrest')
deg_cellage_overlap$group_1=gsub(deg_cellage_overlap$group_1,pattern='_',replacement=' ')
deg_cellage_overlap$Senescence=gsub(deg_cellage_overlap$Senescence,
                                    pattern='Underexpressed',replacement='Underexpressed in RS')
deg_cellage_overlap$Senescence=gsub(deg_cellage_overlap$Senescence,
                                    pattern='Overexpressed',replacement='Overexpressed in RS')
deg_cellage_overlap$Senescence=gsub(deg_cellage_overlap$Senescence,
                                    pattern='Inhibits',replacement='Inhibit CS')
deg_cellage_overlap$Senescence=gsub(deg_cellage_overlap$Senescence,
                                    pattern='Induces',replacement='Induce CS')
recount_cellage_overlap_plot=create_overlap_plot(deg_db_overlap = deg_cellage_overlap,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'group_1',
                    x = 'direction_1',
                    y = 'Senescence',
                    xlab = 'Cell Cycle Arrest DEGs',
                    ylab = 'CellAge',text_size = 15)
save_p(plot = recount_cellage_overlap_plot,
       file_name = 'recount_cellage_overlap',
       save_dir = save_dir_figure,p_width = 10)
#####
#self overlaps
self_overlaps_recount_to_csv=overlap_within_df(dataframe=arrest_degs_merged_sig,
                  group_col='group_1',
                  other_cols='direction_1',
                  gene_col='gene',
                  background=arrest_degs_merged$gene,
                  remove_self = TRUE
)

save_csv(self_overlaps_recount_to_csv,
         file_name = 'self_overlap_recount',
         path = save_dir_csv)

get_pval=self_overlaps_recount_to_csv[,colnames(self_overlaps_recount_to_csv)==c('accession','adj')]

self_overlaps_recount_to_plot=overlap_within_df(dataframe=arrest_degs_merged_sig,
                                               group_col='group_1',
                                               other_cols='direction_1',
                                               gene_col='gene',
                                               background=arrest_degs_merged$gene,
                                               remove_self = FALSE
)
self_overlaps_recount_to_plot=self_overlaps_recount_to_plot[,colnames(self_overlaps_recount_to_plot)!='adj']
self_overlaps_recount_to_plot=merge(self_overlaps_recount_to_plot,get_pval)

self_overlaps_recount_to_plot$outgroup_dir=paste0(self_overlaps_recount_to_plot$group_1_2,' ',self_overlaps_recount_to_plot$direction_1_2)
self_overlaps_recount_to_plot$outgroup_dir=gsub(self_overlaps_recount_to_plot$outgroup_dir,pattern='_',replacement=' ')
self_overlaps_recount_to_plot$direction_1_1[self_overlaps_recount_to_plot$direction_1_1=='up']=
  'Up with Arrest'
self_overlaps_recount_to_plot$direction_1_1[self_overlaps_recount_to_plot$direction_1_1=='down']=
  'Down with Arrest'

self_overlap_p=plot_self_overlaps(self_overlaps_recount_to_plot = self_overlaps_recount_to_plot,
                   odds_column = 'odds',
                   x = 'direction_1_1',
                   y = 'outgroup_dir',
                   facet_col = 'group_1_1',
                   xlab = 'Facet Direction with Arrest',
                   ylab = 'Y-Axis Direction with Arrest',
                   text_size = 15,
                   scale=3,angle_x = TRUE)
save_p(self_overlap_p,file_name = 'self_overlap_recount',
       save_dir = save_dir_figure,p_width=16,p_height = 9)
#table DEGs
DEG_count=data.frame(table(arrest_degs_merged_sig%>%dplyr::select(group_1,direction_1)))
save_csv(DEG_count,file_name = 'deg_count',path = save_dir_csv)

#####
#volcano plot
arrest_degs_merged$group_1=as.character(arrest_degs_merged$group_1)
arrest_degs_merged_cs=arrest_degs_merged[arrest_degs_merged$arrest=='CS',]
volcano_plot_cs=volcano_function(degs=arrest_degs_merged_cs,
                              custom_labs='group_1',
                              gene_col='gene')

save_p(plot = volcano_plot_cs,file_name = 'volcano_plot_cs',save_dir = save_dir_figure)

arrest_degs_merged_cq=arrest_degs_merged[arrest_degs_merged$arrest=='CQ',]
volcano_plot_cq=volcano_function(degs=arrest_degs_merged_cq,
                                 custom_labs='group_1',
                                 gene_col='gene')
save_p(plot = volcano_plot_cq,file_name = 'volcano_plot_cq',save_dir = save_dir_figure)

#indivudal plots
arrest_degs_merged_rs=arrest_degs_merged_cs[arrest_degs_merged_cs$group_1=='RS',]
volcano_plot_rs=volcano_function(degs=arrest_degs_merged_rs,
                                 gene_col='gene',title = 'RS')
# save_p(plot = volcano_plot_rs,file_name = 'volcano_plot_rs',save_dir = save_dir_figure)
arrest_degs_merged_sips=arrest_degs_merged_cs[arrest_degs_merged_cs$group_1=='SIPS',]
volcano_plot_sips=volcano_function(degs=arrest_degs_merged_sips,
                                 gene_col='gene',title='SIPS')
# save_p(plot = volcano_plot_sips,file_name = 'volcano_plot_sips',save_dir = save_dir_figure)
arrest_degs_merged_ois=arrest_degs_merged_cs[arrest_degs_merged_cs$group_1=='OIS',]
volcano_plot_ois=volcano_function(degs=arrest_degs_merged_ois,
                                   gene_col='gene',title='OIS')
# save_p(plot = volcano_plot_ois,file_name = 'volcano_plot_ois',save_dir = save_dir_figure)
save_p(plot=volcano_plot_rs/volcano_plot_sips/volcano_plot_ois,
       file_name = 'volcano_all_cs',save_dir = save_dir_figure,p_height=15)

arrest_degs_merged_contact=arrest_degs_merged_cs[arrest_degs_merged_cq$group_1=='Contact-inhibited CQ',]
volcano_plot_contact=volcano_function(degs=arrest_degs_merged_contact,
                                 gene_col='gene',title = 'Contact-inhibited CQ')

arrest_degs_merged_serum=arrest_degs_merged_cs[arrest_degs_merged_cq$group_1=='Serum-starved CQ',]
volcano_plot_serum=volcano_function(degs=arrest_degs_merged_serum,
                                      gene_col='gene',title = 'Serum-starved CQ')
save_p(plot=volcano_plot_contact/volcano_plot_serum,
       file_name = 'volcano_all_cq',save_dir = save_dir_figure,p_height=10)
#####
#enrichment
human_pc_entrez=getBM(attributes=c('external_gene_name','entrezgene_id'),
                      filters = 'biotype',
                      values = c('protein_coding'),
                      mart = ensembl100)

arrest_degs_merged_up=arrest_degs_merged[arrest_degs_merged$direction_1=='up'&
                                           arrest_degs_merged$sig=='y',]
common_up=find_common_genes(df = arrest_degs_merged_up,
                  group_col = 'group_1_dir',gene_col = 'gene')
common_up_df=data.frame(gene=common_up)
common_up_df$dir='up'

arrest_degs_merged_up=arrest_degs_merged[arrest_degs_merged$direction_1=='up',]
arrest_degs_merged_up_table=table(arrest_degs_merged_up$gene)
arrest_degs_merged_up_background=names(arrest_degs_merged_up_table)[arrest_degs_merged_up_table==5]
arrest_degs_merged_down=arrest_degs_merged[arrest_degs_merged$direction_1=='down',]
arrest_degs_merged_down_table=table(arrest_degs_merged_down$gene)
arrest_degs_merged_down_background=names(arrest_degs_merged_down_table)[arrest_degs_merged_down_table==5]

enrich_up=enrich_genes(gene_list = common_up,
             background = arrest_degs_merged_up_background,
             gene_dictionary = human_pc_entrez,
             use_ensembl = FALSE)

arrest_degs_merged_up=arrest_degs_merged[arrest_degs_merged$direction_1=='up'&
                                           arrest_degs_merged$sig=='y',]
common_up=find_common_genes(df = arrest_degs_merged_up,
                  group_col = 'group_1_dir',gene_col = 'gene')

enrich_up=enrich_genes(gene_list = common_up,
             background = unique(arrest_degs_merged$gene[arrest_degs_merged$direction_1==
                                                           'up']),
             gene_dictionary = human_pc_entrez,
             use_ensembl = FALSE)

shared_up_go=enrichment_dotplot(enrichment_table = 
                     enrich_up$enrichment[enrich_up$enrichment$enrichment=='GO',],
                   enrichment_type = 'GO')

save_p(plot = shared_up_go,file_name = 'shared_up_go',
       save_dir = save_dir_figure,p_width = 11,p_height = 4)

arrest_degs_merged_down=arrest_degs_merged[arrest_degs_merged$direction_1=='down'&
                                           arrest_degs_merged$sig=='y',]
common_down=find_common_genes(df = arrest_degs_merged_down,
                            group_col = 'group_1_dir',gene_col = 'gene')
common_down_df=data.frame(gene=common_down)
common_down_df$dir='down'

enrich_down=enrich_genes(gene_list = common_down,
                       background = arrest_degs_merged_down_background,
                       gene_dictionary = human_pc_entrez,
                       use_ensembl = FALSE)

save_csv(enrich_down$enrichment,file_name = 'common_deg_enrichment',
         path = save_dir_csv)

enrich_down_kegg=enrichment_dotplot(enrich_down$enrichment[enrich_down$enrichment$enrichment=='KEGG',],
                                    text_size = 18)

save_p(enrich_down_kegg,file_name = 'shared_kegg_down',
       save_dir = save_dir_figure,p_width = 10)

enrich_down_go=enrich_down$enrichment[enrich_down$enrichment$enrichment=='GO',]

sem_sim=get_sem_sim()

rrvgo_treemap(ont_sem_sim = sem_sim,
              go_terms = enrich_down_go,
              file_name = 'enriched_GO_down',
              save_dir = save_dir_figure)

common_all=rbind(common_down_df,common_up_df)
save_csv(common_all,file_name = 'common_intercept',path = save_dir_csv)

common_background_down=data.frame(gene=arrest_degs_merged_down_background)
common_background_down$dir='down'
common_background_up=data.frame(gene=arrest_degs_merged_up_background)
common_background_up$dir='up'
common_background_all=rbind(common_background_down,common_background_up)
save_csv(common_background_all,file_name = 'common_background',path = save_dir_csv)
#####
#Overlap with autophagy, lysosome, inflammation, and the SASP atlas
arrest_degs_merged_sig=factor_column_and_modify(df = arrest_degs_merged_sig,
                                            column = 'direction_1',
                                            old_list = c('down','up'),
                                            keyword = '\nin Arrest')

arrest_degs_merged_sig$direction_1=as.character(arrest_degs_merged_sig$direction_1)

autophagy=read.csv('/Users/ravelarvargas/Downloads/post_thesis/autophagy_gene.csv')

autophagy_pc=autophagy[autophagy$Official.Gene.symbol%in%human_pc$external_gene_name,]
autophagy_pc$func=NA
autophagy_pc$func[grepl(autophagy_pc$Biological.Function,pattern = 'Positive')]='Positive'
autophagy_pc$func[grepl(autophagy_pc$Biological.Function,pattern = 'Negative')]='Negative'
auto_lyso=autophagy_pc[autophagy_pc$Group=='lysosome',]
auto_lyso$func='lyso'
save_csv(auto_lyso,file_name = 'lyso_genes',path = save_dir_csv)
autophagy_pc=autophagy_pc[!is.na(autophagy_pc$func),]

save_csv(autophagy_pc,file_name='autophagy_genes',path=save_dir_csv)
autophagy_recount3=overlap_function(df_1=arrest_degs_merged_sig,
                                     df_2=autophagy_pc,
                                     gene_col_1='gene',
                                     gene_col_2='Official.Gene.symbol',
                                     group_col_1=c('group_1','direction_1'),
                                     group_col_2=c('func'),
                                     background=arrest_degs_merged$gene)
save_csv(autophagy_recount3,file_name = 'auto_recount_overlap',path = save_dir_csv)

autophagy_recount3_p=create_overlap_plot(deg_db_overlap = autophagy_recount3,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'group_1',
                    x = 'direction_1',
                    y = 'func',
                    xlab = 'Cell Cycle Arrest DEGs',
                    ylab = 'Autophagy Regulators',
                    ggtitle = 'Arrest DEGs vs Autophagy')
save_p(autophagy_recount3_p,file_name = 'auto_recount_overlap',
       save_dir = save_dir_figure,p_height = 2.5,p_width = 6)

#Lysosome
recount3_lyso=overlap_function(df_1=arrest_degs_merged_sig,
                 df_2 = auto_lyso,
                 gene_col_1='gene',
                 gene_col_2='Official.Gene.symbol',
                 group_col_1=c('group_1','direction_1'),
                 group_col_2=c('func'),
                 background=arrest_degs_merged$gene)

lyso_recount3_p=create_overlap_plot(deg_db_overlap = recount3_lyso,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'group_1',
                    x = 'direction_1',
                    y = 'func',
                    xlab = 'Cell Cycle Arrest DEGs',
                    ylab = NULL,remove_y = TRUE,
                    ggtitle = 'Arrest DEGs vs Lysosome Genes')

save_csv(recount3_lyso,file_name = 'recount_lyso_overlap',path = save_dir_csv)
save_p(lyso_recount3_p,file_name = 'recount_lyso_overlap',
       save_dir = save_dir_figure,p_height = 2.5,p_width = 5)

#Inflammation
inflammation_genes=getBM(attributes=c('external_gene_name','go_id'),
                         filters = 'go',
                         values = c('GO:0050729','GO:0050728'),
                         mart = ensembl100)

inflammation_genes=inflammation_genes[inflammation_genes$external_gene_name%in%human_pc$external_gene_name,]

inflammation_genes=inflammation_genes[inflammation_genes$go_id%in%
                                        c('GO:0050729','GO:0050728'),]

inflammation_genes$go_term=ifelse(inflammation_genes$go_id=='GO:0050729',
                                  'Promotes inflammatory response',
                                  'Inhibits inflammatory response')

save_csv(data = inflammation_genes,file_name = 'go_inflamm',path = save_dir_csv)

inflam_recount3=overlap_function(df_1=arrest_degs_merged_sig,
                                    df_2=inflammation_genes,
                                    gene_col_1='gene',
                                    gene_col_2='external_gene_name',
                                    group_col_1=c('group_1','direction_1'),
                                    group_col_2=c('go_term'),
                                    background=arrest_degs_merged$gene)

inflammation_recount3_p=create_overlap_plot(deg_db_overlap = inflam_recount3,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'group_1',
                    x = 'direction_1',
                    y = 'go_term',
                    xlab = 'Cell Cycle Arrest DEGs',
                    ylab = 'Inflammation GO Term Genes',
                    ggtitle = 'Arrest DEGs vs Inflammation')

save_csv(inflam_recount3,file_name = 'recount_inflamm_overlap',path = save_dir_csv)
save_p(inflammation_recount3_p,file_name = 'recount_inflamm_overlap',
       save_dir = save_dir_figure,p_height = 3)

#SASP
# sasp_dir='/Users/ravelarvargas/Downloads/post_thesis/sasp/'
# 
# sasp_files=list.files(sasp_dir)
# sasp_files=sasp_files[grepl(sasp_files,pattern = '.csv')]
# 
# colnames_temp=c('comparison','eg','genes','go_bp','go_cel','go_mol',
#                 'log2ratio','protein_description','pval','qval','sd',
#                 'uniprot','num_ratio','num_total_eg','num_unique_pep','log2qval')
# 
# all_sasp=c()
# for(i in sasp_files){
#   temp_sasp=read.csv(paste0(sasp_dir,i))
#   temp_sasp=temp_sasp%>%dplyr::select(sort(names(.)))
#   colnames(temp_sasp)=colnames_temp
#   temp_sasp$group=gsub(i,pattern = '.csv',replacement = '')
#   # if(i=='epithelial_irradiated.csv'){
#   #   temp_sasp$uniprot=temp_sasp$genes
#   # }
#   all_sasp=rbind(temp_sasp,all_sasp)
# }
# 
# all_sasp=all_sasp[,colnames(all_sasp)%in%c('comparison','genes','log2ratio','pval','qval','group')]
# all_sasp_fixed=c()
# for(i in 1:nrow(all_sasp)){
#   if(!grepl(all_sasp[i,][['genes']],pattern = ';')){
#     all_sasp_fixed=rbind(all_sasp_fixed,all_sasp[i,])
#   }else{
#     temp_protein=unlist(strsplit(all_sasp[i,][['genes']],split = ';'))
#     temp_merge=merge(all_sasp[i,],temp_protein)
#     temp_merge$genes=temp_merge$y
#     temp_merge=temp_merge[,colnames(temp_merge)!='y']
#     all_sasp_fixed=rbind(all_sasp_fixed,temp_merge)
#   }
# }
# all_sasp_fixed=all_sasp_fixed[all_sasp_fixed$genes%in%human_pc$external_gene_name,]
# 
# all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_ras']='Fibroblast RAS'
# all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_irradiated']='Fibroblast Irradiated'
# all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_atazanavir']='Fibroblast Atazanavir'
# all_sasp_fixed$group[all_sasp_fixed$group=='epithelial_irradiated']='Epithelial Irradiated'
# 
# all_sasp_fixed$arrest=ifelse(all_sasp_fixed$log2ratio>0,'CS Secretion','CQ Secretion')
# 
# all_sasp_fixed$group_arrest=paste0(all_sasp_fixed$group,' ',all_sasp_fixed$arrest)

# save_csv(all_sasp_fixed,file_name = 'sasp_atlas_genes',path = save_dir_csv)
all_sasp_fixed=read.csv(paste0(save_dir_csv,'sasp_atlas_genes.csv'))

all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_ras']='Fibroblast RAS'
all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_irradiated']='Fibroblast Irradiated'
all_sasp_fixed$group[all_sasp_fixed$group=='fibroblast_atazanivir']='Fibroblast Atazanavir'
all_sasp_fixed$group[all_sasp_fixed$group=='epithelial_irradiated']='Epithelial Irradiated'

sasp_recount3=overlap_function(df_1=arrest_degs_merged_sig,
                                 df_2=all_sasp_fixed,
                                 gene_col_1='gene',
                                 gene_col_2='genes',
                                 group_col_1=c('group_1','direction_1'),
                                 group_col_2=c('arrest','group'),
                                 background=arrest_degs_merged$gene)

sasp_recount3$group_1=factor(sasp_recount3$group_1,
                             levels=c('Contact-inhibited CQ',
                                      'Serum-starved CQ',
                                      'RS',
                                      'SIPS',
                                      'OIS'))
sasp_recount3_p=create_overlap_plot(deg_db_overlap = sasp_recount3,
                                            odds_column = 'odds',
                                            pval_col = 'adj',
                                            facet_col = 'arrest',
                                    facet_2 = 'group_1',
                                            x = 'direction_1',
                                            y = 'group',
                                            xlab = 'Cell Cycle Arrest DEGs',
                                            ylab = 'SASP Atlas Profiles',
                                            ggtitle = 'Arrest DEGs vs SASP Atlas')

save_csv(sasp_recount3,file_name = 'sasp_recount3_overlap',path = save_dir_csv)
save_p(sasp_recount3_p,file_name = 'sasp_recount3_overlap',save_dir = save_dir_figure,p_width = 9,
       p_height = 3.5)

#Overlap atlas with inflammation
sasp_inflam=overlap_function(df_1=inflammation_genes,
                               df_2=all_sasp_fixed,
                               gene_col_1='external_gene_name',
                               gene_col_2='genes',
                               group_col_2=c('arrest','group'),
                               group_col_1=c('go_term'),
                               background=arrest_degs_merged$gene)

inflam_vs_sasp=create_overlap_plot(deg_db_overlap = sasp_inflam,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'arrest',
                    # facet_2 = 'group_1',
                    x = 'group',
                    y = 'go_term',
                    xlab = 'SASP Atlas Profiles',
                    ylab = 'Inflammation GO Terms',
                    ggtitle = 'Inflammation vs SASP Atlas',x_tilt = TRUE)
save_p(inflam_vs_sasp,
       file_name = 'inflam_vs_sasp',
       save_dir = save_dir_figure_si,p_width = 6,p_height = 3)
save_csv(sasp_inflam,file_name = 'inflam_sasp',path = save_dir_csv)
#####
compare_logfc_outgroup(degs=arrest_degs_merged,
                       ingroup_col = 'group_1',
                       outgroup_col = 'group_2',
                       degs_gene_col = 'gene',
                       outgroup='Proliferating',
                       # xlab='Day',
                       # facet_order=levels(time_degs$group_1),
                       database=cellage,
                       database_gene_col = 'gene_name',
                       database_group_col='Senescence')
compare_logfc_outgroup(degs=arrest_degs_merged,
                       ingroup_col = 'group_1',
                       outgroup_col = 'group_2',
                       degs_gene_col = 'gene',
                       outgroup='Proliferating',
                       # xlab='Day',
                       # facet_order=levels(time_degs$group_1),
                       database=common_background_all,
                       database_gene_col = 'gene',
                       database_group_col='dir')

#####