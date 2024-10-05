#DEG analysis
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')

ensembl100=useMart(host='https://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)

arrest_degs_merged=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/arrest_degs_final.csv')
arrest_degs_merged=arrest_degs_merged[arrest_degs_merged$gene%in%
                                        human_pc$external_gene_name,]

cs_gene_list=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/1_cellage.csv')

cs_gene_list$accession[cs_gene_list$accession=='Saul et al ']='SenMayo'
cs_gene_list$accession[cs_gene_list$accession=='Driver Inhibits']='CellAge Inhibits CS'
cs_gene_list$accession[cs_gene_list$accession=='Driver Induces']='CellAge Induces CS'
cs_gene_list$accession[cs_gene_list$accession=='Signature Underexpressed']='CellAge Down in RS'
cs_gene_list$accession[cs_gene_list$accession=='Signature Overexpressed']='CellAge Up in RS'
cs_gene_list$accession=
  gsub(cs_gene_list$accession,pattern='down',replacement='Down')
cs_gene_list$accession=
  gsub(cs_gene_list$accession,pattern='up',replacement='Up')

#remove TAO for now
cs_gene_list=cs_gene_list[!grepl(cs_gene_list$source,pattern='Tao'),]
cs_factor_level=c('CellAge Up in RS',
                  'CellAge Down in RS',
                  'CellAge Induces CS',
                  'CellAge Inhibits CS',
                  'SenMayo',
                  'Hernandez-Segura et al Up',
                  'Hernandez-Segura et al Down',
                  'Casella et al Up',
                  'Casella et al Down',
                  'Cherry et al Up',
                  'Cherry et al Down'
                  # 'Tao et al SID1',
                  # 'Tao et al SID2',
                  # 'Tao et al SID3',
                  # 'Tao et al SID4',
                  # 'Tao et al SID5',
                  # 'Tao et al SID6'
                  )
# save_csv(cs_signatures_studies,file_name = '1_cellage',path = save_dir_csv)
arrest_degs_merged$arrest=ifelse(grepl(arrest_degs_merged$group_1,pattern = 'CS'),'CS','CQ')

arrest_degs_merged_sig=arrest_degs_merged[arrest_degs_merged$sig=='y',]
#####
#Upset plots

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
                  facet_col='accession',
                  gene_col='gene',
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
# factor(arrest_degs_merged$accession,
#        levels=c('Contact_inhibited_CQ vs\nProliferating',
#                 'Serum_starved_CQ vs\nProliferating',
#                 'Replicative_CS vs\nProliferating',
#                 'Stress_induced_CS vs\nProliferating',
#                 'Oncogene_induced_CS vs\nProliferating'))

arrest_degs_merged_sig$group_1=
  as.character(arrest_degs_merged_sig$group_1)

# cs_signatures_studies_just_signatures=cs_signatures_studies[!cs_signatures_studies$source%in%
#                                                               c('Driver','Signature'),]
# cs_signatures_studies_just_signatures$source[cs_signatures_studies_just_signatures$source==
#                                                'Saul et al']='SenMayo'
#overlap no direction
# deg_cellage_overlap_no_dir=overlap_function(df_1=arrest_degs_merged_sig,
#                                      df_2=cs_signatures_studies,
#                                      gene_col_1='gene',
#                                      gene_col_2='gene',
#                                      group_col_1=c('group_1'),
#                                      group_col_2=c('source'),
#                                      background=arrest_degs_merged$gene)
# 
# deg_cellage_overlap_no_dir$pseudo=1
# recount_cellage_overlap_plot_no_dir=create_overlap_plot(deg_db_overlap =
#                                                           deg_cellage_overlap_no_dir,
#                                                  odds_column = 'odds',
#                                                  pval_col = 'adj',
#                                                  facet_col = 'group_1',
#                                                  x = 'pseudo',
#                                                  y = 'source',
#                                                  xlab = 'Cell Cycle Arrest DEGs',
#                                                  ylab = 'Gene List',text_size = 15,)
# 
# save_p(recount_cellage_overlap_plot_no_dir,
#        file_name = 'signature_degs_overlap_no_dir',
#        save_dir = '/Users/ravelarvargas/Downloads',p_width = 11)

#overlap direction
arrest_degs_merged_sig$group_1[arrest_degs_merged_sig$group_1=="Replicative_CS"]='RS'
arrest_degs_merged_sig$group_1[arrest_degs_merged_sig$group_1=="Oncogene_induced_CS"]='OIS'
arrest_degs_merged_sig$group_1[arrest_degs_merged_sig$group_1=="Stress_induced_CS"]='SIPS'
arrest_degs_merged_sig$group_1[arrest_degs_merged_sig$group_1=="Contact_inhibited_CQ"]='CICQ'
arrest_degs_merged_sig$group_1[arrest_degs_merged_sig$group_1=="Serum_starved_CQ"]='SSCQ'
arrest_degs_merged_sig$cs_accession=paste0(arrest_degs_merged_sig$group_1,'_',
                                        arrest_degs_merged_sig$direction_1)
arrest_degs_merged_sig$group_1=factor(arrest_degs_merged_sig$group_1,
                                      levels=c('CICQ','OIS','RS','SSCQ','SIPS'))
cs_gene_list$accession=factor(cs_gene_list$accession,levels=
                                    rev(cs_factor_level))
deg_cellage_overlap=overlap_function(df_1=arrest_degs_merged_sig,
                                     df_2=cs_gene_list,
                                     gene_col_1='gene',
                                     gene_col_2='gene',
                                     group_col_1 = 'cs_accession',
                                     carry_col_1 =c('group_1','direction_1'),
                                     carry_col_2=c('dir','source'),
                                     group_col_2='accession',
                                     background=arrest_degs_merged$gene)
deg_cellage_overlap_upset=deg_cellage_overlap

recount_cellage_overlap_plot=create_overlap_plot(deg_db_overlap = deg_cellage_overlap,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'group_1',
                    x = 'direction_1',
                    y = 'accession',
                    xlab = 'Cell Cycle Arrest DEGs',
                    ylab = 'CS Gene Lists',text_size = 15,remove_nonsig = TRUE)
save_csv(data = deg_cellage_overlap,
         file_name = 'recount_databases_overlap',
         path = save_dir_csv)
save_p(plot = recount_cellage_overlap_plot,
       file_name = 'recount_cs_db_overlap',
       save_dir = save_dir_figure,p_width = 12,p_height=10.5)

deg_cellage_overlap$cs_accession=gsub(deg_cellage_overlap$cs_accession,pattern='_',replacement=' ')
deg_cellage_overlap$cs_accession=gsub(deg_cellage_overlap$cs_accession,
                                      pattern='up',replacement='U')
deg_cellage_overlap$cs_accession=gsub(deg_cellage_overlap$cs_accession,
                                      pattern='down',replacement='D')
deg_cellage_overlap$cs_accession=factor(deg_cellage_overlap$cs_accession,
                                        levels=rev(c("RS U","SIPS U","OIS U","CICQ U","SSCQ U",
                                                     "RS D","SIPS D","OIS D",
                                                     "CICQ D","SSCQ D")))

deg_cellage_overlap$accession=gsub(deg_cellage_overlap$accession,pattern='Induces',replacement='Induce')
deg_cellage_overlap$accession=gsub(deg_cellage_overlap$accession,pattern='Inhibits',replacement='Inhibit')
deg_cellage_overlap$accession=factor(deg_cellage_overlap$accession,
                  levels=c('CellAge Induce CS',
                           'CellAge Inhibit CS',
                           'CellAge Up in RS',
                           'CellAge Down in RS',
                           'SenMayo',
                           "Hernandez-Segura et al Up",
                           "Hernandez-Segura et al Down",
                           "Casella et al Up",
                           "Casella et al Down",
                           'Cherry et al Up',
                           'Cherry et al Down'))
simple_deg_v_db=summarise_overlaps(db = deg_cellage_overlap,
               x = 'accession',
               y = 'cs_accession',
               xlab='CS Genelists',
               ylab = 'Arrest-DEGs')

save_p(simple_deg_v_db,file_name = '2b_deg_v_db',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',p_width = 6,p_height = 3.5)
#####
#self overlaps
self_overlaps_recount_to_plot=overlap_function(df_1 = arrest_degs_merged_sig,
                                               df_2 = arrest_degs_merged_sig,
                                               gene_col_1 = 'gene',
                                               gene_col_2 = 'gene',
                                               carry_col_1 = c('group_1','group_1_dir'),
                                               group_col_1 = 'dir_accession',
                                               group_col_2 = 'dir_accession',
                                               carry_col_2 = c('group_1','group_1_dir'),background = arrest_degs_merged$gene
)

self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Oncogene_induced_CS_up']=
  'OIS Up'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Oncogene_induced_CS_down']=
  'OIS Down'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Stress_induced_CS_up']=
  'SIPS Up'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Stress_induced_CS_down']=
  'SIPS Down'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Replicative_CS_up']=
  'RS Up'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Replicative_CS_down']=
  'RS Down'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Contact_inhibited_CQ_up']=
  'CICQ Up'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Contact_inhibited_CQ_down']=
  'CICQ Down'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Serum_starved_CQ_up']=
  'SSCQ Up'
self_overlaps_recount_to_plot$group_1_dir[self_overlaps_recount_to_plot$group_1_dir=='Serum_starved_CQ_down']=
  'SSCQ Down'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Oncogene_induced_CS_up']=
  'OIS Up'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Oncogene_induced_CS_down']=
  'OIS Down'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Stress_induced_CS_up']=
  'SIPS Up'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Stress_induced_CS_down']=
  'SIPS Down'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Replicative_CS_up']=
  'RS Up'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Replicative_CS_down']=
  'RS Down'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Contact_inhibited_CQ_up']=
  'CICQ Up'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Contact_inhibited_CQ_down']=
  'CICQ Down'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Serum_starved_CQ_up']=
  'SSCQ Up'
self_overlaps_recount_to_plot$group_1_dir.1[self_overlaps_recount_to_plot$group_1_dir.1=='Serum_starved_CQ_down']=
  'SSCQ Down'

adj_temp=do.call('rbind',lapply(self_overlaps_recount_to_plot$pval,function(per_pval){
  p.adjust(p = per_pval,method = 'BH',n = 40)
}))
self_overlaps_recount_to_plot$adj=adj_temp
self_overlaps_recount_to_plot$log2odds=log2(self_overlaps_recount_to_plot$odds)
self_overlaps_recount_to_plot=self_overlaps_recount_to_plot%>%rstatix::add_significance('adj')
self_overlaps_recount_to_plot$adj.signif[is.infinite(self_overlaps_recount_to_plot$log2odds)]=''
self_overlaps_recount_to_plot$actual[self_overlaps_recount_to_plot$actual==0]=''
# Determine graph_max and graph_min
graph_max = max(ceiling(self_overlaps_recount_to_plot$log2odds/2)[!is.infinite(ceiling(self_overlaps_recount_to_plot$log2odds))]) * 2
graph_min = min(floor(self_overlaps_recount_to_plot$log2odds/2)[!is.infinite(ceiling(self_overlaps_recount_to_plot$log2odds))])*2
if(graph_min > 0){
  graph_min = 0
}

li <- c(graph_min, graph_max)
la <- c(seq(graph_min,graph_max,2))
br <- c(seq(graph_min,graph_max,2))

self_overlaps_recount_to_plot$group_1_dir=factor(self_overlaps_recount_to_plot$group_1_dir,
                                                 levels=c('CICQ Up',
                                                          'CICQ Down',
                                                          'SSCQ Up',
                                                          'SSCQ Down',
                                                          'RS Up',
                                                          'RS Down',
                                                          'SIPS Up',
                                                          'SIPS Down',
                                                          'OIS Up',
                                                          'OIS Down'))
self_overlaps_recount_to_plot$group_1_dir.1=factor(self_overlaps_recount_to_plot$group_1_dir.1,
                                                 levels=rev(c('CICQ Up',
                                                          'CICQ Down',
                                                          'SSCQ Up',
                                                          'SSCQ Down',
                                                          'RS Up',
                                                          'RS Down',
                                                          'SIPS Up',
                                                          'SIPS Down',
                                                          'OIS Up',
                                                          'OIS Down')))

self_overlaps_recount_to_plot$adj.signif[self_overlaps_recount_to_plot$adj.signif=='ns']=''
self_overlap_p=self_overlaps_recount_to_plot%>%
  ggplot(aes(x=group_1_dir, y=group_1_dir.1)) +
  geom_tile(aes(fill=log2odds),colour='black') +
  scale_fill_gradient2(expression('log'[2]*'(Odds)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks = TRUE, 
                                              ticks.colour = 'black', title.position = 'top', title.hjust=0.5),
                       breaks=br,
                       labels=la,
                       limits=li) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_text(aes(label=adj.signif), nudge_y=0.25, size = 15 / 3) +
  geom_text(aes(label=actual), nudge_y=-0.25, size = 15 / 3) +
  xlab('Condition 1')+
  ylab('Condition 2')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y=element_text(size=12))

save_p(self_overlap_p,file_name = 'self_overlap_recount',
       save_dir = save_dir_figure,p_width=7,p_height = 7)

self_overlaps_recount_to_plot=clean_self_overlaps(df = self_overlaps_recount_to_plot,
                                                  group_1 = 'group_1_dir',
                                                  group_2 = 'group_1_dir.1')
save_csv(self_overlaps_recount_to_plot,
         file_name = 'self_overlap_recount',
         path = save_dir_csv)
#table DEGs
DEG_count=data.frame(table(arrest_degs_merged_sig%>%dplyr::select(group_1,direction_1)))
save_csv(DEG_count,file_name = 'deg_count',path = save_dir_csv)

#simplify
self_overlaps_recount_to_plot$group_1_dir=
  gsub(self_overlaps_recount_to_plot$group_1_dir,pattern='Up',replacement='U')
self_overlaps_recount_to_plot$group_1_dir=
  gsub(self_overlaps_recount_to_plot$group_1_dir,pattern='Down',replacement='D')
self_overlaps_recount_to_plot$group_1_dir.1=
  gsub(self_overlaps_recount_to_plot$group_1_dir.1,pattern='Up',replacement='U')
self_overlaps_recount_to_plot$group_1_dir.1=
  gsub(self_overlaps_recount_to_plot$group_1_dir.1,pattern='Down',replacement='D')

self_overlaps_recount_to_plot$group_1_dir=factor(self_overlaps_recount_to_plot$group_1_dir,
                                        levels=c("CICQ U","SSCQ U","RS U",
                                                 "SIPS U","OIS U",
                                                 "CICQ D","SSCQ D",
                                                 "RS D","SIPS D",
                                                 "OIS D"))
self_overlaps_recount_to_plot$group_1_dir.1=factor(self_overlaps_recount_to_plot$group_1_dir.1,
                                                 levels=rev(c("CICQ U","SSCQ U","RS U",
                                                              "SIPS U","OIS U",
                                                              "CICQ D","SSCQ D",
                                                              "RS D","SIPS D",
                                                              "OIS D")))
self_overlap_simple=summarise_overlaps(db = self_overlaps_recount_to_plot,
                                        x = 'group_1_dir',
                                        y = 'group_1_dir.1',
                                        xlab='Arrest-DEGs 1',
                                        ylab = 'Arrest-DEGs 2',
                                       add_line = TRUE,label_cq = FALSE,
                                       remove_self = TRUE,
                                       self_col = c('group_1','group_1.1'),
                                       legend_position_bottom = TRUE,rotate_x_text = FALSE)
self_overlap_simple=self_overlap_simple+
  geom_vline(xintercept = 5.5,size=1)

save_p(self_overlap_simple,file_name = '2a_self_simple',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',p_width = 6.5,p_height = 2.7)
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

arrest_degs_merged_contact=arrest_degs_merged_cq[arrest_degs_merged_cq$group_1=='Contact-inhibited CQ',]
volcano_plot_contact=volcano_function(degs=arrest_degs_merged_contact,
                                 gene_col='gene',title = 'Contact-inhibited CQ')

arrest_degs_merged_serum=arrest_degs_merged_cq[arrest_degs_merged_cq$group_1=='Serum-starved CQ',]
volcano_plot_serum=volcano_function(degs=arrest_degs_merged_serum,
                                      gene_col='gene',title = 'Serum-starved CQ')
save_p(plot=volcano_plot_contact/volcano_plot_serum,
       file_name = 'volcano_all_cq',save_dir = save_dir_figure,p_height=10)
#####
#intersect enrichment
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
# compare_logfc_outgroup(degs=arrest_degs_merged,
#                        ingroup_col = 'group_1',
#                        outgroup_col = 'group_2',
#                        degs_gene_col = 'gene',
#                        outgroup='Proliferating',
#                        # xlab='Day',
#                        # facet_order=levels(time_degs$group_1),
#                        database=cellage,
#                        database_gene_col = 'gene_name',
#                        database_group_col='Senescence')
# compare_logfc_outgroup(degs=arrest_degs_merged,
#                        ingroup_col = 'group_1',
#                        outgroup_col = 'group_2',
#                        degs_gene_col = 'gene',
#                        outgroup='Proliferating',
#                        # xlab='Day',
#                        # facet_order=levels(time_degs$group_1),
#                        database=common_background_all,
#                        database_gene_col = 'gene',
#                        database_group_col='dir')

#####
# upset test
casella_temp=deg_cellage_overlap_upset[deg_cellage_overlap_upset$source=='Casella et al',]

by_col='dir'

casella_temp$accession=paste0(casella_temp$group_1,'_',casella_temp$direction_1)

overlap=casella_temp
upset_intersects=function(overlap,
         overlap_col='int',
         overlap_sep='/',
         by_col=NULL,
         group_col='accession',
         actual_col='actual',
         title=NULL){
  results=c()
  if(is.null(by_col)){
    overlap$by='test'
    by_col='test'
  }
  all_by=unique(overlap[[by_col]])
  
  lapply(all_by,function(per_by){
    overlap_temp=overlap[overlap[[by_col]]==per_by,]
    overlap_temp=overlap_temp[overlap_temp[[actual_col]]>0,]
    genes_to_upset=do.call('rbind',apply(overlap_temp,1,function(get_genes){
      genes_use=data.frame(gene=unlist(strsplit(get_genes[[overlap_col]],split=overlap_sep)))
      genes_use$group=get_genes[[group_col]]
      return(genes_use)
    }))
    df_wide <- genes_to_upset %>%
      # This will create a new column with all values as 1 (indicating presence)
      mutate(value = 1) %>%
      # Spread the data to wide format, filling absent cases with 0
      spread(key = group, value = value, fill = 0)
    name_use=ifelse(is.null(title),per_by,paste0(title,' ',per_by))
    df_wide[['condition']]=name_use
    results[[paste0(per_by,'_pivot')]]<<-df_wide
    df_wide=df_wide[,!colnames(df_wide)%in%c('gene','condition')]
    
    upset_plot <- upset(
      df_wide,intersect = colnames(df_wide),
      name = name_use
    )
    results[[per_by]]<<-upset_plot
  })
  return(results)
}

upset_df=casella_upset$down_pivot
col_find=c('Oncogene_induced_CS_down',
           'Replicative_CS_down',
           'Stress_induced_CS_down')

col_exclude=c('Contact_inhibited_CQ_down',
              'Serum_starved_CQ_down')

#Function to find common genes, with the ability to exclude genes present in other groups
filter_upset=function(upset_df,
                      gene_col='gene',
                      col_find,
                      col_exclude=NULL){
  if(!is.null(col_exclude)){
    exclude_upset=which(rowSums(upset_df[,c(col_exclude)])>0)
    upset_df=upset_df[-exclude_upset,]
  }
  upset_use=upset_df[,colnames(upset_df)%in%c(gene_col,col_find)]
  upset_use[which(rowSums(upset_use[,col_find])==length(col_find)),][[gene_col]]
}

casella_upset=upset_intersects(overlap = casella_temp,
                 by_col='dir',title = 'Casella')
save_p(casella_upset$down,file_name = 'casella_upset_down',
       save_dir = '/Users/ravelarvargas/Downloads')
save_p(casella_upset$up,file_name = 'casella_upset_up',
       save_dir = '/Users/ravelarvargas/Downloads')


#just up direction
deg_cellage_overlap_upset_up=deg_cellage_overlap_upset[deg_cellage_overlap_upset$direction_1=='up',]
senmayo_temp=deg_cellage_overlap_upset_up[deg_cellage_overlap_upset_up$source=='Saul et al',]
senmayo_temp$accession=paste0(senmayo_temp$group_1,'_',senmayo_temp$direction_1)

senmayo_upset=upset_intersects(overlap = senmayo_temp,
                 title = 'SenMayo',by_col = 'dir')
save_p(senmayo_upset[[2]],file_name = 'senmayo_upset',
       save_dir = save_dir_figure_si,p_width=10)

save_csv(senmayo_upset[[1]],file_name = 'SenMayo_upset',path = save_dir_csv)

upset_filtered=filter_upset(senmayo_upset$up_pivot,
             col_find = c('Oncogene_induced_CS_up',
                                   'Replicative_CS_up',
                                   'Stress_induced_CS_up'),
             col_exclude=c('Contact_inhibited_CQ_up',
                           'Serum_starved_CQ_up'))
