#CQ and CS recount3 studies
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')

exclude_accession=unique(c('SRP113334',
                           'SRP195418',
                           'SRP127595',
                           'SRP050179'))

study_coldata_cq_filtered=read.csv('/Users/ravelarvargas/Downloads/reversible_arrest/cq_study_data_final.csv')
study_coldata_cq_filtered=study_coldata_cq_filtered[!study_coldata_cq_filtered$sra%in%
                                                      exclude_accession,]
recount_cs=read.csv('/Users/ravelarvargas/Downloads/reversible_arrest/cs/recount_final.csv')
save_dir='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/'
save_dir_csv=paste0(save_dir,'SI_tables/')
save_dir_figure=paste0(save_dir,'Figures/')
save_dir_figure_si=paste0(save_dir,'SI_figures/')
#####
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)
#####

cq_studies=c('SRP017378' ,'SRP065206' ,'SRP045867' ,'SRP052706' ,'SRP096629',
             'SRP154577' ,'SRP066947' ,'ERP021140' ,'SRP089801' ,'SRP153205',
             'SRP153724' ,'SRP154382')
#Filter for excluded studies
cq_studies=cq_studies[!cq_studies%in%exclude_accession]

# #Download studies
# cq_studies_full=download_studies(cq_studies)
# # 
# #Process studies
# cq_studies_processed=process_rse(rse=cq_studies_full,
#                                  meta_add=study_coldata_cq_filtered,
#                                  filter_pc=TRUE,
#                                  ensembl_dictionary=human_pc,
#                                  filter_duplicates=TRUE,
#                                  sample_filter=study_coldata_cq_filtered[['sample_ID']],
#                                  low_expression_filter='per',
#                                  scale_counts=TRUE)
# # 
# # # Save study data
# study_coldata_cq=data.frame(colData(cq_studies_processed))
# save_csv(data = study_coldata_cq,file_name = 'study_info_cq',
#          path=save_dir_csv)
# 
# # saveRDS(cq_studies_processed,paste0(save_dir_csv,'cq_samples.rds'))
# cq_studies_processed=readRDS(paste0(save_dir_csv,'cq_samples.rds'))
# 
# #build model matrix
# cq_studies_processed_deseq=deseq_studies(study_obj=cq_studies_processed,
#                                          protect_col='cell_substate',
#                                          fix_col='study')

#####
# #pca
# #No batch correction
# cq_no_correction_pca_counts=normalise_counts(
#   dds_obj_use=cq_studies_processed_deseq[['deseq_obj']],
#   batch_col=NULL,
#   protect_col=NULL,
#   method=NULL,
#   normalisation='vst',
#   vst_n=500,
#   blind=TRUE)
# ##pca
# cq_pca_no_correction=pca_rds(
#   cq_no_correction_pca_counts,
#   pca_plot=c('cell_line',
#              'cell_type','tissue','cell_state','cell_substate',
#              'sra.platform_model'),
#   ntop=500)
# 
# save_p(plot = cq_pca_no_correction,
#        file_name = 'cq_pca_no_correction',
#        save_dir = save_dir_figure_si,p_width = 10,p_height=10)
# 
# cq_pca_no_correction_study=pca_rds(
#   cq_no_correction_pca_counts,
#   pca_plot=c('study'),
#   ntop=500)
# 
# save_p(plot = cq_pca_no_correction_study,
#        file_name = 'cq_pca_no_correction_study',
#        save_dir = save_dir_figure_si,p_width = 6,p_height=8)
# 
# cq_correction_pca_counts=normalise_counts(
#   dds_obj_use=cq_studies_processed_deseq[['deseq_obj']],
#   batch_col='study',
#   protect_col=c('cell_substate'),
#   method='wgcna',
#   normalisation='vst',
#   vst_n=500,
#   blind=TRUE)
# 
# cq_correction_pca=pca_rds(
#   cq_correction_pca_counts,
#   pca_plot=c('cell_line',
#              'cell_type','tissue','cell_state','cell_substate',
#              'sra.platform_model'),
#   ntop=500)
# 
# save_p(plot = cq_correction_pca,
#        file_name = 'cq_pca_correction',
#        save_dir = save_dir_figure_si,p_width = 10,p_height=10)
# 
# cq_pca_correction_study=pca_rds(
#   cq_correction_pca_counts,
#   pca_plot=c('study'),
#   ntop=500)
# 
# save_p(plot = cq_pca_correction_study,
#        file_name = 'cq_pca_correction_study',
#        save_dir = save_dir_figure_si,p_width = 6,p_height=8)
# 
# cq_correction_pca_manuscript=pca_rds(
#   cq_correction_pca_counts,
#   pca_plot = c('cell_state','cell_substate'),
#   ntop=500)
# 
# save_p(plot = cq_correction_pca_manuscript,
#        file_name = 'cq_pca',
#        save_dir = save_dir_figure,p_width = 10,p_height=10)

#####
#DEGs

# #CQ
# cq_studies_processed$cell_substate[cq_studies_processed$cell_substate=='Serum_starved CQ']='Serum-starved CQ'
# cq_studies_processed$cell_substate[cq_studies_processed$cell_substate=='Contact_inhibited CQ']='Contact-inhibited CQ'
# ##contact
# contact_cq_degs=find_degs_new(rse_obj=cq_studies_processed,
#                               group_col='cell_substate',
#                               group_1='Contact-inhibited CQ',
#                               group_2='Proliferating',
#                               batch_col='study')
# 
# cq_contact_heatmap=build_heatmap(degs = contact_cq_degs$degs[contact_cq_degs$degs$sig=='y',],
#                                  dds_obj = contact_cq_degs$dds_deseq$deseq_obj,
#                                  group_col = 'cell_substate',
#                                  groups = c('Contact-inhibited CQ','Proliferating'),
#                                  protect_col = 'cell_substate',
#                                  batch_col = 'study',
#                                  normalisation = 'vst_full',
#                                  group_col_heatmap=c('cell_substate',
#                                                      'study','tissue'),
#                                  gene_col = 'gene',
#                                  rownames_to_symbol = FALSE)
# 
# cq_contact_heatmap_pi=pi_score_degs(sig_degs = contact_cq_degs$degs[contact_cq_degs$degs$sig=='y',])
# save_csv(data = cq_contact_heatmap_pi$all_pi,
#         file_name ='cq_contact_pi',path = save_dir_csv)
# 
# save_pheatmap(pheatmap_to_save=cq_contact_heatmap,
#               file_name='contact_cq_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=3000,
#               p_height=6000)
# 
# ##serum
# serum_cq_degs=find_degs_new(rse_obj=cq_studies_processed,
#                             group_col='cell_substate',
#                             group_1='Serum-starved CQ',
#                             group_2='Proliferating',
#                             batch_col='study')
# 
# cq_serum_heatmap=build_heatmap(degs = serum_cq_degs$degs[serum_cq_degs$degs$sig=='y',],
#                                dds_obj = serum_cq_degs$dds_deseq$deseq_obj,
#                                group_col = 'cell_substate',
#                                groups = c('Serum-starved CQ','Proliferating'),
#                                protect_col = 'cell_substate',
#                                batch_col = 'study',
#                                normalisation = 'vst_full',
#                                gene_col = 'gene',
#                                rownames_to_symbol = FALSE)
# cq_serum_heatmap_pi=pi_score_degs(sig_degs = serum_cq_degs$degs[serum_cq_degs$degs$sig=='y',])
# save_csv(data = cq_serum_heatmap_pi$all_pi,
#          file_name ='cq_serum_pi',path = save_dir_csv)
# 
# save_pheatmap(pheatmap_to_save=cq_serum_heatmap,
#               file_name='serum_cq_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=3000,
#               p_height=6000)
# 
# #all cq
# cq_degs=rbind(contact_cq_degs$degs,serum_cq_degs$degs)
# 
# dds_obj_cq=find_degs_new(rse_obj=cq_studies_processed,
#                          group_col='cell_substate',
#                          group_1=c('Contact-inhibited CQ',
#                                    'Serum-starved CQ'),
#                          group_2='Proliferating',
#                          batch_col='study',
#                          find_degs = FALSE)
# 
# cq_heatmap=build_heatmap(degs = cq_degs[cq_degs$sig=='y',],
#                          dds_obj = dds_obj_cq$deseq_obj,
#                          group_col = 'cell_substate',
#                          groups = c('Serum-starved CQ',
#                                     'Contact-inhibited CQ',
#                                     'Proliferating'),
#                          protect_col = 'cell_substate',
#                          batch_col = 'study',
#                          normalisation = 'vst_full',
#                          group_col_heatmap = c('tissue','cell_substate',
#                                                'cell_state'),
#                          top_x = 15,
#                          gene_col = 'gene',
#                          rownames_to_symbol = FALSE)
# 
# save_pheatmap(pheatmap_to_save=cq_heatmap,
#               file_name='cq_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=2000,
#               p_height=2500)
# 
# save_csv(cq_degs,file_name = 'cq_degs',path = save_dir_csv)

#####
recount_cs=recount_cs[!recount_cs[['sra']]%in%exclude_accession,]
cs_studies=c('SRP017142', 'SRP034163', 'SRP034541', 'SRP040243', 'SRP040745',
             'SRP045867', 'SRP046254', 'SRP050179', 'SRP060598', 'SRP062872',
             'SRP064207', 'SRP065206', 'SRP066917', 'SRP066947', 'SRP069768',
             'SRP070636', 'SRP098713', 'SRP107235', 'SRP113324', 'SRP113329',
             'SRP113334', 'SRP117883', 'SRP121031', 'ERP021140', 'SRP052706',
             'SRP089801', 'SRP096629', 'SRP153205', 'SRP153724', 'SRP154382',
             'SRP154577', 'SRP123346', 'SRP127037', 'SRP127595', 'SRP136071',
             'SRP136727', 'SRP172671', 'SRP195418')
cs_studies=cs_studies[!cs_studies%in%exclude_accession]

# cs_studies_full=download_studies(studies=cs_studies,
#                                          sra_organism='human')
# 
# cs_studies_full_processed=process_rse(rse=cs_studies_full,
#                                               meta_add=recount_cs,
#                                               filter_pc=TRUE,
#                                               ensembl_dictionary=human_pc,
#                                               filter_duplicates=TRUE,
#                                               sample_filter=recount_cs[['sample_ID']],
#                                               low_expression_filter='per',
#                                               scale_counts=TRUE)
# 
# study_coldata_cs=data.frame(colData(cs_studies_full_processed))
# save_csv(data = study_coldata_cs,
#          file_name = 'study_info_cs',
#          path = save_dir_csv)
# 
# cs_studies_full_processed_deseq=deseq_studies(study_obj=cs_studies_full_processed,
#                                                       protect_col='cell_substate',
#                                                       fix_col='study')
# 
# cs_no_correction_pca_counts=normalise_counts(
#   dds_obj_use=cs_studies_full_processed_deseq[['deseq_obj']],
#   batch_col=NULL,
#   protect_col=NULL,
#   method=NULL,
#   normalisation='vst',
#   vst_n=500,
#   blind=TRUE)
# 
# #####
# #pca
# #No batch correction
# cs_no_correction_pca_counts=normalise_counts(
#   dds_obj_use=cs_studies_full_processed_deseq[['deseq_obj']],
#   batch_col=NULL,
#   protect_col=NULL,
#   method=NULL,
#   normalisation='vst',
#   vst_n=500,
#   blind=TRUE)
# ##pca
# cs_pca_no_correction=pca_rds(
#   cs_no_correction_pca_counts,
#   pca_plot=c('cell_line',
#              'cell_type','tissue','cell_state','cell_substate',
#              'sra.platform_model'),
#   ntop=500)
# 
# save_p(plot = cs_pca_no_correction,
#        file_name = 'cs_pca_no_correction',
#        save_dir = save_dir_figure_si,p_width = 10,p_height=10)
# 
# cs_pca_no_correction_study=pca_rds(
#   cs_no_correction_pca_counts,
#   pca_plot=c('study'),
#   ntop=500)
# 
# save_p(plot = cs_pca_no_correction_study,
#        file_name = 'cs_pca_no_correction_study',
#        save_dir = save_dir_figure_si,p_width = 6,p_height=8)
# 
# cs_correction_pca_counts=normalise_counts(
#   dds_obj_use=cs_studies_full_processed_deseq[['deseq_obj']],
#   batch_col='study',
#   protect_col=c('cell_substate'),
#   method='wgcna',
#   normalisation='vst',
#   vst_n=500,
#   blind=TRUE)
# 
# cs_correction_pca=pca_rds(
#   cs_correction_pca_counts,
#   pca_plot=c('cell_line',
#              'cell_type','tissue','cell_state','cell_substate',
#              'sra.platform_model'),
#   ntop=500)
# 
# save_p(plot = cs_correction_pca,
#        file_name = 'cs_pca_correction',
#        save_dir = save_dir_figure_si,p_width = 13,p_height=10)
# 
# cs_pca_correction_study=pca_rds(
#   cs_correction_pca_counts,
#   pca_plot=c('study'),
#   ntop=500)
# 
# save_p(plot = cs_pca_correction_study,
#        file_name = 'cs_pca_correction_study',
#        save_dir = save_dir_figure_si,p_width = 6,p_height=8)
# 
# cs_correction_pca_manuscript=pca_rds(
#   cs_correction_pca_counts,
#   pca_plot = c('cell_state','cell_substate'),
#   ntop=500)
# 
# save_p(plot = cs_correction_pca_manuscript,
#        file_name = 'cs_pca',
#        save_dir = save_dir_figure,p_width = 10,p_height=10)
# 
# #####
# #DEGs
# rs_degs=find_degs_new(rse_obj=cs_studies_full_processed,
#                       group_col='cell_substate',
#                       group_1='Replicative CS',
#                       group_2='Proliferating',
#                       batch_col='study')
# 
# rs_heatmap=build_heatmap(degs = rs_degs$degs[rs_degs$degs$sig=='y',],
#                          dds_obj = rs_degs$dds_deseq$deseq_obj,
#                          group_col = 'cell_substate',
#                          groups = c('Replicative CS','Proliferating'),
#                          protect_col = 'cell_substate',
#                          batch_col = 'study',
#                          normalisation = 'vst_full',
#                          gene_col = 'gene',
#                          rownames_to_symbol = FALSE)
# 
# rs_heatmap_pi=pi_score_degs(sig_degs = rs_degs$degs[rs_degs$degs$sig=='y',])
# save_csv(data =rs_heatmap_pi$all_pi,
#         file_name='rs_contact_pi',path=save_dir_csv)
# 
# save_pheatmap(pheatmap_to_save=rs_heatmap,
#               file_name='rs_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=4000,
#               p_height=5000)
# 
# ois_degs=find_degs_new(rse_obj=cs_studies_full_processed,
#                        group_col='cell_substate',
#                        group_1='Oncogene-induced CS',
#                        group_2='Proliferating',
#                        batch_col='study')
# 
# ois_heatmap=build_heatmap(degs = ois_degs$degs[ois_degs$degs$sig=='y',],
#                           dds_obj = ois_degs$dds_deseq$deseq_obj,
#                           group_col = 'cell_substate',
#                           groups = c('Oncogene-induced CS','Proliferating'),
#                           protect_col = 'cell_substate',
#                           batch_col = 'study',
#                           normalisation = 'vst_full',
#                           gene_col = 'gene',
#                           rownames_to_symbol = FALSE)
# ois_heatmap_pi=pi_score_degs(sig_degs = ois_degs$degs[ois_degs$degs$sig=='y',])
# save_csv(data =ois_heatmap_pi$all_pi,
#         file_name='ois_pi',path=save_dir_csv)
# 
# save_pheatmap(pheatmap_to_save=ois_heatmap,
#               file_name='ois_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=4000,
#               p_height=5000)
# 
# sips_degs=find_degs_new(rse_obj=cs_studies_full_processed,
#                         group_col='cell_substate',
#                         group_1='Stress-induced CS',
#                         group_2='Proliferating',
#                         batch_col='study')
# 
# sips_heatmap=build_heatmap(degs = sips_degs$degs[sips_degs$degs$sig=='y',],
#                            dds_obj = sips_degs$dds_deseq$deseq_obj,
#                            group_col = 'cell_substate',
#                            groups = c('Stress-induced CS','Proliferating'),
#                            protect_col = 'cell_substate',
#                            batch_col = 'study',
#                            normalisation = 'vst_full',
#                            gene_col = 'gene',
#                            rownames_to_symbol = FALSE)
# sips_heatmap_pi=pi_score_degs(sig_degs = sips_degs$degs[sips_degs$degs$sig=='y',])
# save_csv(data =sips_heatmap_pi$all_pi,
#         file_name='sips_pi',path=save_dir_csv)
# 
# save_pheatmap(pheatmap_to_save=sips_heatmap,
#               file_name='sips_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=4000,
#               p_height=5000)
# 
# cs_degs=rbind(rs_degs$degs,ois_degs$degs,sips_degs$degs)
# 
# dds_obj_cs=find_degs_new(rse_obj=cs_studies_full_processed,
#                          group_col='cell_substate',
#                          group_1=c('Stress-induced CS',
#                                    'Oncogene-induced CS',
#                                    'Replicative CS'),
#                          group_2='Proliferating',
#                          batch_col='study',
#                          find_degs = FALSE)
# 
# cs_heatmap=build_heatmap(degs = cs_degs[cs_degs$sig=='y',],
#                          dds_obj = dds_obj_cs$deseq_obj,
#                          group_col = 'cell_substate',
#                          groups = c('Stress-induced CS',
#                                     'Oncogene-induced CS',
#                                     'Replicative CS',
#                                     'Proliferating'),
#                          protect_col = 'cell_substate',
#                          batch_col = 'study',
#                          normalisation = 'vst_full',
#                          group_col_heatmap = c('tissue','cell_substate',
#                                                'cell_state'),
#                          top_x = 15,
#                          gene_col = 'gene',
#                          rownames_to_symbol = FALSE)
# 
# save_pheatmap(pheatmap_to_save=cs_heatmap,
#               file_name='cs_vs_proliferative_heatmap',
#               save_dir=save_dir_figure_si,
#               p_width=4000,
#               p_height=4000)
# 
# save_csv(cs_degs,file_name = 'cs_degs',path = save_dir_csv)

#####
#both
# meta_both=unique(rbind(study_coldata_cq_filtered,
#                 recount_cs))
meta_both=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/study_info_all.csv')
meta_both$cq_test=meta_both$cell_substate
temp_label=unique(meta_both$cq_test)[grepl(unique(meta_both$cq_test),pattern='CQ')]
meta_both$
  cq_test=ifelse(meta_both$cell_substate%in%temp_label,'CQ',meta_both$cell_substate)
# 
# cs_cq=list(c(meta_both$sra),meta_both$sample_ID)
# names(cs_cq)=c('study','sample')

cs_cq_download=download_studies(studies=unique(meta_both$study),
                                sra_organism='human')

cs_cq_all_study_processed=process_rse(rse=cs_cq_download,
                                      meta_add=meta_both,
                                      meta_id_col='external_id',
                                      filter_pc=TRUE,
                                      ensembl_dictionary=human_pc,
                                      filter_duplicates=TRUE,
                                      sample_filter=unique(meta_both[['external_id']]),
                                      low_expression_filter='per',
                                      scale_counts=TRUE)

# saveRDS(cs_cq_all_study_processed,file = '/Users/ravelarvargas/Downloads/marian/cs_cq_all_study_processed_degs.rds')
cs_cq_all_study_processed=readRDS('/Users/ravelarvargas/Downloads/marian/cs_cq_all_study_processed_degs.rds')

cell_counts_table=data.frame(table(data.frame(colData(cs_cq_all_study_processed))%>%dplyr::select(tissue, cell_substate)))
save_csv(cell_counts_table,file_name ='sample_subtype_counts.csv',
         path = save_dir_csv)

temp_meta=data.frame(colData(cs_cq_all_study_processed))
save_csv(temp_meta,file_name ='sample_metadata.csv',
         path = save_dir_csv)

cs_cq_studies_processed_deseq=deseq_studies(study_obj=cs_cq_all_study_processed,
                                            protect_col='cell_substate',
                                            fix_col='study')

# saveRDS(cs_cq_studies_processed_deseq,
#         '/Users/ravelarvargas/Downloads/marian/cs_cq_all_study_processed_deseq2.rds')

#No batch correction
cs_cq_no_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_cq_studies_processed_deseq[['deseq_obj']],
  batch_col=NULL,
  protect_col=NULL,
  method=NULL,
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

cs_cq_no_correction_pca_counts=cs_cq_no_correction_pca_counts
cs_cq_no_correction_pca_counts$cell_substate=as.character(cs_cq_no_correction_pca_counts$cell_substate)
cs_cq_no_correction_pca_counts$cell_substate[cs_cq_no_correction_pca_counts$cell_substate=='Oncogene-induced CS']='OIS'
cs_cq_no_correction_pca_counts$cell_substate[cs_cq_no_correction_pca_counts$cell_substate=='Stress-induced CS']='SIPS'
cs_cq_no_correction_pca_counts$cell_substate[cs_cq_no_correction_pca_counts$cell_substate=='Replicative CS']='RS'
cs_cq_no_correction_pca_counts$cell_substate[cs_cq_no_correction_pca_counts$cell_substate=='Serum-starved CQ']='SSCQ'
cs_cq_no_correction_pca_counts$cell_substate[cs_cq_no_correction_pca_counts$cell_substate=='Contact-inhibited CQ']='CICQ'
cs_cq_no_correction_pca_counts$cell_substate=factor(cs_cq_no_correction_pca_counts$cell_substate,
                                                 levels=c('Proliferating',
                                                          'CICQ','SSCQ','RS','SIPS','OIS'))
##pca
cs_cq_pca_no_correction_1=pca_rds(
  cs_cq_no_correction_pca_counts,
  ntop=500,
  pca_plot=c('cell_state','cell_substate','cell_line',
             'cell_type','tissue','cell_state','cell_state','cell_state',
             'sra.platform_model'),nrow=3)

cs_cq_pca_no_correction_2=pca_rds(
  cs_cq_no_correction_pca_counts,
  ntop=500,
  pca_plot=c('cell_substate'))


cs_cq_pca_no_correction_study=custom_pca(count_data=assay(cs_cq_no_correction_pca_counts),
                                         ntop=500,
                                         pheno_data=colData(cs_cq_no_correction_pca_counts),
                                         intgroup='study',hide_legend = TRUE)

#Batch correction
cs_cq_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_cq_studies_processed_deseq[['deseq_obj']],
  batch_col='study',
  protect_col=c('cell_substate'),
  method='wgcna',
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

cs_cq_correction_pca_counts=cs_cq_correction_pca_counts
cs_cq_correction_pca_counts$cell_substate=as.character(cs_cq_correction_pca_counts$cell_substate)
cs_cq_correction_pca_counts$cell_substate[cs_cq_correction_pca_counts$cell_substate=='Oncogene-induced CS']='OIS'
cs_cq_correction_pca_counts$cell_substate[cs_cq_correction_pca_counts$cell_substate=='Stress-induced CS']='SIPS'
cs_cq_correction_pca_counts$cell_substate[cs_cq_correction_pca_counts$cell_substate=='Replicative CS']='RS'
cs_cq_correction_pca_counts$cell_substate[cs_cq_correction_pca_counts$cell_substate=='Serum-starved CQ']='SSCQ'
cs_cq_correction_pca_counts$cell_substate[cs_cq_correction_pca_counts$cell_substate=='Contact-inhibited CQ']='CICQ'
cs_cq_correction_pca_counts$cell_substate=factor(cs_cq_correction_pca_counts$cell_substate,
                             levels=c('Proliferating',
                                      'CICQ','SSCQ','RS','SIPS','OIS'))

cs_cq_correction_pca_1=pca_rds(
  cs_cq_correction_pca_counts,
  pca_plot=c('cell_state','cell_substate','cell_line',
             'cell_type','tissue','cell_state','cell_state','cell_state',
             'sra.platform_model'),
  ntop=500,nrow = 3)

cs_cq_correction_pca_2=pca_rds(
  cs_cq_correction_pca_counts,
  pca_plot=c('cell_substate'),
  ntop=500)

study_correction=custom_pca(count_data=assay(cs_cq_correction_pca_counts),
           ntop=500,
           pheno_data=colData(cs_cq_correction_pca_counts),
           intgroup='study',hide_legend = TRUE)

save_p(plot = study_correction,
       file_name = 'cs_cq_batch_correction_study',save_dir = save_dir_figure_si,
       p_width = 2.5,p_height = 2.5)

save_p(plot = cs_cq_pca_no_correction_1,
       file_name = 'cs_cq_no_batch_correction_1',save_dir = save_dir_figure_si,
       p_width = 10,p_height = 10)
save_p(plot = cs_cq_pca_no_correction_2,
       file_name = 'cs_cq_no_batch_correction_2',
       save_dir = save_dir_figure_si,p_height=10)
save_p(plot = cs_cq_correction_pca_1,
       file_name = 'cs_cq_batch_correction_wgcna_1',save_dir = save_dir_figure_si,
       p_width = 13,p_height = 13)
save_p(plot = cs_cq_correction_pca_2,
       file_name = 'cs_cq_batch_correction_wgcna_2',save_dir = save_dir_figure_si,
       p_height = 10)
save_p(plot = cs_cq_correction_pca_manuscript,
       file_name = 'cs_cq_batch_correction_wgcna_all',save_dir = save_dir_figure_si,
       p_width = 10,p_height = 10)
save_p(plot = cs_cq_pca_no_correction_study,
       file_name = 'cs_cq_no_batch_correction_study',save_dir = save_dir_figure_si,
       p_width = 2.5,p_height = 2.5)

#find degs
rs_degs=find_degs_new(rse_obj=cs_cq_all_study_processed,
                      group_col='cell_substate',
                      group_1='Replicative CS',
                      group_2='Proliferating',
                      batch_col='study')
ois_degs=find_degs_new(rse_obj=cs_cq_all_study_processed,
                      group_col='cell_substate',
                      group_1='Oncogene-induced CS',
                      group_2='Proliferating',
                      batch_col='study')
sips_degs=find_degs_new(rse_obj=cs_cq_all_study_processed,
                      group_col='cell_substate',
                      group_1='Stress-induced CS',
                      group_2='Proliferating',
                      batch_col='study')
cq_ci_degs=find_degs_new(rse_obj=cs_cq_all_study_processed,
                      group_col='cell_substate',
                      group_1='Contact_inhibited CQ',
                      group_2='Proliferating',
                      batch_col='study')
cq_ss_degs=find_degs_new(rse_obj=cs_cq_all_study_processed,
                         group_col='cell_substate',
                         group_1='Serum_starved CQ',
                         group_2='Proliferating',
                         batch_col='study')

arrest_degs=rbind(rs_degs$degs,
                  ois_degs$degs,
                  sips_degs$degs,
                  cq_ci_degs$degs,
                  cq_ss_degs$degs)

save_csv(arrest_degs,file_name ='arrest_degs_final.csv',
         path = save_dir_csv)

arrest_degs=read.csv(paste0(save_dir_csv,'arrest_degs_final.csv'))
#heatmap
cs_cq_studies_processed_deseq$deseq_obj$cell_substate=as.character(cs_cq_studies_processed_deseq$deseq_obj$cell_substate)
cs_cq_studies_processed_deseq$deseq_obj$cell_substate=gsub(cs_cq_studies_processed_deseq$deseq_obj$cell_substate,pattern = '_',replacement = '-')
all_heatmap_15=build_heatmap(degs = arrest_degs,
                             dds_obj = cs_cq_studies_processed_deseq[['deseq_obj']],
                             group_col = 'cell_substate',
                             groups = c('Stress-induced CS',
                                        'Oncogene-induced CS',
                                        'Replicative CS',
                                        'Serum-starved CQ',
                                        'Contact-inhibited CQ',
                                        'Proliferating'),
                             protect_col = 'cell_substate',
                             batch_col = 'study',
                             group_col_heatmap = c('tissue','cell_substate',
                                                   'cell_state'),
                             normalisation = 'vst_full',
                             top_x = 15,
                             gene_col = 'gene',
                             rownames_to_symbol = FALSE)
# 
save_pheatmap(pheatmap_to_save=all_heatmap_15,
              file_name='all_arrest_condition_top_15',
              save_dir=save_dir_figure_si,
              p_width=5000,
              p_height=5000)

all_heatmap=build_heatmap(degs = arrest_degs,
                          dds_obj = cs_cq_studies_processed_deseq[['deseq_obj']],
                          group_col = 'cell_substate',
                          groups = c('Stress-induced CS',
                                     'Oncogene-induced CS',
                                     'Replicative CS',
                                     'Serum-starved CQ',
                                     'Contact-inhibited CQ',
                                     'Proliferating'),
                          protect_col = 'cell_substate',
                          batch_col = 'study',
                          group_col_heatmap = c('tissue','cell_substate',
                                                'cell_state'),
                          normalisation = 'vst_full',
                          show_rownames = FALSE,top_x = NULL,
                          gene_col = 'gene',
                          rownames_to_symbol = FALSE)
# 
save_pheatmap(pheatmap_to_save=all_heatmap,
              file_name='all_arrest_condition_all_deg',
              save_dir=save_dir_figure_si,
              p_width=5000,
              p_height=5000)

#volcano plots
rs=arrest_degs[arrest_degs$group_1%in%c('Replicative_CS'),]
sips=arrest_degs[arrest_degs$group_1%in%c('Stress_induced_CS'),]
ois=arrest_degs[arrest_degs$group_1%in%c('Oncogene_induced_CS'),]
cq_serum=arrest_degs[arrest_degs$group_1%in%c('Serum_starved_CQ'),]
cq_contact=arrest_degs[arrest_degs$group_1%in%c('Contact_inhibited_CQ'),]

rs_plot=EnhancedVolcano(rs,
                        lab = rs$gene,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = 'Replicative CS',subtitle = NULL)+
  theme(legend.position = "none")
sips_plot=EnhancedVolcano(sips,
                          lab = sips$gene,
                          x = 'log2FoldChange',
                          y = 'padj',
                          title = 'Stress-induced CS',subtitle = NULL)+
  theme(legend.position = "none")
ois_plot=EnhancedVolcano(ois,
                         lab = ois$gene,
                         x = 'log2FoldChange',
                         y = 'padj',
                         title = 'Oncogene-induced CS',
                         subtitle = NULL,
                         legendPosition = 'right')

cq_serum_plot=EnhancedVolcano(cq_serum,
                              lab = cq_serum$gene,
                              x = 'log2FoldChange',
                              y = 'padj',
                              title = 'Serum-starved CQ',
                              subtitle = NULL)+
  theme(legend.position = "none")
cq_contact_plot=EnhancedVolcano(cq_contact,
                                lab = cq_contact$gene,
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Contact-inhibited CQ',
                                subtitle = NULL)+
  theme(legend.position = "none")

volcano_plots=(rs_plot+sips_plot+ois_plot)/
  (cq_serum_plot+cq_contact_plot)

ggsave(volcano_plots,filename = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/volcano_plots.png',
       device='png',width=15,height=15)

#pi scores
arrest_degs_sig=arrest_degs[arrest_degs$sig=='y',]

all_pi_scores_sig=pi_score_filter_by_group(sig_degs = arrest_degs_sig,
                         group_col = 'dir_accession')

save_csv(all_pi_scores_sig,file_name ='pi_scores_sig_final',
         path = save_dir_csv)
#####
save_dir_csv_cq=paste0(save_dir_csv,'cq_test/')
dir.create(save_dir_csv_cq)
save_dir_figure_si_cq=paste0(save_dir_figure_si,'cq_test/')
dir.create(save_dir_figure_si_cq)
save_dir_figure_cq=paste0(save_dir_figure,'cq_test/')
dir.create(save_dir_figure_cq)
# cs vs cq
meta_both_cq_test=meta_both
meta_both_cq_test=meta_both_cq_test[meta_both_cq_test$cq_test!='Proliferating',]

filter_dataframe <- function(dataframe,
                             iterate_col,
                             filter_col,
                             col_label_1,
                             col_label_2,
                             col_1_n_filter = 0,
                             col_2_n_filter = 0) {
  dataframe=dataframe[grepl(dataframe[[filter_col]],pattern = col_label_1)|
                        grepl(dataframe[[filter_col]],pattern = col_label_2),]
  # Function to apply to each unique iterate_col value
  filter_val=do.call(c,lapply(unique(dataframe[[iterate_col]]),function(value) {
    subset_df <- dataframe[dataframe[[iterate_col]] == value,]
    
    counts_col_label_1 <- sum(grepl(subset_df[[filter_col]], pattern=col_label_1))>col_1_n_filter
    counts_col_label_2 <- sum(grepl(subset_df[[filter_col]], pattern=col_label_2))>col_2_n_filter
    
    if(counts_col_label_1&
       counts_col_label_2){
      return(value)
    }else{
      return(NULL)
    }
  }))
  
  # Apply filter_function to each unique value in iterate_col and combine results
  filtered_df <- dataframe[dataframe[[iterate_col]]%in%filter_val,]
  
  return(filtered_df)
}

cq_full_rank=filter_dataframe(dataframe=meta_both_cq_test,
                              iterate_col='sra',
                              filter_col='cq_test',
                              col_label_1='CQ',
                              col_label_2='CS')

cs_cq_all_study_processed_cq_test=process_rse(rse=cs_cq_download,
                                      meta_add=meta_both,
                                      meta_id_col = 'external_id',
                                      filter_pc=TRUE,
                                      ensembl_dictionary=human_pc,
                                      filter_duplicates=TRUE,
                                      sample_filter=cq_full_rank[['external_id']],
                                      low_expression_filter='per',
                                      scale_counts=TRUE)

study_table_cq=data.frame(table(data.frame(colData(cs_cq_all_study_processed_cq_test))%>%dplyr::select(tissue,
                                                                                            cq_test,
                                                                                            study)))

cell_counts_table_cq_test=data.frame(table(data.frame(colData(cs_cq_all_study_processed_cq_test))%>%dplyr::select(tissue, cq_test)))
save_csv(cell_counts_table_cq_test,file_name ='sample_subtype_counts.csv',
         path = save_dir_csv_cq)

temp_meta=data.frame(colData(cs_cq_all_study_processed_cq_test))
save_csv(temp_meta,file_name ='sample_metadata.csv',
         path = save_dir_csv_cq)

cs_cq_studies_processed_deseq_cq_test=deseq_studies(study_obj=cs_cq_all_study_processed_cq_test,
                                            protect_col='cq_test',
                                            fix_col='study')

#No batch correction
cs_cq_no_correction_pca_counts_cq_test=normalise_counts(
  dds_obj_use=cs_cq_studies_processed_deseq_cq_test[['deseq_obj']],
  batch_col=NULL,
  protect_col=NULL,
  method=NULL,
  normalisation='vst',
  vst_n=500,
  blind=TRUE)
##pca
cs_cq_pca_no_correction_cq_test=pca_rds(
  cs_cq_no_correction_pca_counts_cq_test,
  ntop=500,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cq_test',
             'sra.platform_model','study'))

#Batch correction
cs_cq_correction_pca_counts_cq_test=normalise_counts(
  dds_obj_use=cs_cq_studies_processed_deseq_cq_test[['deseq_obj']],
  batch_col='study',
  protect_col=c('cq_test'),
  method='wgcna',
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

cs_cq_correction_pca_cq_test=pca_rds(
  cs_cq_correction_pca_counts_cq_test,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cq_test',
             'sra.platform_model','study'),
  ntop=500)

save_p(plot = cs_cq_pca_no_correction_cq_test,
       file_name = 'cs_cq_no_batch_correction',save_dir = save_dir_figure_si_cq,
       p_width = 10,p_height = 10)
save_p(plot = cs_cq_correction_pca_cq_test,
       file_name = 'cs_cq_batch_correction_wgcna',save_dir = save_dir_figure_cq,
       p_width = 10,p_height = 10)
# save_p(plot = cs_cq_correction_pca_manuscript,
#        file_name = 'cs_cq_batch_correction_wgcna_all',save_dir = save_dir_figure_si_cq,
#        p_width = 10,p_height = 10)

cq_test_sips_degs=find_degs_new(rse_obj=cs_cq_all_study_processed_cq_test,
                         group_col='cq_test',
                         group_1=c('Stress-induced CS'),
                         group_2='CQ',
                         batch_col='study',
                         find_degs = TRUE)

cq_test_ois_degs=find_degs_new(rse_obj=cs_cq_all_study_processed_cq_test,
                                group_col='cq_test',
                                group_1=c('Oncogene-induced CS'),
                                group_2='CQ',
                                batch_col='study',
                                find_degs = TRUE)

cq_vs_cs_degs=rbind(cq_test_ois_degs$degs,
              cq_test_sips_degs$degs)

table(cq_vs_cs_degs%>%dplyr::filter(sig=='y')%>%dplyr::select(group_1,direction_1))

save_csv(cq_vs_cs_degs,file_name ='cq_vs_cs_degs.csv',
         path = save_dir_csv_cq)

#heatmap
arrest_degs=cq_vs_cs_degs
cs_cq_studies_processed_deseq_cq_test$deseq_obj$cell_substate=as.character(cs_cq_studies_processed_deseq_cq_test$deseq_obj$cell_substate)
cs_cq_studies_processed_deseq_cq_test$deseq_obj$cell_substate=gsub(cs_cq_studies_processed_deseq_cq_test$deseq_obj$cell_substate,pattern = '_',replacement = '-')
all_heatmap_15=build_heatmap(degs = arrest_degs,
                             dds_obj = cs_cq_studies_processed_deseq_cq_test[['deseq_obj']],
                             group_col = 'cell_substate',
                             groups = c('Stress-induced CS',
                                        'Oncogene-induced CS',
                                        'Serum-starved CQ',
                                        'Contact-inhibited CQ'),
                             protect_col = 'cell_substate',
                             batch_col = 'study',
                             group_col_heatmap = c('tissue','cell_substate',
                                                   'cell_state'),
                             normalisation = 'vst_full',
                             top_x = 15,
                             gene_col = 'gene',
                             rownames_to_symbol = FALSE)
# 
save_pheatmap(pheatmap_to_save=all_heatmap_15,
              file_name='all_arrest_condition_top_15',
              save_dir=save_dir_figure_cq,
              p_width=5000,
              p_height=5000)

all_heatmap=build_heatmap(degs = arrest_degs,
                          dds_obj = cs_cq_studies_processed_deseq_cq_test[['deseq_obj']],
                          group_col = 'cell_substate',
                          groups = c('Stress-induced CS',
                                     'Oncogene-induced CS',
                                     'Serum-starved CQ',
                                     'Contact-inhibited CQ'),
                          protect_col = 'cell_substate',
                          batch_col = 'study',
                          group_col_heatmap = c('tissue','cell_substate',
                                                'cell_state'),
                          normalisation = 'vst_full',
                          show_rownames = FALSE,top_x = NULL,
                          gene_col = 'gene',
                          rownames_to_symbol = FALSE)
# 
save_pheatmap(pheatmap_to_save=all_heatmap,
              file_name='all_arrest_condition_all_deg',
              save_dir=save_dir_figure_si_cq,
              p_width=5000,
              p_height=5000)
#####