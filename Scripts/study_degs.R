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

#Download studies
# cq_studies_full=download_studies(cq_studies)
# 
#Process studies
# cq_studies_processed=process_rse(rse=cq_studies_full,
#                                  meta_add=study_coldata_cq_filtered,
#                                  filter_pc=TRUE,
#                                  ensembl_dictionary=human_pc,
#                                  filter_duplicates=TRUE,
#                                  sample_filter=study_coldata_cq_filtered[['sample_ID']],
#                                  low_expression_filter=TRUE,
#                                  scale_counts=TRUE)
# 
# # Save study data
# study_coldata_cq=data.frame(colData(cq_studies_processed))
# save_csv(data = study_coldata_cq,file_name = 'study_info_cq',
#          path=save_dir_csv)
# 
# saveRDS(cq_studies_processed,paste0(save_dir_csv,'cq_samples.rds'))
cq_studies_processed=readRDS(paste0(save_dir_csv,'cq_samples.rds'))

#build model matrix
cq_studies_processed_deseq=deseq_studies(study_obj=cq_studies_processed,
                                         protect_col='cell_substate',
                                         fix_col='study')

#####
#pca
#No batch correction
cq_no_correction_pca_counts=normalise_counts(
  dds_obj_use=cq_studies_processed_deseq[['deseq_obj']],
  batch_col=NULL,
  protect_col=NULL,
  method=NULL,
  normalisation='vst',
  vst_n=500,
  blind=TRUE)
##pca
cq_pca_no_correction=pca_rds(
  cq_no_correction_pca_counts,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cell_substate',
             'sra.platform_model'),
  ntop=500)

save_p(plot = cq_pca_no_correction,
       file_name = 'cq_pca_no_correction',
       save_dir = save_dir_figure_si,p_width = 10,p_height=10)

cq_pca_no_correction_study=pca_rds(
  cq_no_correction_pca_counts,
  pca_plot=c('study'),
  ntop=500)

save_p(plot = cq_pca_no_correction_study,
       file_name = 'cq_pca_no_correction_study',
       save_dir = save_dir_figure_si,p_width = 6,p_height=8)

cq_correction_pca_counts=normalise_counts(
  dds_obj_use=cq_studies_processed_deseq[['deseq_obj']],
  batch_col='study',
  protect_col=c('cell_substate'),
  method='wgcna',
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

cq_correction_pca=pca_rds(
  cq_correction_pca_counts,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cell_substate',
             'sra.platform_model'),
  ntop=500)

save_p(plot = cq_correction_pca,
       file_name = 'cq_pca_correction',
       save_dir = save_dir_figure_si,p_width = 10,p_height=10)

cq_pca_correction_study=pca_rds(
  cq_correction_pca_counts,
  pca_plot=c('study'),
  ntop=500)

save_p(plot = cq_pca_correction_study,
       file_name = 'cq_pca_correction_study',
       save_dir = save_dir_figure_si,p_width = 6,p_height=8)

cq_correction_pca_manuscript=pca_rds(
  cq_correction_pca_counts,
  pca_plot = c('cell_state','cell_substate'),
  ntop=500)

save_p(plot = cq_correction_pca_manuscript,
       file_name = 'cq_pca',
       save_dir = save_dir_figure,p_width = 10,p_height=10)

#####
#DEGs

#CQ
cq_studies_processed$cell_substate[cq_studies_processed$cell_substate=='Serum_starved CQ']='Serum-starved CQ'
cq_studies_processed$cell_substate[cq_studies_processed$cell_substate=='Contact_inhibited CQ']='Contact-inhibited CQ'
##contact
contact_cq_degs=find_degs_new(rse_obj=cq_studies_processed,
                              group_col='cell_substate',
                              group_1='Contact-inhibited CQ',
                              group_2='Proliferating',
                              batch_col='study')

cq_contact_heatmap=build_heatmap(degs = contact_cq_degs$degs[contact_cq_degs$degs$sig=='y',],
                                 dds_obj = contact_cq_degs$dds_deseq$deseq_obj,
                                 group_col = 'cell_substate',
                                 groups = c('Contact-inhibited CQ','Proliferating'),
                                 protect_col = 'cell_substate',
                                 batch_col = 'study',
                                 normalisation = 'vst_full',
                                 group_col_heatmap=c('cell_substate',
                                                     'study','tissue'),
                                 gene_col = 'gene',
                                 rownames_to_symbol = FALSE)

cq_contact_heatmap_pi=pi_score_degs(sig_degs = contact_cq_degs$degs[contact_cq_degs$degs$sig=='y',])
save_csv(data = cq_contact_heatmap_pi$all_pi,
        file_name ='cq_contact_pi',path = save_dir_csv)

save_pheatmap(pheatmap_to_save=cq_contact_heatmap,
              file_name='contact_cq_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=3000,
              p_height=6000)

##serum
serum_cq_degs=find_degs_new(rse_obj=cq_studies_processed,
                            group_col='cell_substate',
                            group_1='Serum-starved CQ',
                            group_2='Proliferating',
                            batch_col='study')

cq_serum_heatmap=build_heatmap(degs = serum_cq_degs$degs[serum_cq_degs$degs$sig=='y',],
                               dds_obj = serum_cq_degs$dds_deseq$deseq_obj,
                               group_col = 'cell_substate',
                               groups = c('Serum-starved CQ','Proliferating'),
                               protect_col = 'cell_substate',
                               batch_col = 'study',
                               normalisation = 'vst_full',
                               gene_col = 'gene',
                               rownames_to_symbol = FALSE)
cq_serum_heatmap_pi=pi_score_degs(sig_degs = serum_cq_degs$degs[serum_cq_degs$degs$sig=='y',])
save_csv(data = cq_serum_heatmap_pi$all_pi,
         file_name ='cq_serum_pi',path = save_dir_csv)

save_pheatmap(pheatmap_to_save=cq_serum_heatmap,
              file_name='serum_cq_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=3000,
              p_height=6000)

#all cq
cq_degs=rbind(contact_cq_degs$degs,serum_cq_degs$degs)

dds_obj_cq=find_degs_new(rse_obj=cq_studies_processed,
                         group_col='cell_substate',
                         group_1=c('Contact-inhibited CQ',
                                   'Serum-starved CQ'),
                         group_2='Proliferating',
                         batch_col='study',
                         find_degs = FALSE)

cq_heatmap=build_heatmap(degs = cq_degs[cq_degs$sig=='y',],
                         dds_obj = dds_obj_cq$deseq_obj,
                         group_col = 'cell_substate',
                         groups = c('Serum-starved CQ',
                                    'Contact-inhibited CQ',
                                    'Proliferating'),
                         protect_col = 'cell_substate',
                         batch_col = 'study',
                         normalisation = 'vst_full',
                         group_col_heatmap = c('tissue','cell_substate',
                                               'cell_state'),
                         top_x = 15,
                         gene_col = 'gene',
                         rownames_to_symbol = FALSE)

save_pheatmap(pheatmap_to_save=cq_heatmap,
              file_name='cq_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=2000,
              p_height=2500)

save_csv(cq_degs,file_name = 'cq_degs',path = save_dir_csv)

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

cs_studies_full=download_studies(studies=cs_studies,
                                         sra_organism='human')

cs_studies_full_processed=process_rse(rse=cs_studies_full,
                                              meta_add=recount_cs,
                                              filter_pc=TRUE,
                                              ensembl_dictionary=human_pc,
                                              filter_duplicates=TRUE,
                                              sample_filter=recount_cs[['sample_ID']],
                                              low_expression_filter=TRUE,
                                              scale_counts=TRUE)

study_coldata_cs=data.frame(colData(cs_studies_full_processed))
save_csv(data = study_coldata_cs,
         file_name = 'study_info_cs',
         path = save_dir_csv)

cs_studies_full_processed_deseq=deseq_studies(study_obj=cs_studies_full_processed,
                                                      protect_col='cell_substate',
                                                      fix_col='study')

cs_no_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_studies_full_processed_deseq[['deseq_obj']],
  batch_col=NULL,
  protect_col=NULL,
  method=NULL,
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

#####
#pca
#No batch correction
cs_no_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_studies_full_processed_deseq[['deseq_obj']],
  batch_col=NULL,
  protect_col=NULL,
  method=NULL,
  normalisation='vst',
  vst_n=500,
  blind=TRUE)
##pca
cs_pca_no_correction=pca_rds(
  cs_no_correction_pca_counts,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cell_substate',
             'sra.platform_model'),
  ntop=500)

save_p(plot = cs_pca_no_correction,
       file_name = 'cs_pca_no_correction',
       save_dir = save_dir_figure_si,p_width = 10,p_height=10)

cs_pca_no_correction_study=pca_rds(
  cs_no_correction_pca_counts,
  pca_plot=c('study'),
  ntop=500)

save_p(plot = cs_pca_no_correction_study,
       file_name = 'cs_pca_no_correction_study',
       save_dir = save_dir_figure_si,p_width = 6,p_height=8)

cs_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_studies_full_processed_deseq[['deseq_obj']],
  batch_col='study',
  protect_col=c('cell_substate'),
  method='wgcna',
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

cs_correction_pca=pca_rds(
  cs_correction_pca_counts,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cell_substate',
             'sra.platform_model'),
  ntop=500)

save_p(plot = cs_correction_pca,
       file_name = 'cs_pca_correction',
       save_dir = save_dir_figure_si,p_width = 13,p_height=10)

cs_pca_correction_study=pca_rds(
  cs_correction_pca_counts,
  pca_plot=c('study'),
  ntop=500)

save_p(plot = cs_pca_correction_study,
       file_name = 'cs_pca_correction_study',
       save_dir = save_dir_figure_si,p_width = 6,p_height=8)

cs_correction_pca_manuscript=pca_rds(
  cs_correction_pca_counts,
  pca_plot = c('cell_state','cell_substate'),
  ntop=500)

save_p(plot = cs_correction_pca_manuscript,
       file_name = 'cs_pca',
       save_dir = save_dir_figure,p_width = 10,p_height=10)

#####
#DEGs
rs_degs=find_degs_new(rse_obj=cs_studies_full_processed,
                      group_col='cell_substate',
                      group_1='Replicative CS',
                      group_2='Proliferating',
                      batch_col='study')

rs_heatmap=build_heatmap(degs = rs_degs$degs[rs_degs$degs$sig=='y',],
                         dds_obj = rs_degs$dds_deseq$deseq_obj,
                         group_col = 'cell_substate',
                         groups = c('Replicative CS','Proliferating'),
                         protect_col = 'cell_substate',
                         batch_col = 'study',
                         normalisation = 'vst_full',
                         gene_col = 'gene',
                         rownames_to_symbol = FALSE)

rs_heatmap_pi=pi_score_degs(sig_degs = rs_degs$degs[rs_degs$degs$sig=='y',])
save_csv(data =rs_heatmap_pi$all_pi,
        file_name='rs_contact_pi',path=save_dir_csv)

save_pheatmap(pheatmap_to_save=rs_heatmap,
              file_name='rs_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=4000,
              p_height=5000)

ois_degs=find_degs_new(rse_obj=cs_studies_full_processed,
                       group_col='cell_substate',
                       group_1='Oncogene-induced CS',
                       group_2='Proliferating',
                       batch_col='study')

ois_heatmap=build_heatmap(degs = ois_degs$degs[ois_degs$degs$sig=='y',],
                          dds_obj = ois_degs$dds_deseq$deseq_obj,
                          group_col = 'cell_substate',
                          groups = c('Oncogene-induced CS','Proliferating'),
                          protect_col = 'cell_substate',
                          batch_col = 'study',
                          normalisation = 'vst_full',
                          gene_col = 'gene',
                          rownames_to_symbol = FALSE)
ois_heatmap_pi=pi_score_degs(sig_degs = ois_degs$degs[ois_degs$degs$sig=='y',])
save_csv(data =ois_heatmap_pi$all_pi,
        file_name='ois_pi',path=save_dir_csv)

save_pheatmap(pheatmap_to_save=ois_heatmap,
              file_name='ois_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=4000,
              p_height=5000)

sips_degs=find_degs_new(rse_obj=cs_studies_full_processed,
                        group_col='cell_substate',
                        group_1='Stress-induced CS',
                        group_2='Proliferating',
                        batch_col='study')

sips_heatmap=build_heatmap(degs = sips_degs$degs[sips_degs$degs$sig=='y',],
                           dds_obj = sips_degs$dds_deseq$deseq_obj,
                           group_col = 'cell_substate',
                           groups = c('Stress-induced CS','Proliferating'),
                           protect_col = 'cell_substate',
                           batch_col = 'study',
                           normalisation = 'vst_full',
                           gene_col = 'gene',
                           rownames_to_symbol = FALSE)
sips_heatmap_pi=pi_score_degs(sig_degs = sips_degs$degs[sips_degs$degs$sig=='y',])
save_csv(data =sips_heatmap_pi$all_pi,
        file_name='sips_pi',path=save_dir_csv)

save_pheatmap(pheatmap_to_save=sips_heatmap,
              file_name='sips_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=4000,
              p_height=5000)

cs_degs=rbind(rs_degs$degs,ois_degs$degs,sips_degs$degs)

dds_obj_cs=find_degs_new(rse_obj=cs_studies_full_processed,
                         group_col='cell_substate',
                         group_1=c('Stress-induced CS',
                                   'Oncogene-induced CS',
                                   'Replicative CS'),
                         group_2='Proliferating',
                         batch_col='study',
                         find_degs = FALSE)

cs_heatmap=build_heatmap(degs = cs_degs[cs_degs$sig=='y',],
                         dds_obj = dds_obj_cs$deseq_obj,
                         group_col = 'cell_substate',
                         groups = c('Stress-induced CS',
                                    'Oncogene-induced CS',
                                    'Replicative CS',
                                    'Proliferating'),
                         protect_col = 'cell_substate',
                         batch_col = 'study',
                         normalisation = 'vst_full',
                         group_col_heatmap = c('tissue','cell_substate',
                                               'cell_state'),
                         top_x = 15,
                         gene_col = 'gene',
                         rownames_to_symbol = FALSE)

save_pheatmap(pheatmap_to_save=cs_heatmap,
              file_name='cs_vs_proliferative_heatmap',
              save_dir=save_dir_figure_si,
              p_width=4000,
              p_height=4000)

save_csv(cs_degs,file_name = 'cs_degs',path = save_dir_csv)

#####
#both
meta_both=rbind(study_coldata_cq_filtered,
                recount_cs)
# 
# cs_cq=list(c(meta_both$sra),meta_both$sample_ID)
# names(cs_cq)=c('study','sample')

cs_cq_download=download_studies(studies=unique(c(cs_studies,
                                          cq_studies)),
                                sra_organism='human')

cs_cq_all_study_processed=process_rse(rse=cs_cq_download,
                                      meta_add=meta_both,
                                      filter_pc=TRUE,
                                      ensembl_dictionary=human_pc,
                                      filter_duplicates=TRUE,
                                      sample_filter=meta_both[['sample_ID']],
                                      low_expression_filter=TRUE,
                                      scale_counts=TRUE)

cell_counts_table=data.frame(table(data.frame(colData(cs_cq_all_study_processed))%>%dplyr::select(tissue, cell_substate)))
save_csv(cell_counts_table,file_name ='sample_subtype_counts.csv',
         path = save_dir_csv)

temp_meta=data.frame(colData(cs_cq_all_study_processed))
save_csv(temp_meta,file_name ='sample_metadata.csv',
         path = save_dir_csv)

cs_cq_studies_processed_deseq=deseq_studies(study_obj=cs_cq_all_study_processed,
                                            protect_col='cell_substate',
                                            fix_col='study')

#No batch correction
cs_cq_no_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_cq_studies_processed_deseq[['deseq_obj']],
  batch_col=NULL,
  protect_col=NULL,
  method=NULL,
  normalisation='vst',
  vst_n=500,
  blind=TRUE)
##pca
cs_cq_pca_no_correction=pca_rds(
  cs_cq_no_correction_pca_counts,
  ntop=500,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cell_substate',
             'sra.platform_model'))

cs_cq_pca_no_correction_study=pca_rds(
  cs_cq_no_correction_pca_counts,
  ntop=500,
  pca_plot=c('study'))

#Batch correction
cs_cq_correction_pca_counts=normalise_counts(
  dds_obj_use=cs_cq_studies_processed_deseq[['deseq_obj']],
  batch_col='study',
  protect_col=c('cell_substate'),
  method='wgcna',
  normalisation='vst',
  vst_n=500,
  blind=TRUE)

cs_cq_correction_pca=pca_rds(
  cs_cq_correction_pca_counts,
  pca_plot=c('cell_line',
             'cell_type','tissue','cell_state','cell_substate',
             'sra.platform_model'),
  ntop=500)

cs_cq_correction_pca_manuscript=pca_rds(
  cs_cq_correction_pca_counts,
  pca_plot = c('cell_state','cell_substate'),
  ntop=500)

save_p(plot = cs_cq_pca_no_correction,
       file_name = 'cs_cq_no_batch_correction',save_dir = save_dir_figure_si,
       p_width = 10,p_height = 10)
save_p(plot = cs_cq_correction_pca,
       file_name = 'cs_cq_batch_correction_wgcna',save_dir = save_dir_figure,
       p_width = 10,p_height = 10)
save_p(plot = cs_cq_correction_pca_manuscript,
       file_name = 'cs_cq_batch_correction_wgcna_all',save_dir = save_dir_figure_si,
       p_width = 10,p_height = 10)

#heatmap
arrest_degs=rbind(cs_degs,cq_degs)
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
              save_dir=save_dir_figure,
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
