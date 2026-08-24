#post_publication analysis stress modules
source("R/config.R")
source("R/functions.R")

arrest_temp=read.csv(file.path(RERUN_DIR, 'arrest_degs_final_RERUN.csv'))
arrest_temp_sig=arrest_temp[arrest_temp$sig=='y',]

stress_new=read.csv('/Volumes/exist_image/datasets/stress/stress_pathways_39408822_human.csv')

overlap_stress_test=overlap_function(df_1 = arrest_temp_sig,
                 df_2 = stress_new,
                 gene_col_1 = 'gene',
                 gene_col_2 = 'gene',
                 group_col_1 = 'dir_accession',
                 group_col_2 = 'accession',
                 background = unique(arrest_temp$gene),
                 carry_col_1 = c('group_1','direction_1'),
                 carry_col_2=c('pathway','pathway_type'))

overlap_stress_test_relevant_stress=overlap_stress_test[overlap_stress_test$pathway_type=='stress',]
overlap_stress_test_relevant_aging=overlap_stress_test[overlap_stress_test$pathway_type=='longevity',]

sig_stress=unique(overlap_stress_test_relevant_stress$pathway[overlap_stress_test_relevant_stress$adj<0.05])
sig_age=unique(overlap_stress_test_relevant_aging$pathway[overlap_stress_test_relevant_aging$adj<0.05])

create_overlap_plot(overlap_stress_test_relevant_aging,
                    x = 'group_1',
                    y = 'pathway',
                    facet_col = 'direction_1',remove_nonsig = TRUE,x_tilt = 90)
