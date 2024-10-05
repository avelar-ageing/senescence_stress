#msigdb
# Load necessary libraries
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')

#protein-coding genes
ensembl100=useMart(host='http://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)

#####
# Define the path to MSigDB
file_path <- "/Users/ravelarvargas/Downloads/h.all.v2023.2.Hs.symbols.gmt"

# Read the file into a dataframe
data <- readLines(file_path)

# Create an empty dataframe to store the results
result <- data.frame(pathway = character(), genes = character(), stringsAsFactors = FALSE)

# Loop through each line of the file and process it
for (line in data) {
  # Split the line by tab
  parts <- strsplit(line, "\t")[[1]]
  
  # Extract the pathway and the genes
  pathway <- parts[1]
  genes <- parts[-c(1,2)]
  
  # Create a temporary dataframe for the current line
  temp_df <- data.frame(pathway = rep(pathway, length(genes)), genes = genes, stringsAsFactors = FALSE)
  
  # Bind the temporary dataframe to the result dataframe
  result <- bind_rows(result, temp_df)
}

pathway_of_interest=c(
  "HALLMARK TNFA SIGNALING VIA NFKB",
  "HALLMARK P53 PATHWAY",
  "HALLMARK MYC TARGETS",
  "HALLMARK MTORC1 SIGNALING",
  "HALLMARK MITOTIC SPINDLE",
  # "HALLMARK KRAS SIGNALING UP",
  # "HALLMARK KRAS SIGNALING DN",
  "HALLMARK INTERFERON GAMMA RESPONSE",
  "HALLMARK INTERFERON ALPHA RESPONSE",
  "HALLMARK INFLAMMATORY RESPONSE",
  "HALLMARK IL6 JAK STAT3 SIGNALING",
  "HALLMARK HYPOXIA",
  "HALLMARK G2M CHECKPOINT",
  "HALLMARK E2F TARGETS",
  "HALLMARK APOPTOSIS",
  "HALLMARK DNA REPAIR",
  # "HALLMARK PI3K AKT MTOR SIGNALING",
  # "HALLMARK PROTEIN SECRETION",
  # "HALLMARK REACTIVE OXYGEN SPECIES PATHWAY",
  # "HALLMARK UNFOLDED PROTEIN RESPONSE",
  # "HALLMARK UV RESPONSE DN",
  # "HALLMARK UV RESPONSE UP",
  'Lysosomal Genes'
)

result$pathway=gsub(result$pathway,pattern='_V1',replacement='')
result$pathway=gsub(result$pathway,pattern='_V2',replacement='')
result$pathway=gsub(result$pathway,pattern='_',replacement=' ')

#add lysosomal genes
autophagy=read.csv('/Users/ravelarvargas/Downloads/marian/temp-autophagy.csv')
lysosome_genes=autophagy[autophagy$class=='Lysosome',]
lysosome_genes$pathway='Lysosomal Genes'
lysosome_genes=lysosome_genes%>%dplyr::select(pathway,gene)
colnames(lysosome_genes)[2]='genes'

result=rbind(result,
             lysosome_genes)

result=result[result$genes%in%
                human_pc$external_gene_name,]

result=unique(result)

save_csv(result,file_name = 'stress_response_pathways',
         path = '/Users/ravelarvargas/Downloads/marian/stress_results')
#####
#read in senescence data
arrest_degs_merged=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/arrest_degs_final.csv')
arrest_degs_merged=factor_column_and_modify(df = arrest_degs_merged,
                                            column = 'group_1',
                                            old_list = c('Contact_inhibited_CQ','Serum_starved_CQ',
                                                         'Replicative_CS','Stress_induced_CS','Oncogene_induced_CS'),
                                            keyword = NULL)
arrest_degs_merged_sig=arrest_degs_merged[arrest_degs_merged$sig=='y',]

temporal_degs=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/all_time_degs.csv')

temporal_degs[['cell_type_time_point']]=paste0(temporal_degs[['cell_type']],'_',
                                               temporal_degs[['group_1']])

temporal_degs=factor_column_and_modify(temporal_degs,
                                       column = 'group_1',
                                       old_list = rev(c('4_days',
                                                        '10_days',
                                                        '20_days')),
                                       keyword = NULL)

temporal_degs_sig=temporal_degs[temporal_degs$sig=='y',]

temporal_degs_sig$dir_accession=paste0(temporal_degs_sig$cell_type,'_',
                                       temporal_degs_sig$dir_accession)
#####
#main overlaps
ora_main_stress=overlap_function(df_1=result,
                                 df_2=arrest_degs_merged_sig,
                                 gene_col_1='genes',
                                 gene_col_2='gene',
                                 group_col_2=c('dir_accession'),
                                 group_col_1=c('pathway'),
                                 # carry_col_2=c('group','arrest'),
                                 carry_col_2 = c('group_1','direction_1'),
                                 background=arrest_degs_merged$gene)

ora_main_stress_plot=ora_main_stress[ora_main_stress$pathway%in%pathway_of_interest,]

ora_main_stress_plot$direction_1[ora_main_stress_plot$direction_1=='up']='Up in Arrest'
ora_main_stress_plot$direction_1[ora_main_stress_plot$direction_1=='down']='Down in Arrest'
ora_main_stress_plot$group_1=as.character(ora_main_stress_plot$group_1)

ora_main_stress_plot$group_1[ora_main_stress_plot$group_1=="Replicative CS"]='RS'
ora_main_stress_plot$group_1[ora_main_stress_plot$group_1=="Oncogene Induced CS"]='OIS'
ora_main_stress_plot$group_1[ora_main_stress_plot$group_1=="Stress Induced CS"]='SIPS'
ora_main_stress_plot$group_1[ora_main_stress_plot$group_1=="Contact Inhibited CQ"]='CICQ'
ora_main_stress_plot$group_1[ora_main_stress_plot$group_1=="Serum Starved CQ"]='SSCQ'
ora_main_stress_plot$group_1=factor(ora_main_stress_plot$group_1,
                                    levels=c('CICQ','SSCQ','RS','SIPS','OIS'))

#simple
ora_main_stress_plot$accession=apply(ora_main_stress_plot,1,function(x){
  temp_dir=ifelse(x[['direction_1']]=='Up in Arrest','Up','Down')
  paste0(x[['group_1']],' ',temp_dir)
})
ora_main_stress_plot$accession=gsub(ora_main_stress_plot$accession,pattern='Down',
                                    replacement='D')
ora_main_stress_plot$accession=gsub(ora_main_stress_plot$accession,pattern='Up',
                                    replacement='U')
ora_main_stress_plot$accession=factor(ora_main_stress_plot$accession,
                                      levels=c("SSCQ D",
                                               "CICQ D",
                                               "OIS D",
                                               "SIPS D",
                                               "RS D",
                                               "SSCQ U",
                                               "CICQ U",
                                               "OIS U",
                                               "SIPS U",
                                               "RS U"))

simple_deg_v_pathway=summarise_overlaps(db = ora_main_stress_plot,
               x = 'pathway',
               y = 'accession',
               xlab='MSigDB Pathways',
               ylab = 'Arrest-DEGs')

save_p(simple_deg_v_pathway,file_name = '3a_deg_v_pathway',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
       p_width = 6.5,p_height = 4)

ora_main_stress_plot$direction_1[ora_main_stress_plot$direction_1=='Down in Arrest']=
  'Down in Arrest   '
ora_main_stress_plot$direction_1[ora_main_stress_plot$direction_1=='Up in Arrest']=
  'Up in Arrest   '
ora_main_stress_p=create_overlap_plot(deg_db_overlap = ora_main_stress_plot,
                                      odds_column = 'odds',
                                      pval_col = 'adj',
                                      facet_col = 'direction_1',
                                      x = 'group_1',
                                      y = 'pathway',
                                      x_tilt = 45,
                                      remove_nonsig = TRUE,xlab='Cell Cycle Arrest DEGs',
                                      ylab = 'MSigDB Pathways')

save_p(save_dir = '/Users/ravelarvargas/Downloads/marian/stress_results',
       ora_main_stress_p,
       file_name = 'hallmark_vs_arrest_degs',p_height = 5.5,p_width = 7.5)

save_csv(path = '/Users/ravelarvargas/Downloads/marian/stress_results',
         file_name = 'hallmark_vs_arrest_degs',
         data = ora_main_stress)
#####
#temporal overlap
temporal_degs_sig$cell_type_time_point_dir=paste0(temporal_degs_sig$cell_type_time_point,
                                                  '_',temporal_degs_sig$direction_1)
ora_temporal_stress=overlap_function(df_1=result,
                                     df_2=temporal_degs_sig,
                                     gene_col_1='genes',
                                     gene_col_2='ensembl',
                                     group_col_2 = 'cell_type_time_point_dir',
                                     group_col_1=c('pathway'),
                                     carry_col_2=c('cell_type','group_1','direction_1'),
                                     background=unique(temporal_degs$ensembl))
ora_temporal_stress_plot=ora_temporal_stress
ora_temporal_stress_plot$cell_type_time=paste0(ora_temporal_stress_plot$cell_type,'_',
                                               ora_temporal_stress_plot$group_1)

ora_temporal_stress_plot=factor_column_and_modify(ora_temporal_stress_plot,
                                                  column = 'cell_type_time',
                                                  old_list = c('fibroblast_4_days',
                                                               'fibroblast_10_days',
                                                               'fibroblast_20_days',
                                                               'keratinocyte_4_days',
                                                               'keratinocyte_10_days',
                                                               'keratinocyte_20_days',
                                                               'melanocyte_4_days',
                                                               'melanocyte_10_days',
                                                               'melanocyte_20_days'),
                                                  keyword = NULL)

save_csv(data = ora_temporal_stress,
         file_name = 'temporal_hallmarks_overlap',
         path = '/Users/ravelarvargas/Downloads/marian/stress_results')

ora_temporal_stress_plot=ora_temporal_stress_plot[ora_temporal_stress_plot$pathway%in%pathway_of_interest,]
ora_temporal_stress_plot$direction_1[ora_temporal_stress_plot$direction_1=='up']='Up in Arrest'
ora_temporal_stress_plot$direction_1[ora_temporal_stress_plot$direction_1=='down']='Down in Arrest'
temporal_stress_p=create_overlap_plot(deg_db_overlap = ora_temporal_stress_plot,
                                      odds_column = 'odds',
                                      pval_col = 'adj',
                                      facet_col = 'direction_1',
                                      x = 'cell_type_time',
                                      y = 'pathway',
                                      x_tilt = 45,remove_nonsig = TRUE)
save_p(temporal_stress_p,
       file_name = 'temporal_hallmarks_overlap',
       save_dir = '/Users/ravelarvargas/Downloads/marian/stress_results',
       p_width = 10,p_height = 6)

ora_temporal_stress_plot$summary_y=apply(ora_temporal_stress_plot,MARGIN = 1,function(x){
  temp_x=ifelse(grepl(x[['direction_1']],pattern='Up'),'Up','Down')
  paste0(x[['group_1']],' ',temp_x)
})

ora_temporal_stress_plot$summary_y=factor(ora_temporal_stress_plot$summary_y,
                                          levels=c("20 Days Down",
                                                   "10 Days Down",
                                                   "4 Days Down",
                                                   "20 Days Up",
                                                   "10 Days Up",
                                                   "4 Days Up"))

ora_temporal_stress_plot$cell_type[ora_temporal_stress_plot$cell_type=='fibroblast']='Fibroblast'
ora_temporal_stress_plot$cell_type[ora_temporal_stress_plot$cell_type=='keratinocyte']='Keratinocyte'
ora_temporal_stress_plot$cell_type[ora_temporal_stress_plot$cell_type=='melanocyte']='Melanocyte'

ora_temporal_stress_plot$summary_y=as.character(ora_temporal_stress_plot$summary_y)
ora_temporal_stress_plot$summary_y=gsub(ora_temporal_stress_plot$summary_y,pattern='Up',replacement='U')
ora_temporal_stress_plot$summary_y=gsub(ora_temporal_stress_plot$summary_y,pattern='Down',replacement='D')
ora_temporal_stress_plot$summary_y=factor(ora_temporal_stress_plot$summary_y,
                                          levels=c("20 Days D",
                                                   "10 Days D",
                                                   "4 Days D",
                                                   "20 Days U",
                                                   "10 Days U",
                                                   "4 Days U"))

summary_p_temporal_stress=summarise_overlaps_facet(db=ora_temporal_stress_plot,
                         x='pathway',
                         y='summary_y',
                         odds_col='odds',
                         pval_col='adj',
                         xlab='MSigDB Pathways',
                         ylab='Temporal DEGs',
                         add_line=TRUE,
                         legend_title=NULL,
                         label_cq=TRUE,
                         facet_col='cell_type',angle = 45)

save_p(summary_p_temporal_stress,file_name = '7a_simplified_overlaps_temporal_deg_v_pathway',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
       p_width = 6,p_height = 5.5)
#####
#gene lists
cs_gene_list=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/1_cellage.csv')

cs_gene_list$accession[cs_gene_list$accession=='Saul et al ']='SenMayo'
cs_gene_list$accession[cs_gene_list$accession=='Driver Inhibits']='CellAge Inhibit CS'
cs_gene_list$accession[cs_gene_list$accession=='Driver Induces']='CellAge Induce CS'
cs_gene_list$accession[cs_gene_list$accession=='Signature Underexpressed']='CellAge Down in RS'
cs_gene_list$accession[cs_gene_list$accession=='Signature Overexpressed']='CellAge Up in RS'
cs_gene_list$accession=
  gsub(cs_gene_list$accession,pattern='down',replacement='Down')
cs_gene_list$accession=
  gsub(cs_gene_list$accession,pattern='up',replacement='Up')

cs_gene_list=cs_gene_list[!grepl(cs_gene_list$accession,pattern='Tao'),]
cs_factor_level=c('CellAge Induce CS',
                  'CellAge Inhibit CS',
                  'CellAge Up in RS',
                  'CellAge Down in RS',
                  'SenMayo',
                  'Hernandez-Segura et al Up',
                  'Hernandez-Segura et al Down',
                  'Casella et al Up',
                  'Casella et al Down',
                  'Cherry et al Up',
                  'Cherry et al Down')

#CS gene list overlap
msig_cs=overlap_function(df_1=result,
                         df_2=cs_gene_list,
                         gene_col_1='genes',
                         gene_col_2='gene',
                         group_col_2=c('accession'),
                         group_col_1=c('pathway'),
                         background=arrest_degs_merged$gene)

msig_cs_plot=msig_cs[msig_cs$pathway%in%pathway_of_interest,]
msig_cs_plot$accession=factor(msig_cs_plot$accession,levels = rev(cs_factor_level))
stress_vs_genelists_p=simple_overlap_plot(df = msig_cs_plot,
                                          x = 'accession',
                                          y = 'pathway',
                                          add_almost_sig = FALSE,
                                          remove_nonsig = TRUE,
                                          x_tilt = 45,
                                          xlab = 'CS Genelists',
                                          ylab = 'MSigDB Pathways')

save_p(stress_vs_genelists_p,
       file_name = 'cs_genelist_hallmarks_overlap',
       save_dir = '/Users/ravelarvargas/Downloads/marian/stress_results',
       p_width = 10,p_height = 6)
save_csv(msig_cs,
         file_name = 'cs_genelist_hallmarks_overlap',
         path = '/Users/ravelarvargas/Downloads/marian/stress_results')

database_v_pathway=summarise_overlaps(db = msig_cs_plot,
                                        x = 'pathway',
                                        y = 'accession',
                                        xlab='MSigDB Pathways',
                                        ylab = 'CS Genelists',
                                      add_line = FALSE,label_cq = FALSE)
save_p(database_v_pathway,file_name = '3b_simplified_db_v_pathway',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
       p_width = 6.5,p_height = 4)
#####
#stress vs sasp
all_sasp_fixed=read.csv(paste0(save_dir_csv,'sasp_atlas_genes.csv'))
#overlap by accession with background depending on stress condition
all_sasp_fixed$accession=paste0(all_sasp_fixed$arrest,' ',all_sasp_fixed$sasp_class)
sasp_accessions=unique(all_sasp_fixed$accession)
#filter for pc-genes in ensembl v100

sasp_stress=do.call('rbind',lapply(sasp_accessions,function(per_accession){
  temp_sasp=all_sasp_fixed[all_sasp_fixed$accession==per_accession,]
  temp_background=temp_sasp$gene
  temp_overlap=do.call('rbind',lapply(c('CS Secretion','CQ Secretion'),function(per_direction){
    temp_sasp_dir=temp_sasp[temp_sasp$direction==per_direction,]
    temp_sasp_dir_sig=temp_sasp_dir[temp_sasp_dir$sig,]
    temp_overlap=overlap_function(df_1 = temp_sasp_dir_sig,
                                  df_2 = result,
                                  gene_col_1 = 'gene',
                                  gene_col_2 = 'genes',
                                  group_col_1 = 'accession',
                                  carry_col_1 = c('sasp_class','direction'),
                                  group_col_2 = 'pathway',
                                  background = temp_background)
    return(temp_overlap)
  }))
}))
sasp_stress$adj=p.adjust(sasp_stress$pval,method='BH')

pathway_of_interest_sasp=c(
  pathway_of_interest,
  "HALLMARK EPITHELIAL MESENCHYMAL TRANSITION",
  "HALLMARK ANGIOGENESIS",
  "HALLMARK COAGULATION",
  "HALLMARK COMPLEMENT"
)

sasp_stress_plot=sasp_stress[sasp_stress$pathway%in%pathway_of_interest_sasp,]

sasp_v_stress=create_overlap_plot(deg_db_overlap = sasp_stress_plot,
                                  odds_column = 'odds',
                                  pval_col = 'adj',
                                  facet_col = 'direction',
                                  # facet_2 = 'group_1',
                                  x = 'accession',
                                  xlab = 'SASP Atlas Profiles',
                                  ylab = 'Stress Pathways',
                                  y = 'pathway',x_tilt = 45,
                                  remove_nonsig = TRUE)

save_p(sasp_v_stress,file_name = 'sasp_vs_hallmark',
       save_dir = '/Users/ravelarvargas/Downloads/marian/stress_results',
       p_width = 8,p_height = 6.5)
save_csv(path = '/Users/ravelarvargas/Downloads/marian/stress_results',
         data = sasp_stress,file_name='sasp_vs_hallmark')

#simplified
sasp_vs_stress_simplified=summarise_overlaps_facet(db=sasp_stress_plot,
                                                   x='pathway',
                                                   y='accession',
                                                   odds_col='odds',
                                                   pval_col='adj',
                                                   xlab='MSigDB Pathways',
                                                   ylab='SASP Profiles',
                                                   add_line=FALSE,
                                                   legend_title=NULL,
                                                   label_cq=TRUE,
                                                   facet_col='direction',angle = 45)

save_p(sasp_vs_stress_simplified,file_name = '5b_sasp_v_pathway',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',
       p_height = 4.5,p_width=7.5)
#####
#self overlaps stress
self_overlap_temp=self_overlaps_full_grid(df = result,
                                          gene_col = 'genes',
                                          background = arrest_degs_merged$gene,
                                          group_col = 'pathway')

self_overlap_temp_p=self_overlap_temp[self_overlap_temp$pathway%in%pathway_of_interest&
                                        self_overlap_temp$pathway.1%in%pathway_of_interest,]

stress_self_p=simple_overlap_plot(df = self_overlap_temp_p,
                                  x = 'pathway',
                                  y = 'pathway.1',
                                  self=TRUE,
                                  x_tilt=45,remove_nonsig = TRUE)

save_p(stress_self_p,file_name = 'stress_self_overlaps',
       save_dir = '/Users/ravelarvargas/Downloads/marian/stress_results',
       p_width = 8,p_height = 6.5)

self_overlap_temp=self_overlap_temp[self_overlap_temp$pathway!=self_overlap_temp$pathway.1,]
save_csv(self_overlap_temp,file_name = 'stress_self_overlaps',
         path = '/Users/ravelarvargas/Downloads/marian/stress_results')

#####
arrest_degs_merged=factor_column_and_modify(df = arrest_degs_merged,
                                            column = 'group_1',
                                            old_list = c('Contact_inhibited_CQ','Serum_starved_CQ',
                                                         'Replicative_CS','Stress_induced_CS','Oncogene_induced_CS'),
                                            keyword = NULL)
arrest_degs_merged$group=NA
#key genes
key_inflammation=c(
  'IL6',
  'CXCL8',
  'PTGS2',
  'IL1B'
  # 'NFKB1'
)
key_apoptosis=c(
  # 'BAX',
  'CASP3',
  # 'FAS',
  'BAK1',
  'BBC3',
  'BCL2L1',
  'PMAIP1'
  # 'APAF1'
)

key_cell_cycle=c(
  'CDKN1A',
  'CDKN2A',
  # 'RB1',
  # 'CCNA2',
  # 'CCNB1',
  'CCNE1',
  'CDK2',
  'CDK4',
  'CDK6',
  'MDM2',
  'CCND1'
)

key_autophagy=c(
  'MAP1LC3B',
  # 'BECN1',
  # 'ATG5',
  'ULK1',
  'LAMP1'
)

stress_tfs=c('TP53',
             # 'FOXO1',
             'FOXO3',
             'FOXO4',
             # 'NFE2L2',
             'CEBPB',
             'EGR2',
             'MITF',
             'ATF4')

chromatin_architecture=c('LMNB1',
                         'HMGA1',
                         'HMGA2',
                         'SUV39H1',
                         'LMNA'
                         # 'HIRA',
                         # 'ASF1A',
                         # 'H2AFY',
                         # 'HP1'
)

# arrest_degs_merged=label_genes(arrest_degs_merged,
#                                new_column = 'group',
#                                gene_vector = c(hsps),
#                                gene_column = 'gene',
#                                gene_label = 'hsps')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(chromatin_architecture),
                               gene_column = 'gene',
                               gene_label = 'Chromatin Architecture')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(key_inflammation),
                               gene_column = 'gene',
                               gene_label = 'Inflammation')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(key_apoptosis),
                               gene_column = 'gene',
                               gene_label = 'Apoptosis')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(key_cell_cycle),
                               gene_column = 'gene',
                               gene_label = 'Cell Cycle')

# arrest_degs_merged=label_genes(arrest_degs_merged,
#                                new_column = 'group',
#                                gene_vector = c(key_dna_damage),
#                                gene_column = 'gene',
#                                gene_label = 'key_dna')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(key_autophagy),
                               gene_column = 'gene',
                               gene_label = 'Lysosome/Autophagy')

arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(stress_tfs),
                               gene_column = 'gene',
                               gene_label = 'TFs')

# hif_components=c('HIF1A',
#                  'ARNT',
#                  'EPAS1',
#                  'ARNT2',)
# arrest_degs_merged=label_genes(arrest_degs_merged,
#                      new_column = 'group',
#                      gene_vector = c(hif_components),
#                      gene_column = 'gene',
#                      gene_label = 'HIF')

arrest_degs_GOI=compare_expression(degs=arrest_degs_merged,
                   gene_col = 'gene',
                   group = 'group_1',
                   genes=c(key_inflammation,
                           key_apoptosis,
                           key_cell_cycle,
                           # key_dna_damage,
                           key_autophagy,
                           stress_tfs,
                           chromatin_architecture
                           # hif_components
                   ),
                   facet = 'group',
                   flip_axis = TRUE,
                   scale_temp = 2,fill_cap = 6,
                   tilt_x = TRUE,
                   text_size = 13,
                   legend_bottom = TRUE,
                   hide_legend=FALSE)

save_p(arrest_degs_GOI,
       file_name = 'arrest_degs_goi',
       save_dir = '/Users/ravelarvargas/Downloads/marian/final',p_width = 10,p_height = 6.5)

#temporal
temporal_degs$group=NA
# temporal_degs=label_genes(temporal_degs,
#                                new_column = 'group',
#                                gene_vector = c(hsps),
#                                gene_column = 'ensembl',
#                                gene_label = 'hsps')

temporal_degs=label_genes(temporal_degs,
                          new_column = 'group',
                          gene_vector = c(chromatin_architecture),
                          gene_column = 'ensembl',
                          gene_label = 'Chromatin Architecture')

temporal_degs=label_genes(temporal_degs,
                          new_column = 'group',
                          gene_vector = c(key_inflammation),
                          gene_column = 'ensembl',
                          gene_label = 'Inflammation')

temporal_degs=label_genes(temporal_degs,
                          new_column = 'group',
                          gene_vector = c(key_apoptosis),
                          gene_column = 'ensembl',
                          gene_label = 'Apoptosis')

temporal_degs=label_genes(temporal_degs,
                          new_column = 'group',
                          gene_vector = c(key_cell_cycle),
                          gene_column = 'ensembl',
                          gene_label = 'Cell Cycle')

# temporal_degs=label_genes(temporal_degs,
#                           new_column = 'group',
#                           gene_vector = c(key_dna_damage),
#                           gene_column = 'ensembl',
#                           gene_label = 'key_dna')

temporal_degs=label_genes(temporal_degs,
                          new_column = 'group',
                          gene_vector = c(key_autophagy),
                          gene_column = 'ensembl',
                          gene_label = 'Autophagy')

temporal_degs=label_genes(temporal_degs,
                          new_column = 'group',
                          gene_vector = c(stress_tfs),
                          gene_column = 'ensembl',
                          gene_label = 'TFs')

temporal_goi=compare_expression(degs=temporal_degs,
                   gene_col = 'ensembl',
                   group = 'group_1',
                   genes=c(key_inflammation,
                           key_apoptosis,
                           key_cell_cycle,
                           # key_dna_damage,
                           key_autophagy,
                           stress_tfs,
                           chromatin_architecture
                           # hif_components
                   ),
                   facet = 'cell_type',
                   facet2 = 'group',
                   flip_axis = FALSE,
                   scale_temp = 2,fill_cap = 6,
                   tilt_x = TRUE,
                   text_size = 13,
                   legend_bottom = TRUE,
                   hide_legend=FALSE)

save_p(temporal_goi,
       file_name = 'temporal_goi',
       save_dir = '/Users/ravelarvargas/Downloads/marian/final',
       p_width = 11,p_height = 5.5)

#check unique inflammation genes
result=read.csv('/Users/ravelarvargas/Downloads/marian/stress_results/stress_response_pathways.csv')
cs_gene_list_induces=cs_gene_list[cs_gene_list$accession=="Driver Inhibits",]
arrest_degs_merged_sig_up_cs=arrest_degs_merged[arrest_degs_merged$sig=='y'&
                                               arrest_degs_merged$direction_1=='down'&
                                                 arrest_degs_merged$group_1%in%
                                                 c('Replicative CS',
                                                   'Oncogene Induced CS',
                                                   'Stress Induced CS'),]
arrest_degs_merged_sig_up_relevant=arrest_degs_merged_sig_up_cs[arrest_degs_merged_sig_up_cs$gene%in%cs_gene_list_induces$gene,]
temp_table=data.frame(table(result$genes)[table(result$genes)==1])

temp_table$Var1=as.character(temp_table$Var1)
list_induces=temp_table[temp_table$Var1%in%arrest_degs_merged_sig_up_relevant$gene,]
result_induces=result[result$genes%in%list_induces$Var1,]
View(result_induces)
