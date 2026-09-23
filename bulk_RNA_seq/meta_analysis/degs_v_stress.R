#overlap with stress responses
#####
# cs_degs=read.csv(file.path(SAVE_DIR_CSV, 'cs_degs.csv'))
# cq_degs=read.csv(file.path(SAVE_DIR_CSV, 'cq_degs.csv'))
cs_signatures_studies=read.csv(file.path(SAVE_DIR_CSV, '1_cellage.csv'))
# cs_signatures_studies=read.csv('/Users/ravelarvargas/Downloads/CS_signatures_studies.csv',
#                                skip = 1)
# 
# cs_signatures_studies=cs_signatures_studies[,colnames(cs_signatures_studies)%in%
#                                               c('gene','dir','source')]
# colnames(cellage)=c('gene','dir','source')
# cs_signatures_studies=rbind(cs_signatures_studies,cellage)
# cs_signatures_studies$dir[is.na(cs_signatures_studies$dir)]=''
source("R/config.R")
source("R/functions.R")

ensembl100=useMart(host='https://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)

human_pc_entrez=getBM(attributes=c('external_gene_name','ensembl_gene_id',
                                   'entrezgene_id'),
                      filters = 'biotype',
                      values = c('protein_coding'),
                      mart = ensembl100)

# arrest_degs_merged=rbind(cs_degs,
#                          cq_degs)
arrest_degs_merged=read.csv(file.path(SAVE_DIR_CSV, 'arrest_degs_final.csv'))

# save_csv(cellage,file_name = '1_cellage',path = save_dir_csv)
arrest_degs_merged$arrest=ifelse(grepl(arrest_degs_merged$group_1,pattern = 'CS'),'CS','CQ')

arrest_degs_merged_sig=arrest_degs_merged[arrest_degs_merged$sig=='y',]
#####
# gene lists
# ##autophagy
# autophagy=read.csv('/Users/ravelarvargas/Downloads/post_thesis/autophagy_gene.csv')
# autophagy_pc=autophagy[autophagy$Official.Gene.symbol%in%human_pc$external_gene_name,]
# autophagy_pc$func=NA
# autophagy_pc$func[grepl(autophagy_pc$Biological.Function,pattern = 'Positive')]='Positive'
# autophagy_pc$func[grepl(autophagy_pc$Biological.Function,pattern = 'Negative')]='Negative'
# auto_lyso=autophagy_pc[autophagy_pc$Group=='lysosome',]
# auto_lyso$func='lyso'
# save_csv(auto_lyso,file_name = 'lyso_genes',path = save_dir_csv)
# autophagy_pc=autophagy_pc[!is.na(autophagy_pc$func),]
# 
# save_csv(autophagy_pc,file_name='autophagy_genes',path=save_dir_csv)
# 
# #Inflammation
# inflammation_genes=getBM(attributes=c('external_gene_name','go_id'),
#                          filters = 'go',
#                          values = c('GO:0050729','GO:0050728'),
#                          mart = ensembl100)
# 
# inflammation_genes=inflammation_genes[inflammation_genes$external_gene_name%in%human_pc$external_gene_name,]
# 
# inflammation_genes=inflammation_genes[inflammation_genes$go_id%in%
#                                         c('GO:0050729','GO:0050728'),]
# 
# inflammation_genes$go_term=ifelse(inflammation_genes$go_id=='GO:0050729',
#                                   'Promotes inflammatory response',
#                                   'Inhibits inflammatory response')
# 
# save_csv(data = inflammation_genes,file_name = 'go_inflamm',path = save_dir_csv)
# 
# #DNA damage
dna_damage=read.csv('/Users/ravelarvargas/Downloads/marian/DNA_repair.csv')
# 
# #apoptosis
# apoptosis_regultors=read.csv('/Users/ravelarvargas/Downloads/marian/apoptosis_regulators.csv')
# apoptosis_regultors=apoptosis_regultors[,colnames(apoptosis_regultors)!='source']
# 
# #scaps
# scaps=read.csv('/Users/ravelarvargas/Downloads/SCAPs.csv')

#####
arrest_degs_merged_sig$direction_1=as.character(arrest_degs_merged_sig$direction_1)
#Overlap with autophagy, lysosome, inflammation, and the SASP atlas
arrest_degs_merged_sig=factor_column_and_modify(df = arrest_degs_merged_sig,
                                                column = 'direction_1',
                                                old_list = c('down','up'),
                                                keyword = '\nin Arrest')

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
                                   ggtitle = 'Inflammation vs SASP Atlas',x_tilt = 45)
save_p(inflam_vs_sasp,
       file_name = 'inflam_vs_sasp',
       save_dir = save_dir_figure_si,p_width = 6,p_height = 3)
save_csv(sasp_inflam,file_name = 'inflam_sasp',path = save_dir_csv)

#DNA damage
dna_damage=dna_damage[dna_damage$keep,]
dna_damage$regulation='test'
dna_damage_recount3=overlap_function(df_1=arrest_degs_merged_sig,
                                    df_2=dna_damage,
                                    gene_col_1='gene',
                                    gene_col_2='gene',
                                    group_col_1=c('dir_accession'),
                                    carry_col_1 = c('group_1','direction_1'),
                                    group_col_2=c('regulation'),
                                    background=arrest_degs_merged$gene)
save_csv(dna_damage_recount3,file_name = 'dna_damage_recount_overlap',path = save_dir_csv)

dna_damage_recount3_p=create_overlap_plot(deg_db_overlap = dna_damage_recount3,
                                         odds_column = 'odds',
                                         pval_col = 'adj',
                                         facet_col = 'group_1',
                                         x = 'direction_1',
                                         y = 'regulation',
                                         xlab = 'Cell Cycle Arrest DEGs',
                                         ylab = 'DNA Damage Regulators',
                                         ggtitle = 'Arrest DEGs vs DNA Damage')
save_p(dna_damage_recount3_p,file_name = 'dna_damage_recount_overlap',
       save_dir = save_dir_figure,p_height = 2.5,p_width = 6)

#apoptosis
apoptosis_recount3=overlap_function(df_1=arrest_degs_merged_sig,
                                     df_2=apoptosis_regultors,
                                     gene_col_1='gene',
                                     gene_col_2='gene',
                                     group_col_1=c('group_1','direction_1'),
                                     group_col_2=c('regulation'),
                                     background=arrest_degs_merged$gene)
save_csv(apoptosis_recount3,file_name = 'apoptosis_recount_overlap',path = save_dir_csv)

apoptosis_recount3_p=create_overlap_plot(deg_db_overlap = apoptosis_recount3,
                                          odds_column = 'odds',
                                          pval_col = 'adj',
                                          facet_col = 'group_1',
                                          x = 'direction_1',
                                          y = 'regulation',
                                          xlab = 'Cell Cycle Arrest DEGs',
                                          ylab = 'DNA Damage Regulators',
                                          ggtitle = 'Arrest DEGs vs DNA Damage')
save_p(dna_damage_recount3_p,file_name = 'dna_damage_recount_overlap',
       save_dir = save_dir_figure,p_height = 2.5,p_width = 6)

#scaps
scap_recount3=overlap_function(df_1=arrest_degs_merged_sig,
                                    df_2=scaps,
                                    gene_col_1='gene',
                                    gene_col_2='gene',
                                    group_col_1=c('group_1','direction_1'),
                                    group_col_2=c('type'),
                                    background=arrest_degs_merged$gene)
#####
#stress v gene lists
# senmayo=cs_signatures_studies[cs_signatures_studies$source=='Saul et al',]

cs_signatures_studies_use=cs_signatures_studies[!cs_signatures_studies$source%in%
                                                  c('Driver','Signature'),]
cs_signatures_studies_use$dir[cs_signatures_studies_use$dir=='']='up'

autophagy_senmayo=overlap_function(df_1=cs_signatures_studies_use,
                                    df_2=autophagy_pc,
                                    gene_col_1='gene',
                                    gene_col_2='Official.Gene.symbol',
                                    group_col_1=c('source','dir'),
                                    group_col_2=c('func'),
                                    background=arrest_degs_merged$gene)
cs_list_v_autophagy=create_overlap_plot(deg_db_overlap = autophagy_senmayo,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'dir',
                    # facet_2 = 'group_1',
                    x = 'source',
                    y = 'func',
                    xlab = 'Gene List',
                    ylab = 'Autophagy Regulators',
                    ggtitle = 'Autophagy vs CS Gene List',x_tilt = 45)

save_p(cs_list_v_autophagy,
       file_name = 'auto_vs_cs_list',
       save_dir = save_dir_figure_si,p_width = 7,p_height = 3)
save_csv(autophagy_senmayo,file_name = 'auto_vs_cs_list',path = save_dir_csv)

#lyso
senmayo_lyso=overlap_function(df_1=cs_signatures_studies_use,
                               df_2 = auto_lyso,
                               gene_col_1='gene',
                               gene_col_2='Official.Gene.symbol',
                               group_col_1=c('source','dir'),
                               group_col_2=c('func'),
                               background=arrest_degs_merged$gene)

cs_list_vs_lyso=create_overlap_plot(deg_db_overlap = senmayo_lyso,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'dir',
                    # facet_2 = 'group_1',
                    x = 'source',
                    y = 'func',
                    xlab = 'Gene List',
                    ylab = 'Lysosome Regulators',
                    ggtitle = 'Lysosome vs CS Gene List',x_tilt = 45)
save_p(cs_list_vs_lyso,
       file_name = 'lyso_vs_cs_list',
       save_dir = save_dir_figure_si,p_width = 7,p_height = 3)
save_csv(senmayo_lyso,file_name = 'lyso_vs_cs_list',path = save_dir_csv)

#inflammation
senmayo_inflammation=overlap_function(df_1=cs_signatures_studies_use,
                                      df_2=inflammation_genes,
                                      gene_col_1='gene',
                                      gene_col_2='external_gene_name',
                                      group_col_1=c('source','dir'),
                                      group_col_2=c('go_term'),
                                      background=arrest_degs_merged$gene)

cs_list_vs_inflamm=create_overlap_plot(deg_db_overlap = senmayo_inflammation,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'dir',
                    # facet_2 = 'group_1',
                    x = 'source',
                    y = 'go_term',
                    xlab = 'Gene List',
                    ylab = 'Inflammation Regulators',
                    ggtitle = 'Inflammation vs CS Gene List',x_tilt = 45)

save_p(cs_list_vs_inflamm,
       file_name = 'inflam_vs_cs_list',
       save_dir = save_dir_figure_si,p_width = 7,p_height = 3)
save_csv(senmayo_inflammation,file_name = 'inflam_vs_cs_list',path = save_dir_csv)

#DNA damage
senmayo_dna_damage=overlap_function(df_1=cs_signatures_studies_use,
                                      df_2=dna_damage,
                                      gene_col_1='gene',
                                      gene_col_2='gene',
                                      group_col_1=c('source','dir'),
                                      group_col_2=c('regulation'),
                                      background=arrest_degs_merged$gene)

senmayo_dna_damage$x=paste0(senmayo_dna_damage$source,' ',
                            senmayo_dna_damage$dir)
senmayo_dna_damage$y='dna'
cs_list_vs_dna_damage=create_overlap_plot(deg_db_overlap = senmayo_dna_damage,
                                       odds_column = 'odds',
                                       pval_col = 'adj',
                                       facet_col = 'regulation',
                                       # facet_2 = 'group_1',
                                       x = 'x',
                                       y = 'y',
                                       xlab = 'Gene List',
                                       ylab = 'Inflammation Regulators',
                                       ggtitle = 'Inflammation vs CS Gene List',
                                       x_tilt = 45,remove_y = TRUE)

#apoptosis
senmayo_apoptosis=overlap_function(df_1=cs_signatures_studies_use,
                                    df_2=apoptosis_regultors,
                                    gene_col_1='gene',
                                    gene_col_2='gene',
                                    group_col_1=c('source','dir'),
                                    group_col_2=c('regulation'),
                                    background=arrest_degs_merged$gene)
senmayo_apoptosis$x=paste0(senmayo_apoptosis$source,' ',
                           senmayo_apoptosis$dir)
senmayo_apoptosis$y='dna'

cs_list_vs_apoptosis=create_overlap_plot(deg_db_overlap = senmayo_apoptosis,
                                          odds_column = 'odds',
                                          pval_col = 'adj',
                                          facet_col = 'regulation',
                                          # facet_2 = 'group_1',
                                          x = 'x',
                                          y = 'y',
                                          xlab = 'Gene List',
                                          ylab = 'Inflammation Regulators',
                                          ggtitle = 'Inflammation vs CS Gene List',x_tilt = 45,
                                         remove_y = TRUE)

#scaps
senmayo_scaps=overlap_function(df_1=cs_signatures_studies_use,
                                   df_2=scaps,
                                   gene_col_1='gene',
                                   gene_col_2='gene',
                                   group_col_1=c('source','dir'),
                                   group_col_2=c('type'),
                                   background=arrest_degs_merged$gene)
senmayo_scaps$x=paste0(senmayo_scaps$source,' ',
                           senmayo_scaps$dir)
senmayo_scaps$y='dna'

cs_list_vs_scaps=create_overlap_plot(deg_db_overlap = senmayo_scaps,
                                         odds_column = 'odds',
                                         pval_col = 'adj',
                                         facet_col = 'type',
                                         # facet_2 = 'group_1',
                                         x = 'x',
                                         y = 'y',
                                         xlab = 'Gene List',
                                         ylab = 'SCAPs',
                                         ggtitle = 'SCAPs vs CS Gene List',x_tilt = 45,
                                         remove_y = TRUE)

#Enrich senmayo
senmayo_temp=cs_signatures_studies_use[cs_signatures_studies_use$source=='Saul et al',]
# senmayo_temp=human_pc_entrez$entrezgene_id[human_pc_entrez$external_gene_name%in%
#                                              senmayo_temp$gene]
senmayo_enrich=enrich_genes(gene_list=senmayo_temp$gene,
                            background = unique(human_pc_entrez$external_gene_name),
             gene_dictionary=human_pc_entrez,
             use_ensembl=FALSE)

save_csv(senmayo_enrich$enrichment,'senmayo_enrichment',
         path =  RERUN_DIR)

senmayo_enrich$go_plot

#overlap SIDs
sids=read.csv('/Users/ravelarvargas/Downloads/sig_group.csv')

sid_degs_overlap=overlap_function(df_1 = sids,
                 df_2 = arrest_degs_merged_sig,
                 gene_col_1 = 'gene',
                 gene_col_2 = 'gene',
                 carry_col_2 =c('group_1','direction_1'),
                 group_col_2 = c('group_1','direction_1'),
                 group_col_1=c('SID'),
                 background = arrest_degs_merged_sig$gene)

create_overlap_plot(sid_degs_overlap,x = 'direction_1',y = 'SID',
                    facet_col = 'group_1')

sid_autophagy=overlap_function(df_1 = sids,
                                  df_2 = autophagy_pc,
                                  gene_col_1 = 'gene',
                                  gene_col_2 = 'Official.Gene.symbol',
                               group_col_1 = 'SID',
                                  group_col_2 ='func',
                                  background = arrest_degs_merged_sig$gene)

create_overlap_plot(sid_autophagy,x = 'SID',y = 'func',
                    facet_col = 'func',remove_y = TRUE)

sid_scaps=overlap_function(df_1 = sids,
                               df_2 = scaps[scaps$type=='SCAP',],
                               gene_col_1 = 'gene',
                               gene_col_2 = 'gene',
                               group_col_1 = 'SID',
                               group_col_2 ='type',
                               background = arrest_degs_merged_sig$gene)

sid_scaps=overlap_function(df_1 = sids,
                           df_2 = auto_lyso,
                           gene_col_1 = 'gene',
                           gene_col_2 = 'Official.Gene.symbol',
                           group_col_1 = 'SID',
                           group_col_2 ='func',
                           background = arrest_degs_merged_sig$gene)

sid_apoptosis=overlap_function(df_1 = sids,
                           df_2 = apoptosis_regultors,
                           gene_col_1 = 'gene',
                           gene_col_2 = 'gene',
                           group_col_1 = 'SID',
                           group_col_2 ='regulation',
                           background = arrest_degs_merged_sig$gene)

sid_inflam=overlap_function(df_1 = sids,
                               df_2 = inflammation_genes,
                               gene_col_1 = 'gene',
                               gene_col_2 = 'external_gene_name',
                               group_col_1 = 'SID',
                               group_col_2 ='go_term',
                               background = arrest_degs_merged_sig$gene)
