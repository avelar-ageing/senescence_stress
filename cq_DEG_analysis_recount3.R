#analyse CQ DEGs
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')
save_dir='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/'
save_dir_csv=paste0(save_dir,'SI_tables/')
save_dir_csv_cq=paste0(save_dir_csv,'cq_test/')
save_dir_figure_cq=paste0(save_dir,'SI_figures/cq_test/')
cq_vs_cs_degs=read.csv(paste0(save_dir_csv_cq,'cq_vs_cs_degs.csv'))
sips_sig_cq=cq_vs_cs_degs[cq_vs_cs_degs$group_1=='Stress_induced_CS'&
                            cq_vs_cs_degs$sig=='y',]
ois_sig_cq=cq_vs_cs_degs[cq_vs_cs_degs$group_1=='Oncogene_induced_CS'&
                            cq_vs_cs_degs$sig=='y',]

cq_vs_cs_degs$accession=paste0(cq_vs_cs_degs$group_1,' ',
                               cq_vs_cs_degs$direction_1)
temp_overlap=self_overlaps_full_grid(df = cq_vs_cs_degs[cq_vs_cs_degs$sig=='y',],
                        gene_col = 'gene',
                        background = unique(c(cq_vs_cs_degs$gene)),
                        group_col = 'accession',
                        carry_col = c('group_1','direction_1'))

temp_overlap$accession=gsub(temp_overlap$accession,pattern='Oncogene_induced_CS',replacement='OIS')
temp_overlap$accession.1=gsub(temp_overlap$accession.1,pattern='Oncogene_induced_CS',replacement='OIS')
temp_overlap$accession=gsub(temp_overlap$accession,pattern='Stress_induced_CS',replacement='SIPS')
temp_overlap$accession.1=gsub(temp_overlap$accession.1,pattern='Stress_induced_CS',replacement='SIPS')
temp_overlap$accession=gsub(temp_overlap$accession,pattern='up',replacement='Up vs CQ')
temp_overlap$accession.1=gsub(temp_overlap$accession.1,pattern='up',replacement='Up vs CQ')
temp_overlap$accession=gsub(temp_overlap$accession,pattern='down',replacement='Down vs CQ')
temp_overlap$accession.1=gsub(temp_overlap$accession.1,pattern='down',replacement='Down vs CQ')

temp_overlap$pval[temp_overlap$actual==0]=1
temp_overlap$adj=unlist(lapply(temp_overlap$pval,function(x)p.adjust(x,method = 'BH',n = 4)))
cq_self_overlaps_p=simple_overlap_plot(df = temp_overlap,
                                  x = 'accession',
                                  y = 'accession.1',
                                  self=TRUE,
                                  x_tilt=45,remove_nonsig = TRUE)

save_p(plot = cq_self_overlaps_p,
       file_name = 'cs_cq_degs_self',
       save_dir = save_dir_figure_cq,p_height = 3,p_width=6)
save_csv(data = temp_overlap,
         file_name = 'cs_cq_degs_self',
         path = save_dir_csv_cq)
#####
#vs sasp atlas
all_sasp_fixed=read.csv(paste0(save_dir_csv,'sasp_atlas_genes.csv'))
all_sasp_fixed$accession=paste0(all_sasp_fixed$arrest,' ',all_sasp_fixed$sasp_class)
sasp_accessions=unique(all_sasp_fixed$accession)
##filter for SASP genes within DEG data
all_sasp_fixed=all_sasp_fixed[all_sasp_fixed$gene%in%
                                unique(cq_vs_cs_degs$gene),]

sasp_v_cs=do.call('rbind',lapply(sasp_accessions,function(per_accession){
  temp_sasp=all_sasp_fixed[all_sasp_fixed$accession==per_accession,]
  temp_background=temp_sasp$gene
  temp_overlap=do.call('rbind',lapply(c('CS Secretion','CQ Secretion'),function(per_direction){
    temp_sasp_dir=temp_sasp[temp_sasp$direction==per_direction,]
    temp_sasp_dir_sig=temp_sasp_dir[temp_sasp_dir$sig,]
    temp_overlap=overlap_function(df_1 = temp_sasp_dir_sig,
                                  df_2 = cq_vs_cs_degs,
                                  gene_col_1 = 'gene',
                                  gene_col_2 = 'gene',
                                  group_col_1 = 'accession',
                                  carry_col_1 = c('sasp_class','direction'),
                                  group_col_2 = 'accession',
                                  carry_col_2 = c('group_1','direction_1'),
                                  background = temp_background)
    return(temp_overlap)
  }))
}))

sasp_v_cs$adj=p.adjust(sasp_v_cs$pval,method='BH')

sasp_v_cs$accession.1[sasp_v_cs$accession.1=="Oncogene_induced_CS up"]=
  'OIS vs CQ Up'
sasp_v_cs$accession.1[sasp_v_cs$accession.1=="Oncogene_induced_CS down"]=
  'OIS vs CQ Down'
sasp_v_cs$accession.1[sasp_v_cs$accession.1=="Stress_induced_CS up"]=
  'SIPS vs CQ Up'
sasp_v_cs$accession.1[sasp_v_cs$accession.1=="Stress_induced_CS down"]=
  'SIPS vs CQ Down'
sasp_recount3_p=create_overlap_plot(deg_db_overlap = sasp_v_cs,
                                    odds_column = 'odds',
                                    pval_col = 'adj',
                                    facet_col = 'direction',
                                    x = 'accession',
                                    y = 'accession.1',
                                    ylab = 'CS DEGs direction',
                                    xlab = 'SASP Atlas Profiles',x_tilt = 45,
                                    remove_nonsig = TRUE)

save_p(plot = sasp_recount3_p,
       file_name = 'cs_cq_degs_vs_sasp_atlas',
       save_dir = save_dir_figure_cq,p_height = 3,p_width=6)
save_csv(data = sasp_v_cs,
         file_name = 'cs_cq_degs_vs_sasp_atlas',
         path = save_dir_csv_cq)

#simplify
df=sasp_v_cs
df$log2odds=log2(df$odds)
df$fill=ifelse(df$adj>0.05,'Not Sig',
               ifelse(df$log2odds<0,'Sig Under-represented','Sig Over-represented'))

fill_manual=c("Not Sig" = "#f1f1f1", "Sig Over-represented" = "#f8766d",
              "Sig Under-represented" = "#00bfc4")

df$accession.1[df$accession.1=='OIS vs CQ Up']='OIS  U  CQ  D'
df$accession.1[df$accession.1=='OIS vs CQ Down']='OIS  D  CQ  U'
df$accession.1[df$accession.1=='SIPS vs CQ Down']='SIPS  D  CQ  U'
df$accession.1[df$accession.1=='SIPS vs CQ Up']='SIPS  U  CQ  D'
summary_cq=df%>%ggplot(aes(x=accession,y=accession.1,fill=fill))+
  geom_tile(colour='black')+
  scale_fill_manual(name=legend_title,
                    values = fill_manual)+
  facet_grid(~direction)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  # geom_hline(yintercept = 3.5,size=1)+
  # geom_vline(xintercept = 3.5,size=1)+
  ylab('CS DEGs direction')+
  xlab('SASP Atlas Profiles')+
  theme(
    # legend.position='bottom',
        strip.text = element_text(size = 12),   # Facet title text size
        axis.text = element_text(size = 10),    # Axis text size
        axis.title = element_text(size = 12),   # Axis title text size
        legend.text = element_text(size = 10),  # Legend text size
        legend.title = element_text(size = 12),  # Legend title text size
        axis.text.x=element_text(angle=45,hjust=1)
  )
save_p(summary_cq,
       file_name = 'cq_sasp_summary',
       save_dir = '/Users/ravelarvargas/Downloads/marian/simplified_overlaps',p_width = 7,
       p_height = 3.5)
#####
#read in stress responses
# Define the path to your file
result=read.csv('/Users/ravelarvargas/Downloads/marian/stress_results/stress_response_pathways.csv')
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
#overlap with hallmark stress
cs_cq_v_stress=overlap_function(df_1=result,
                 df_2=cq_vs_cs_degs[cq_vs_cs_degs$sig=='y',],
                 gene_col_1='genes',
                 gene_col_2='gene',
                 group_col_2=c('dir_accession'),
                 group_col_1=c('pathway'),
                 # carry_col_2=c('group','arrest'),
                 carry_col_2 = c('group_1','direction_1'),
                 background=unique(cq_vs_cs_degs$gene))

cs_cq_v_stress_plot=cs_cq_v_stress[cs_cq_v_stress$pathway%in%pathway_of_interest,]

save_csv(cs_cq_v_stress,'cs_cq_v_stress',
         path = save_dir_csv_cq)

cs_cq_v_stress_plot$group_1=gsub(cs_cq_v_stress_plot$group_1,pattern='_',
                                 replacement=' ')
cs_cq_v_stress_plot$pathway=gsub(cs_cq_v_stress_plot$pathway,pattern='_',
                                 replacement=' ')
cs_cq_v_stress_plot$direction_1=gsub(cs_cq_v_stress_plot$direction_1,pattern='up',
                                 replacement='Up in CS vs CQ')
cs_cq_v_stress_plot$direction_1=gsub(cs_cq_v_stress_plot$direction_1,pattern='down',
                                 replacement='Down in CS vs CQ')
cs_cq_v_stress_plot$group_1=gsub(cs_cq_v_stress_plot$group,pattern='Oncogene induced CS',
                                     replacement='OIS')
cs_cq_v_stress_plot$group_1=gsub(cs_cq_v_stress_plot$group,pattern='Stress induced CS',
                               replacement='SIPS')
stress_cq=create_overlap_plot(deg_db_overlap = cs_cq_v_stress_plot,
                    odds_column = 'odds',
                    pval_col = 'adj',
                    facet_col = 'direction_1',
                    x = 'group_1',
                    y = 'pathway',
                    x_tilt = 0,
                    xlab = 'CS DEGs',
                    ylab = 'Stress Pathways',
                    remove_nonsig = TRUE)

save_p(stress_cq,'cs_cq_v_stress',
         save_dir = save_dir_figure_cq,p_height = 4,p_width = 7)
