library(biomaRt)
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')

dir_use='/Users/ravelarvargas/Downloads/marian/final/macrophage/'
senmayo_mouse=read.csv(paste0(dir_use,'senmayo_mouse.csv'))

macrophage_temp=read.csv(paste0(dir_use,'activated_macrophage_degs.csv'))

macrophage_temp$dir=ifelse(macrophage_temp$logFC>0,'Up','Down')

macrophage_temp_filtered=macrophage_temp[abs(macrophage_temp$logFC)>log2(1.5),]

ensembl100=useMart(host='http://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='mmusculus_gene_ensembl')

mouse_pc=getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
               filters = 'biotype',
               values = c('protein_coding'),
               mart = ensembl100)

mouse_pc_homology=getBM(attributes=c('external_gene_name','ensembl_gene_id','hsapiens_homolog_associated_gene_name',
                                     'hsapiens_homolog_orthology_type','hsapiens_homolog_orthology_confidence'),
                                     filters = 'biotype',
                                     values = c('protein_coding'),
                                     mart = ensembl100)

mouse_pc_homology_full=mouse_pc_homology[!is.na(mouse_pc_homology$hsapiens_homolog_orthology_confidence)&
                                           mouse_pc_homology$hsapiens_homolog_orthology_type=='ortholog_one2one'&
                                           mouse_pc_homology$hsapiens_homolog_orthology_confidence==1,]

senmayo_mouse_homology=merge(senmayo_mouse,mouse_pc_homology_full,by.x='gene',by.y='external_gene_name')

genelists_human=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/1_cellage.csv')
genelists_human_senmayo=genelists_human[genelists_human$source=='Saul et al',]

sum(genelists_human_senmayo$gene%in%senmayo_mouse_homology$hsapiens_homolog_associated_gene_name)/nrow(genelists_human_senmayo)*100

macrophage_temp_filtered$cond_2=
  ifelse(grepl(pattern='LPS',macrophage_temp_filtered$condition),'Classically Activated (LPS+IFNγ)',
         'Alternatively Activated (IL-4)')

macrophage_temp_filtered$condition_dir=paste0(macrophage_temp_filtered$cond_2,' ',
                                              macrophage_temp_filtered$dir)

senmayo_mouse$db='SenMayo Mouse'
macrophage_vs_senmayo=overlap_function(
  df_1 = senmayo_mouse,
  df_2 = macrophage_temp_filtered,
  gene_col_1 = 'gene',
  gene_col_2 = 'Gene.name',
  carry_col_1 = NULL,
  carry_col_2 = c('cond_2','dir'),
  group_col_1 = 'db',
  group_col_2 = 'condition_dir',
  background = mouse_pc$external_gene_name)

macro_vs_senmayo=simple_overlap_plot(df = macrophage_vs_senmayo,
                    x = 'db',
                    y = 'condition_dir')

save_p(plot = macro_vs_senmayo,
       file_name = 'macro_vs_senmayo',
       save_dir = dir_use,
       p_width = 4,p_height = 2)

save_csv(data = macrophage_vs_senmayo,
         file_name = 'macro_vs_senmayo',
         path = dir_use)
