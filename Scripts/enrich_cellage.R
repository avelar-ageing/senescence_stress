
source('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/Scripts/for_github/general_functions.R')
#Read in CellAge
cellage=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/1_cellage.csv')
cellage_driver=cellage[cellage$Database=='Driver',]
#####
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')

human_pc_entrez=getBM(attributes=c('external_gene_name','entrezgene_id'),
                      filters = 'biotype',
                      values = c('protein_coding'),
                      mart = ensembl100)
human_pc_entrez=human_pc_entrez[!is.na(human_pc_entrez$entrezgene_id),]

save_csv(human_pc_entrez,file_name = 'enrichment_background',
         path = save_dir_csv)

cellage_enrichment_KEGG=enrich_genes(gene_list = cellage_driver$gene_name,
             background = human_pc_entrez$external_gene_name,
             gene_dictionary = human_pc_entrez,
             graph_enrichment = TRUE,use_ensembl = FALSE)

write.csv(cellage_enrichment_KEGG$enrichment,
          paste0(save_dir_csv,'enriched_kegg.csv'),row.names=FALSE)

save_p(cellage_enrichment_KEGG$kegg_plot,
       file_name = 'enriched_kegg',save_dir = save_dir_figure,
       p_width = 11,p_height = 10)

cellage_enrichment_GO=enrich_genes(gene_list = cellage_driver$gene_name,
                                     background = human_pc_entrez$external_gene_name,
                                     gene_dictionary = human_pc_entrez,
                                     graph_enrichment = FALSE)
save_csv(cellage_enrichment_GO,
         path = save_dir_csv,file_name = 'enriched_GO.csv',row.names=FALSE)

sem_sim=get_sem_sim()

rrvgo_treemap(ont_sem_sim = sem_sim,
              go_terms = cellage_enrichment_GO,
              file_name = 'enriched_GO',
              save_dir = save_dir_figure)

#####
#Enrich individual cellage databases
cellage_promotes=cellage[cellage$Senescence=='Induces',]
cellage_inhibits=cellage[cellage$Senescence=='Inhibits',]

cellage_promotes_enrichment=enrich_genes(gene_list = cellage_promotes$gene_name,
             background = cellage$gene_name,
             gene_dictionary = human_pc_entrez,
             graph_enrichment = FALSE,use_ensembl = FALSE)
cellage_promotes_enrichment=cellage_promotes_enrichment$enrichment
cellage_promotes_enrichment$group='Induce_cs'
cellage_promotes_enrichment_kegg=cellage_promotes_enrichment[cellage_promotes_enrichment$enrichment==
                                                               'KEGG',]
cellage_promotes_enrichment_kegg_plot=
  enrichment_dotplot(cellage_promotes_enrichment_kegg,plot_max = 40)
save_p(cellage_promotes_enrichment_kegg_plot,
       file_name = 'enriched_kegg_induce_cs',save_dir = save_dir_figure,
       p_width = 11,p_height = 9)

cellage_promotes_enrichment_go=cellage_promotes_enrichment[cellage_promotes_enrichment$enrichment=='GO',]
rrvgo_treemap(ont_sem_sim = sem_sim,
              go_terms = cellage_promotes_enrichment_go,
              file_name = 'enriched_GO_induce_cs',
              save_dir = save_dir_figure)

cellage_inhibits_enrichment=enrich_genes(gene_list = cellage_inhibits$gene_name,
                                         background = cellage$gene_name,
                                         gene_dictionary = human_pc_entrez,
                                         graph_enrichment = FALSE)
cellage_inhibits_enrichment$group='Inhibit_cs'
rrvgo_treemap(ont_sem_sim = sem_sim,
              go_terms = cellage_inhibits_enrichment,
              file_name = 'enriched_GO_inhibit_cs',
              save_dir = save_dir_figure)

cellage_enrichment_both=rbind(cellage_promotes_enrichment,
                             cellage_inhibits_enrichment)

save_csv(data = cellage_enrichment_both,
         file_name = 'cellage_induce_inhibit_enrichment',
         path = save_dir_csv)
#####