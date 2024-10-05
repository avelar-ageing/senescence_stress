#DEG analysis
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

human_pc_entrez=getBM(attributes=c('external_gene_name','ensembl_gene_id',
                                   'entrezgene_id'),
                      filters = 'biotype',
                      values = c('protein_coding'),
                      mart = ensembl100)
human_pc_entrez=human_pc_entrez[!is.na(human_pc_entrez$entrezgene_id),]
#####
arrest_degs_merged=read.csv('/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/arrest_degs_final.csv')
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
#####
#Marian pathways
marian_csv=read.csv('/Users/ravelarvargas/Downloads/marian/pathways_check.csv')

library(rWikiPathways)
arrest_degs_merged_plot=arrest_degs_merged
marian_wiki=marian_csv[marian_csv$type=='wiki',]

pathway_get=
  lapply(1:nrow(marian_wiki), function(i) c(marian_wiki$code[i], marian_wiki$pathway[i]))

genes_plot=c()

pathway_get=c(pathway_get,
              list(c('WP111','Oxidative phosphorylation')))
pathway_genes=do.call('rbind',lapply(pathway_get,function(get_pathway){
  message(get_pathway)
  pathway_accession=get_pathway[1]
  
  gene_list_pathway=getXrefList(pathway_accession, 'L')
  # getXrefList(pathway_accession, 'En')
  # getPathwayInfo(pathway_accession)
  # findPathwaysByXref('TCA','H')
  # intersect(test_1,test_2)
  if(sum(grepl(gene_list_pathway,pattern='ENSG'))>0){
    print('TRUE')
    temp_grepl=gene_list_pathway[grepl(gene_list_pathway,pattern='ENSG')]
    temp_add=human_pc_entrez[human_pc_entrez$ensembl_gene_id%in%temp_grepl,][['entrezgene_id']]
    temp_remove=gene_list_pathway[!grepl(gene_list_pathway,pattern='ENSG')]
    gene_list_pathway=c(temp_remove,temp_add)
  }
  message(length(gene_list_pathway))
  human_pc_entrez_pathway=data.frame(gene=human_pc_entrez[human_pc_entrez$entrezgene_id%in%gene_list_pathway,][['external_gene_name']])
  message(nrow(human_pc_entrez_pathway))
  human_pc_entrez_pathway$code=get_pathway[1]
  human_pc_entrez_pathway$pathway=get_pathway[2]
  
  if(is.null(genes_plot)){
    # arrest_degs_merged_plot<<-extract_label_genes(df = arrest_degs_merged,
    #                                          new_column = 'group',
    #                                          gene_vector = c(human_pc_entrez_pathway$gene),
    #                                          gene_column = 'gene',
    #                                          gene_label = get_pathway[2])
    # time_degs_plot<<-extract_label_genes(df = time_degs,
    #                                      new_column = 'group',
    #                                      gene_vector = c(human_pc_entrez_pathway$gene),
    #                                      gene_column = 'ensembl',
    #                                      gene_label = get_pathway[2])
  }else{
    # arrest_degs_merged_plot<<-extract_label_genes(df = arrest_degs_merged,
    #                                          new_column = 'group',
    #                                          gene_vector = c(human_pc_entrez_pathway$gene),
    #                                          gene_column = 'gene',
    #                                          gene_label = get_pathway[2],
    #                                          new_df = arrest_degs_merged_plot)
    # time_degs_plot<<-extract_label_genes(df = time_degs,
    #                                      new_column = 'group',
    #                                      gene_vector = c(human_pc_entrez_pathway$gene),
    #                                      gene_column = 'ensembl',
    #                                      gene_label = get_pathway[2],
    #                                      new_df=time_degs_plot)
  }
  
  genes_plot<<-c(genes_plot,human_pc_entrez_pathway$gene)
  return(human_pc_entrez_pathway)
}))
pathway_genes=pathway_genes%>%dplyr::arrange(code,gene)
pathway_genes_unique=unique(pathway_genes)
# save_csv(pathway_genes,'wikipathway_genes',path = '/Users/ravelarvargas/Downloads/marian')

#KEGG
library(KEGGREST)
marian_kegg=marian_csv[marian_csv$type=='kegg',]
all_kegg=marian_kegg$code

kegg_modules=download_KEGG('hsa', keggType = "MKEGG", keyType = "kegg")$KEGGPATHID2EXTID
kegg_modules=kegg_modules[kegg_modules$from%in%all_kegg,]
kegg_modules=merge(kegg_modules,human_pc_entrez,by.x='to',by.y='entrezgene_id')
kegg_modules=unique(kegg_modules[,2:(ncol(kegg_modules)-1)])
colnames(kegg_modules)=c('code','gene')
kegg_modules$pathway='Pentose Phosphate Pathway + PRPP biosynthesis'
kegg_modules=kegg_modules%>%dplyr::select(gene,code,pathway)
kegg_modules=kegg_modules%>%dplyr::arrange(code,gene)
# save_csv(kegg_modules,'mkegg_genes',path = '/Users/ravelarvargas/Downloads/marian')

#biocyc
##purine
purine_temp=data.frame(fread('/Users/ravelarvargas/Downloads/marian/pathway-genes-PWY-841.txt'))
colnames(purine_temp)[1]='gene'
purine_temp$code='PWY-841'
purine_temp$pathway='Superpathway of purine nucleotides de novo biosynthesis I'
purine_temp=unique(purine_temp[purine_temp$Organism=='Homo sapiens',]%>%dplyr::select(Gene.name,code,pathway))
purine_temp$type='metabolic enzyme'
purine_temp$keep=TRUE

##pyrimidine
pyrimidine_temp=data.frame(fread('/Users/ravelarvargas/Downloads/marian/pathway-genes-PWY-7211.txt'))
colnames(pyrimidine_temp)[1]='gene'
pyrimidine_temp$code='PWY-7211'
pyrimidine_temp$pathway='Superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis'
pyrimidine_temp=unique(pyrimidine_temp[pyrimidine_temp$Organism=='Homo sapiens',]%>%dplyr::select(Gene.name,code,pathway))
pyrimidine_temp$type='metabolic enzyme'
pyrimidine_temp$keep=TRUE

#TCA
TCA_temp=data.frame(fread('/Users/ravelarvargas/Downloads/marian/pathway-genes-PWY66-398.txt'))
colnames(TCA_temp)[1]='gene'
TCA_temp$code='PWY66-398'
TCA_temp$pathway='TCA cycle III (animals)'
TCA_temp=unique(TCA_temp[TCA_temp$Organism=='Homo sapiens',]%>%dplyr::select(Gene.name,code,pathway))
TCA_temp$type='metabolic enzyme'
TCA_temp$keep=TRUE

#glycolysis
glycolysis_temp=data.frame(fread('/Users/ravelarvargas/Downloads/marian/pathway-genes-ANAGLYCOLYSIS-PWY.txt'))
colnames(glycolysis_temp)[1]='gene'
glycolysis_temp$code='ANAGLYCOLYSIS-PWY'
glycolysis_temp$pathway='glycolysis III (from glucose)'
glycolysis_temp=unique(glycolysis_temp[glycolysis_temp$Organism=='Homo sapiens',]%>%dplyr::select(Gene.name,code,pathway))
glycolysis_temp$type='metabolic enzyme'
glycolysis_temp$keep=TRUE

biocyc_pathways=rbind(purine_temp,pyrimidine_temp, TCA_temp,glycolysis_temp)
# save_csv(biocyc_pathways,'biocyc_genes',path = '/Users/ravelarvargas/Downloads/marian')
#####

glut=c('GSS', 'GSR', 'GCLC', 'GCLM', 'GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX5',
       'GPX6', 'GPX7', 'GPX8', 'GSTA1', 'GSTM1', 'GSTP1', 'GSTT1', 'GSTO1')
glut=glut[!glut%in%c('GCLM', 'GPX2', 'GPX4', 'GPX7', 'GSR', 'GSS', 'GSTM1',
                     'GSTP1')]

sphingo=c('SPTLC1', 'SPTLC2', 'SPTLC3', 'SGPL1', 'SGPP1', 'SGPP2', 'SMPD1')
sphingo=sphingo[!sphingo%in%c('SGPL1', 'SPTLC1', 'SPTLC2', 'SPTLC3')]

glyco=c('HK1', 'HK2', 'HK3', 'PFKM', 'PFKP', 'PFKL', 'GAPDH', 'PGK1',
        'PGK2', 'PKM', 'PKLR')
glyco=glyco[!glyco%in%c('PFKM', 'PFKL', 'HK1')]

pentose=c('G6PD', 'PGD', 'RPE', 'RPIA')
rm(pentose)

tca=c('CS', 'IDH1', 'IDH2', 'IDH3A', 'IDH3B', 'IDH3G', 'OGDH',
      'SDHA', 'SDHB', 'SDHC', 'SDHD', 'FH')
tca=tca[!tca%in%c('IDH3A','IDH3B','IDH3G','SDHA','SDHC','SDHD','CS','SDHB')]

nad=c('NAMPT', 'NMNAT1', 'NMNAT2', 'NMNAT3', 'SIRT1', 'SIRT2',
      'SIRT3', 'SIRT4', 'SIRT5', 'SIRT6', 'SIRT7', 'PARP1', 'PARP2')
nad=nad[!nad%in%c('SIRT5', 'SIRT6', 'SIRT7', 'SIRT3', 'SIRT1', 'NMNAT1')]

kynu=c('IDO1', 'IDO2', 'TDO2', 'KAT1', 'KAT2', 'KAT3', 'KAT4', 'KYNU')

leuk=c('ALOX5', 'ALOX5AP', 'LTA4H', 'LTC4S','PLA2G4A','CYSLTR1')

prostag=c('PTGS1', 'PTGS2', 'PTGDS', 'PTGES', 'PTGFR')

#leuk
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(leuk),
                               gene_column = 'gene',
                               gene_label = 'Leukotrienes')

#prostag
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(prostag),
                               gene_column = 'gene',
                               gene_label = 'Prostaglandins')

#glut
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(glut),
                               gene_column = 'gene',
                               gene_label = 'Glutathione')

#sphingo
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(sphingo),
                               gene_column = 'gene',
                               gene_label = 'Sphingolipids')

#glyco
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(glyco),
                               gene_column = 'gene',
                               gene_label = 'Glycolysis')
#pentose
# arrest_degs_merged=label_genes(arrest_degs_merged,
#                                new_column = 'group',
#                                gene_vector = c(pentose),
#                                gene_column = 'gene',
#                                gene_label = 'Pentose Phosphate')
#tca
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(tca),
                               gene_column = 'gene',
                               gene_label = 'TCA')
#nad
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(nad),
                               gene_column = 'gene',
                               gene_label = 'NAD')
#kynu
arrest_degs_merged=label_genes(arrest_degs_merged,
                               new_column = 'group',
                               gene_vector = c(kynu),
                               gene_column = 'gene',
                               gene_label = 'Kynurenine')

arrest_heterogeneity_metabolism=compare_expression(degs=arrest_degs_merged,
                                         gene_col = 'gene',
                                         group = 'group_1',
                                         genes=c(glut,
                                                 sphingo,
                                                 glyco,
                                                 # pentose,
                                                 tca,
                                                 nad,
                                                 kynu,
                                                 leuk,
                                                 prostag),
                                         facet = 'group',
                                         # facet2 = 'arrest',
                                         flip_axis = TRUE,
                                         scale_temp = 2,
                                         fill_cap = 6,
                                         tilt_x = 90
)
save_p(arrest_heterogeneity_metabolism,
       file_name = 'heterogeneity_metabolism',
       save_dir = '/Volumes/GoogleDrive/My Drive/Grants/2023/Impetus_metabolism',
       p_width = 7.5,p_height=4.5)
#####
#temporal

time_degs=compile_time_files(dir_use='/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/ERP021140/',
                             deg_file='degs.csv')
time_degs$group_1[time_degs$group_1=='4_days']='4 days'
time_degs$group_1[time_degs$group_1=='10_days']='10 days'
time_degs$group_1[time_degs$group_1=='20_days']='20 days'
time_degs$cell_type[time_degs$cell_type=='fibroblast']='Fibroblast'
time_degs$cell_type[time_degs$cell_type=='keratinocyte']='Keratinocyte'
time_degs$cell_type[time_degs$cell_type=='melanocyte']='Melanocyte'

tca=tca[!tca%in%c('OGDH')]
prostag=prostag[!prostag%in%c('PTGFR')]

#glut
time_degs=label_genes(time_degs,
                               new_column = 'group',
                               gene_vector = c(glut),
                               gene_column = 'ensembl',
                               gene_label = 'Glutathione')

#leuk
time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = c(leuk),
                      gene_column = 'ensembl',
                      gene_label = 'Leukotrienes')

#prostag
time_degs=label_genes(time_degs,
                      new_column = 'group',
                      gene_vector = c(prostag),
                      gene_column = 'ensembl',
                      gene_label = 'Prostaglandins')

#sphingo
time_degs=label_genes(time_degs,
                               new_column = 'group',
                               gene_vector = c(sphingo),
                               gene_column = 'ensembl',
                               gene_label = 'Sphingolipids')

#glyco
time_degs=label_genes(time_degs,
                               new_column = 'group',
                               gene_vector = c(glyco),
                               gene_column = 'ensembl',
                               gene_label = 'Glycolysis')
#pentose
# time_degs=label_genes(time_degs,
#                                new_column = 'group',
#                                gene_vector = c(pentose),
#                                gene_column = 'ensembl',
#                                gene_label = 'Pentose Phosphate')
#tca
time_degs=label_genes(time_degs,
                               new_column = 'group',
                               gene_vector = c(tca),
                               gene_column = 'ensembl',
                               gene_label = 'TCA')
#nad
time_degs=label_genes(time_degs,
                               new_column = 'group',
                               gene_vector = c(nad),
                               gene_column = 'ensembl',
                               gene_label = 'NAD')
#kynu
time_degs=label_genes(time_degs,
                               new_column = 'group',
                               gene_vector = c(kynu),
                               gene_column = 'ensembl',
                               gene_label = 'Kynurenine')

compare_expression=function(degs,
         genes,
         gene_col='ensembl',
         group,
         facet=NULL,
         facet2=NULL,
         fill_col='log2FoldChange',
         pval_col='padj',
         scale_temp=1,
         log2_cutoff=log2(1.5),
         flip_axis=FALSE,
         fill_cap=NULL,
         tilt_x=FALSE,
         hide_ns=TRUE,
         text_size=10, 
         facet_text_size=10,
         legend_bottom=FALSE,
         gene_lab='Gene',
         condition_lab='Condition',
         graph_limits=NULL){
  
  degs_relevant=degs[degs[[gene_col]]%in%genes,]
  
  # Cap the fill values
  if (!is.null(fill_cap)) {
    degs_relevant[[fill_col]] = pmin(degs_relevant[[fill_col]], fill_cap)
  }
  
  degs_relevant[['x']]=degs_relevant[[gene_col]]
  degs_relevant[['y']]=degs_relevant[[group]]
  degs_relevant[['fill']]=degs_relevant[[fill_col]]
  degs_relevant=degs_relevant%>%rstatix::add_significance(pval_col)
  
  if(!is.null(facet)){
    degs_relevant[['facet_temp']]=degs_relevant[[facet]]
  }
  if(!is.null(facet2)){
    degs_relevant[['facet_temp2']]=degs_relevant[[facet2]]
  }
  
  graph_max=max(ceiling(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  graph_min=min(floor(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  if(graph_min>0){
    graph_min=0
  }
  
  if(is.null(graph_limits)){
    li <- c(graph_min, graph_max)
    la <- c(seq(graph_min,graph_max,scale_temp))
    br <- c(seq(graph_min,graph_max,scale_temp))
  }else{
    li <- c(graph_limits[1], graph_limits[2])
    la <- c(seq(graph_limits[1],graph_limits[2],scale_temp))
    br <- c(seq(graph_limits[1],graph_limits[2],scale_temp))
  }
  
  
  # Modify labels if fill_cap is not NULL
  if (!is.null(fill_cap)) {
    la[la == fill_cap] = paste0(">", fill_cap)
  }
  
  if(!is.null(log2_cutoff)){
    degs_relevant$padj.signif[abs(degs_relevant[[fill_col]])<log2_cutoff]='ns'
  }
  
  if(hide_ns){
    degs_relevant$padj.signif[degs_relevant$padj.signif=='ns']=''
    # degs_relevant$padj.signif=paste0('\n',degs_relevant$padj.signif)
    nudge_y=-0.2
  }else{
    nudge_y=0
  }
  
  if(!flip_axis){
    p_temp=degs_relevant%>%ggplot(aes(x=x,y=y,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), 
                           low = "blue", 
                           mid = "white", 
                           high = "red2", 
                           midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',
                                                  title.position = 'top',
                                                  title.hjust=0.5,
                                                  barwidth = unit(2, "inches")), # increase this for a longer legend
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(gene_lab)+
      ylab(condition_lab)
  }else{
    p_temp=degs_relevant%>%ggplot(aes(x=y,y=x,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), 
                           low = "blue", 
                           mid = "white", 
                           high = "#ff0000", 
                           midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',
                                                  title.position = 'top',
                                                  title.hjust=0.5,
                                                  barwidth = unit(2, "inches")), # increase this for a longer legend
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(condition_lab)+
      ylab(gene_lab)
  }
  if (tilt_x) {
    p_temp = p_temp + theme(axis.text.x=element_text(angle=45,hjust=1, size=text_size))
  }
  
  # add size control for axis labels, axis titles, legend and geom_text
  p_temp = p_temp + theme(
    text = element_text(size=text_size),
    axis.title = element_text(size=text_size),
    axis.text = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    strip.text = element_text(size=facet_text_size),
    # modify geom_text size
    plot.subtitle = element_text(size=text_size),
    plot.caption = element_text(size=text_size)
  )
  
  # move the legend to the bottom if legend_bottom is TRUE
  if (legend_bottom) {
    p_temp = p_temp + theme(legend.position = "bottom")
  }
  
  if(!is.null(facet)&is.null(facet2)){
    if(!flip_axis){
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_x')
    }else{
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_y')
    }
  }else if(!is.null(facet)&!is.null(facet2)){
    p_temp=p_temp+
      facet_grid(facet_temp~facet_temp2,scales = 'free',space='free')
  }
  return(p_temp)
}

library(patchwork)
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

temp_legend=extract_legend(compare_expression(degs=time_degs,
                                              gene_col = 'ensembl',
                                              group = 'group_1',
                                              genes=c(glut,
                                                      sphingo,
                                                      glyco,
                                                      leuk,
                                                      # pentose,
                                                      prostag,
                                                      nad,
                                                      tca,
                                                      kynu),
                                              facet = 'cell_type',
                                              facet2 = 'group',
                                              flip_axis = FALSE,
                                              scale_temp = 2,fill_cap = 6,
                                              tilt_x = TRUE,
                                              text_size = 13,
                                              legend_bottom = TRUE,
                                              condition_lab='Days Post-Irradiation',
                                              graph_limits = c(-2,4)))

time_degs$group_1=factor(time_degs$group_1,
                                 levels=rev(c('4 days','10 days','20 days')))
p_1=compare_expression(degs=time_degs,
                                         gene_col = 'ensembl',
                                         group = 'group_1',
                                       genes=c(glut,
                                               sphingo,
                                               glyco,
                                               # pentose,
                                               prostag),
                                         facet = 'cell_type',
                                         facet2 = 'group',
                                         flip_axis = FALSE,
                                         scale_temp = 9,fill_cap = 6,
                                         tilt_x = TRUE,
                                         text_size = 13,
                                         legend_bottom = TRUE,
                                         condition_lab='Days Post-Irradiation',
                                       graph_limits = c(-9,9))

save_p(p_1,file_name = 'get_legend',
       save_dir = '/Users/ravelarvargas/Downloads')
  # theme(legend.position = "none")
p_2=compare_expression(degs=time_degs,
                     gene_col = 'ensembl',
                     group = 'group_1',
                     genes=c(tca,
                             nad,
                             kynu,
                             leuk),
                     facet = 'cell_type',
                     facet2 = 'group',
                     flip_axis = FALSE,
                     scale_temp = 2,fill_cap = 6,
                     tilt_x = TRUE,
                     text_size = 13,
                     legend_bottom = TRUE,
                     condition_lab='Days Post-Irradiation',
                     graph_limits = c(-4,4)
  )+
  theme(legend.position = "none")

library("gridExtra")
temporal_metabolism=grid.arrange(arrangeGrob(p_1, p_2, ncol = 1),
             temp_legend, nrow = 2, heights = c(10, 1))
save_p(temporal_metabolism,
       file_name = 'temporal_metabolism',
       save_dir = '/Volumes/GoogleDrive/My Drive/Grants/2023/Impetus_metabolism',
       p_width = 10,p_height = 9)

temporal_metabolism_full=compare_expression(degs=time_degs,
                   gene_col = 'ensembl',
                   group = 'group_1',
                   genes=c(glut,
                           sphingo,
                           glyco,
                           # pentose,
                           prostag,
                           tca,
                           nad,
                           # kynu,
                           leuk),
                   facet = 'cell_type',
                   facet2 = 'group',
                   flip_axis = FALSE,
                   scale_temp = 2,fill_cap = 6,
                   tilt_x = 45,
                   text_size = 13,
                   legend_bottom = TRUE,
                   condition_lab='Days Post-Irradiation')

save_p(temporal_metabolism_full,
       file_name = 'temporal_metabolism_full',
       save_dir = '/Volumes/GoogleDrive/My Drive/Grants/2023/Impetus_metabolism',
       p_width = 14,p_height = 5)

#####
#metabolism pathway overlaps
##read in metabolism genes
metabolism_background=unique(read.table('/Users/ravelarvargas/Downloads/genes.tsv',header = TRUE)[['geneSymbols']])

metabolism_genes=read.csv('/Users/ravelarvargas/Downloads/marian/pathways_genes.csv')
metabolism_background=unique(c(metabolism_background,metabolism_genes$gene))
metabolism_genes=metabolism_genes[metabolism_genes$keep,]
metabolism_genes=metabolism_genes[metabolism_genes$gene%in%arrest_degs_merged$gene,]
arrest_metabolism=arrest_degs_merged[arrest_degs_merged$gene%in%
                                       metabolism_background,]
metabolism_background=metabolism_background[metabolism_background%in%
                                              arrest_degs_merged$gene]

# table(arrest_metabolism%>%dplyr::filter(sig=='y')%>%dplyr::select(group_1_dir))
#get just oxphos
oxo_temp=metabolism_genes[metabolism_genes$class=='Oxidative phosphorylation',]
oxo_degs=arrest_degs_merged[arrest_degs_merged$gene%in%oxo_temp$gene,]
oxo_degs=merge(oxo_degs,human_pc,by.x='gene',by.y='external_gene_name')
oxo_degs=oxo_degs%>%dplyr::select(ensembl_gene_id,
                                           group_1,
                                           log2FoldChange
                                            )

all_group=unique(oxo_degs$group_1)
oxo_degs=oxo_degs%>%dplyr::arrange(ensembl_gene_id)
genes_return=unique(oxo_degs$ensembl_gene_id)
oxo_degs=do.call('cbind',lapply(all_group,function(per_group){
  oxo_degs_temp=oxo_degs[oxo_degs$group_1==per_group,]
  oxo_degs_temp=oxo_degs_temp[,colnames(oxo_degs_temp)!='group_1']
  colnames(oxo_degs_temp)[2]=paste0(per_group,'_FC')
  return(oxo_degs_temp%>%dplyr::select(paste0(per_group,'_FC')))
}))
oxo_degs=cbind(genes_return,oxo_degs)

colnames(oxo_degs)[1]='Identifier'
oxo_degs=oxo_degs%>%dplyr::select(Identifier,
                                  'Contact-inhibited CQ_FC',
                                  'Serum-starved CQ_FC',
                                  RS_FC,
                                  SIPS_FC,
                                  OIS_FC)
oxo_degs[['System Code']]='En'
write.table(oxo_degs,file = 
              '/Users/ravelarvargas/Downloads/marian/oxos.txt',
            quote = FALSE,sep='\t',row.names = FALSE)

#overlaps
save_csv(metabolism_background,
         file_name = 'background_pc_arrest_filtered',
         path = '/Users/ravelarvargas/Downloads/marian')

save_csv(metabolism_genes,
         file_name = 'genelist_pc_arrest_filtered',
         path = '/Users/ravelarvargas/Downloads/marian')

#overlap with self
meta_self_overlap=overlap_within_df(dataframe = metabolism_genes,
                  group_col = 'class',
                  gene_col = 'gene',
                  background = metabolism_background,
                  other_cols = NULL,remove_self = TRUE)

temp_table_x=rev(names(sort(table(meta_self_overlap$class_2_1))))
meta_self_overlap$class_2_1=
  factor(meta_self_overlap$class_2_1,levels=temp_table_x)
temp_table_y=rev(names(sort(table(meta_self_overlap$class_2_2))))
meta_self_overlap$class_2_2=
  factor(meta_self_overlap$class_2_2,levels=temp_table_y)

#grouping some pathways
metabolism_self_overlaps_p_2=plot_self_overlaps(meta_self_overlap,
                   x = 'class_2_1',
                   y = 'class_2_2',
                   facet_col = NULL,angle_x = TRUE)
save_p(metabolism_self_overlaps_p_2,
       file_name = 'metabolism_self_overlaps_synthase',
       save_dir='/Users/ravelarvargas/Downloads/marian')

#overlap with degs
arrest_degs_merged_sig=arrest_degs_merged[arrest_degs_merged$sig=='y',]

temp_filter=names(table(metabolism_genes$class_3)[table(metabolism_genes$class_3)>=15])
metabolism_genes_filtered=metabolism_genes[metabolism_genes$class_3%in%temp_filter,]
metabolism_deg_overlap=overlap_function(df_1 = arrest_degs_merged_sig,
                 gene_col_1 = 'gene',
                 group_col_1 = 'dir_accession',
                 carry_col_1 = c('group_1','direction_1'),
                 df_2 = metabolism_genes_filtered,
                 gene_col_2 = 'gene',
                 group_col_2 = 'class_3',
                 background = metabolism_background)

metabolism_deg_overlap$direction_1[metabolism_deg_overlap$direction_1=='up']=
  'Up in Arrest'
metabolism_deg_overlap$direction_1[metabolism_deg_overlap$direction_1=='down']=
  'Down in Arrest'
overlaps_metabolism=create_overlap_plot(metabolism_deg_overlap,
                                        x = 'group_1',
                                        y='direction_1',
                                        facet_col = 'class_3',
                                        x_tilt = 90,xlab = 'Arrest Class',
                                        ylab = 'Arrest Direction')

save_p(overlaps_metabolism,
       file_name = 'degs_vs_metabolism',
       save_dir = '/Users/ravelarvargas/Downloads/marian/',
       p_width = 11,p_height=5.3)

##with synthase
metabolism_deg_overlap_grouped=overlap_function(df_1 = arrest_degs_merged_sig,
                                        gene_col_1 = 'gene',
                                        group_col_1 = 'dir_accession',
                                        carry_col_1 = c('group_1','direction_1'),
                                        df_2 = metabolism_genes,
                                        gene_col_2 = 'gene',
                                        group_col_2 = 'class_2',
                                        background = metabolism_background)

metabolism_deg_overlap_grouped$group_1=factor(metabolism_deg_overlap_grouped$group_1,
                                      levels=c('Contact-inhibited CQ','Serum-starved CQ',
                                               'RS','OIS','SIPS'))
metabolism_deg_overlap_grouped$direction_1[metabolism_deg_overlap_grouped$direction_1=='up']=
  'Up in Arrest'
metabolism_deg_overlap_grouped$direction_1[metabolism_deg_overlap_grouped$direction_1=='down']=
  'Down in Arrest'
overlaps_metabolism_2=create_overlap_plot(metabolism_deg_overlap_grouped,
                                        x = 'group_1',
                                        y='direction_1',
                                        facet_col = 'class_2',
                                        x_tilt = 90,xlab = 'Arrest Class',
                                        ylab = 'Arrest Direction')

save_p(overlaps_metabolism_2,
       file_name = 'degs_vs_metabolism_synthase',
       save_dir = '/Users/ravelarvargas/Downloads/marian/',
       p_width = 11,p_height=5.3)

#Function to make heatmaps from subset of genes
degs=arrest_degs_merged

metabolism_heatmap=build_heatmap(degs = arrest_degs_merged[arrest_degs_merged$gene%in%
                                               metabolism_genes$gene,],
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
              normalisation = 'vst_full',top_x=NULL)

save_pheatmap(pheatmap_to_save=metabolism_heatmap,
              file_name='metabolism_all_samples',
              save_dir='/Users/ravelarvargas/Downloads/marian/',
              p_width=3000,
              p_height=3000)

genelist_heatmap=function(degs,
                        gene_col='gene',
                        group_col='group_1',
                        plot_col='log2FoldChange',
                        scale='row'
                        # gene_list,
                        # gene_list_gene_col='gene'
                        ){
  # degs_temp=degs[degs[[gene_col]]%in%gene_list[[gene_list_gene_col]],]
  degs_temp=degs
  degs_temp_pheat=degs_temp[,colnames(degs_temp)%in%c(gene_col,
                                                      group_col,
                                                      plot_col)]
  all_condition=unique(degs_temp_pheat[[group_col]])
  to_pheat=do.call('cbind',lapply(all_condition,function(per_group){
    degs_temp_pheat_get=degs_temp_pheat[degs_temp_pheat[[group_col]]==per_group,]
    degs_temp_pheat_get=degs_temp_pheat_get%>%dplyr::arrange_(gene_col)
    rownames(degs_temp_pheat_get)=NULL
    rows_temp=degs_temp_pheat_get[[gene_col]]
    get_genes=degs_temp_pheat_get[,!colnames(degs_temp_pheat_get)%in%
                                    c(group_col,gene_col)]
    get_genes=matrix(get_genes,ncol=1)
    rownames(get_genes)=rows_temp
    colnames(get_genes)=per_group
    return(get_genes)
  }))
  
  pheatmap(to_pheat, cluster_rows=TRUE, show_rownames=FALSE,
           show_colnames=TRUE,cluster_cols=TRUE,
           scale=scale,
           border_color =NA,fontsize = 12,color = greenred(75))
}

deg_heatmap=genelist_heatmap(degs = arrest_degs_merged[arrest_degs_merged$gene%in%
                                             metabolism_genes$gene,])

save_pheatmap(pheatmap_to_save=deg_heatmap,
              file_name='metabolism_deg_heatmap',
              save_dir='/Users/ravelarvargas/Downloads/marian/',
              p_width=3000,
              p_height=3000)

# pathway_dir='/Users/ravelarvargas/Downloads/wikipathways/'

# temp_files=gsub(list.files(pathway_dir),pattern='.csv',replacement='')

pathway_get=list(c('WP5171','Leukotriene'),
                 c('WP534','Glycolysis'),
                 c('WP100','Glutathione'),
                 c('WP465','Tryptophan'),
                 c('WP3630','NAD_ageing'),
                 c('WP98','Prostaglandin'),
                 c('WP78','TCA'),
                 c('WP5424','Omega3'))

genes_plot=c()
arrest_degs_merged_plot=arrest_degs_merged

pathway_genes=do.call('rbind',lapply(pathway_get,function(get_pathway){
  pathway_accession=get_pathway[1]
  
  # gene_list_pathway=getXrefList(pathway_accession, 'L')
  getXrefList(pathway_accession, 'En')
  getPathwayInfo(pathway_accession)
  findPathwaysByXref('TCA','H')
  intersect(test_1,test_2)
  if(sum(grepl(gene_list_pathway,pattern='ENSG'))>0){
    temp_grepl=gene_list_pathway[grepl(gene_list_pathway,pattern='ENSG')]
    temp_add=human_pc_entrez[human_pc_entrez$ensembl_gene_id%in%temp_grepl,][['entrezgene_id']]
    temp_remove=gene_list_pathway[!grepl(gene_list_pathway,pattern='ENSG')]
    gene_list_pathway=c(temp_remove,temp_add)
  }
  human_pc_entrez_pathway=data.frame(gene=human_pc_entrez[human_pc_entrez$entrezgene_id%in%gene_list_pathway,][['external_gene_name']])
  human_pc_entrez_pathway$code=get_pathway[1]
  human_pc_entrez_pathway$pathway=get_pathway[2]
  
  if(is.null(genes_plot)){
    # arrest_degs_merged_plot<<-extract_label_genes(df = arrest_degs_merged,
    #                                          new_column = 'group',
    #                                          gene_vector = c(human_pc_entrez_pathway$gene),
    #                                          gene_column = 'gene',
    #                                          gene_label = get_pathway[2])
    # time_degs_plot<<-extract_label_genes(df = time_degs,
    #                                      new_column = 'group',
    #                                      gene_vector = c(human_pc_entrez_pathway$gene),
    #                                      gene_column = 'ensembl',
    #                                      gene_label = get_pathway[2])
  }else{
    # arrest_degs_merged_plot<<-extract_label_genes(df = arrest_degs_merged,
    #                                          new_column = 'group',
    #                                          gene_vector = c(human_pc_entrez_pathway$gene),
    #                                          gene_column = 'gene',
    #                                          gene_label = get_pathway[2],
    #                                          new_df = arrest_degs_merged_plot)
    # time_degs_plot<<-extract_label_genes(df = time_degs,
    #                                      new_column = 'group',
    #                                      gene_vector = c(human_pc_entrez_pathway$gene),
    #                                      gene_column = 'ensembl',
    #                                      gene_label = get_pathway[2],
    #                                      new_df=time_degs_plot)
  }
  
  genes_plot<<-c(genes_plot,human_pc_entrez_pathway$gene)
  return(human_pc_entrez_pathway)
}))

compare_expression=function(degs,
                            genes,
                            gene_col='ensembl',
                            group,
                            facet=NULL,
                            facet2=NULL,
                            fill_col='log2FoldChange',
                            pval_col='padj',
                            scale_temp=1,
                            log2_cutoff=log2(1.5),
                            flip_axis=FALSE,
                            fill_cap=NULL,
                            tilt_x=90,
                            hide_ns_lab=TRUE,
                            text_size=10, 
                            facet_text_size=10,
                            legend_bottom=FALSE,
                            gene_lab='Gene',
                            condition_lab='Condition',
                            hide_ns=TRUE,
                            facet_n=3){
  
  degs_relevant=degs[degs[[gene_col]]%in%genes,]
  
  if(hide_ns){
    all_gene=unique(degs_relevant[[gene_col]])
    genes_filter=do.call(c,lapply(all_gene,function(get_gene){
      degs_relevant_use=degs_relevant[degs_relevant[[gene_col]]==get_gene,]
      degs_relevant_use=degs_relevant_use[degs_relevant_use$sig=='y',]
      if(nrow(degs_relevant_use)>0){
        return(get_gene)
      }else{
        return(NULL)
      }
    }))
    degs_relevant=degs_relevant[degs_relevant[[gene_col]]%in%genes_filter,]
  }
  
  # Cap the fill values
  if (!is.null(fill_cap)) {
    degs_relevant[[fill_col]] = pmin(degs_relevant[[fill_col]], fill_cap)
  }
  
  degs_relevant[['x']]=degs_relevant[[gene_col]]
  degs_relevant[['y']]=degs_relevant[[group]]
  degs_relevant[['fill']]=degs_relevant[[fill_col]]
  degs_relevant=degs_relevant%>%rstatix::add_significance(pval_col)
  
  if(!is.null(facet)){
    degs_relevant[['facet_temp']]=degs_relevant[[facet]]
  }
  if(!is.null(facet2)){
    degs_relevant[['facet_temp2']]=degs_relevant[[facet2]]
  }
  
  graph_max=max(ceiling(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  graph_min=min(floor(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  if(graph_min>0){
    graph_min=0
  }
  
  li <- c(graph_min, graph_max)
  la <- c(seq(graph_min,graph_max,scale_temp))
  br <- c(seq(graph_min,graph_max,scale_temp))
  
  # Modify labels if fill_cap is not NULL
  if (!is.null(fill_cap)) {
    la[la == fill_cap] = paste0(">", fill_cap)
  }
  
  if(!is.null(log2_cutoff)){
    degs_relevant$padj.signif[abs(degs_relevant[[fill_col]])<log2_cutoff]='ns'
  }
  
  if(hide_ns_lab){
    degs_relevant$padj.signif[degs_relevant$padj.signif=='ns']=''
    # degs_relevant$padj.signif=paste0('\n',degs_relevant$padj.signif)
    nudge_y=-0.2
  }else{
    nudge_y=0
  }
  
  if(!flip_axis){
    p_temp=degs_relevant%>%ggplot(aes(x=x,y=y,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',title.position = 'top',title.hjust=0.5),
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(gene_lab)+
      ylab(condition_lab)
  }else{
    p_temp=degs_relevant%>%ggplot(aes(x=y,y=x,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), low = "blue", mid = "white", high = "red2", midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',title.position = 'top',title.hjust=0.5),
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(condition_lab)+
      ylab(gene_lab)
  }
  if (tilt_x==45) {
    p_temp = p_temp + theme(axis.text.x=element_text(angle=45,hjust=1, size=text_size))
  }else if(tilt_x==90){
    p_temp = p_temp + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=text_size))
    
  }else if(tilt_x==0){
    p_temp + theme(axis.text.x=element_text(size=text_size))
  }
  
  # add size control for axis labels, axis titles, legend and geom_text
  p_temp = p_temp + theme(
    text = element_text(size=text_size),
    axis.title = element_text(size=text_size),
    axis.text = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    strip.text = element_text(size=facet_text_size),
    # modify geom_text size
    plot.subtitle = element_text(size=text_size),
    plot.caption = element_text(size=text_size)
  )
  
  # move the legend to the bottom if legend_bottom is TRUE
  if (legend_bottom) {
    p_temp = p_temp + theme(legend.position = "bottom")
  }
  
  if(!is.null(facet)&is.null(facet2)){
    if(!flip_axis){
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_x',ncol = facet_n)
    }else{
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_y',ncol=facet_n)
    }
  }else if(!is.null(facet)&!is.null(facet2)){
    p_temp=p_temp+
      facet_grid(facet_temp~facet_temp2,scales = 'free',space='free',cols = facet_n)
  }
  return(p_temp)
}

arrest_degs_merged_plot$group_1=factor(arrest_degs_merged_plot$group_1,
                                          levels=c('Contact-inhibited CQ',
                                                   'Serum-starved CQ',
                                                   'RS','SIPS','OIS'))
arrest_heterogeneity_metabolism=compare_expression(degs=arrest_degs_merged_plot,
                                                   gene_col = 'gene',
                                                   group = 'group_1',
                                                   genes=unique(genes_plot),
                                                   facet = 'group',
                                                   # facet2 = 'arrest',
                                                   flip_axis = TRUE,
                                                   scale_temp = 2,
                                                   fill_cap = 6,
                                                   tilt_x = 90,facet_n = 4
)
save_p(arrest_heterogeneity_metabolism,
       file_name = 'metabolism_arrest_genes',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_figures/',
       p_width = 10,p_height=10)

#Overlap
arrest_degs_merged$group_1=as.character(arrest_degs_merged$group_1)
arrest_degs_merged_sig=arrest_degs_merged[arrest_degs_merged$sig=='y',]
metabolism_arrest_overlaps=overlap_function(df_1 = arrest_degs_merged_sig,
                 df_2 = pathway_genes,
                 gene_col_1 = 'gene',
                 gene_col_2 = 'gene',
                 group_col_1 = 'dir_accession',
                 carry_col_1 = c('group_1','direction_1'),
                 group_col_2 = 'pathway',
                 background = unique(arrest_degs_merged$gene))

save_csv(metabolism_arrest_overlaps,file_name = 'metabolism_arrest_fishers',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/')

metabolism_arrest_overlaps$group_1[metabolism_arrest_overlaps$group_1=='Contact-inhibited CQ']='Contact-inhibited\nCQ'
metabolism_arrest_overlaps$group_1[metabolism_arrest_overlaps$group_1=='Serum-starved CQ']='Serum-starved\nCQ'
metabolism_arrest_overlaps$group_1=factor(metabolism_arrest_overlaps$group_1,
                                          levels=c('Contact-inhibited\nCQ',
                                                   'Serum-starved\nCQ',
                                                   'RS','SIPS','OIS'))
metabolism_arrest_overlaps$direction_1[metabolism_arrest_overlaps$direction_1=='up']='Up\nin Arrest'
metabolism_arrest_overlaps$direction_1[metabolism_arrest_overlaps$direction_1=='down']='Down\nin Arrest'
overlaps_metabolism=create_overlap_plot(metabolism_arrest_overlaps,
                    x = 'group_1',
                    y='direction_1',
                    facet_col = 'pathway',x_tilt = 90)

save_p(overlaps_metabolism,
       file_name = 'metabolism_arrest',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_figures/')

#####
#Temporal

#compare expression
compare_expression_temporal=function(degs,
                            genes,
                            gene_col='ensembl',
                            group,
                            facet=NULL,
                            facet2=NULL,
                            fill_col='log2FoldChange',
                            pval_col='padj',
                            scale_temp=1,
                            log2_cutoff=log2(1.5),
                            flip_axis=FALSE,
                            fill_cap=NULL,
                            tilt_x=0,
                            hide_ns=TRUE,
                            text_size=10, 
                            facet_text_size=10,
                            legend_bottom=FALSE,
                            gene_lab='Gene',
                            condition_lab='Condition',
                            graph_limits=NULL){
  
  degs_relevant=degs[degs[[gene_col]]%in%genes,]
  
  # Cap the fill values
  if (!is.null(fill_cap)) {
    degs_relevant[[fill_col]] = pmin(degs_relevant[[fill_col]], fill_cap)
  }
  
  degs_relevant[['x']]=degs_relevant[[gene_col]]
  degs_relevant[['y']]=degs_relevant[[group]]
  degs_relevant[['fill']]=degs_relevant[[fill_col]]
  degs_relevant=degs_relevant%>%rstatix::add_significance(pval_col)
  
  if(!is.null(facet)){
    degs_relevant[['facet_temp']]=degs_relevant[[facet]]
  }
  if(!is.null(facet2)){
    degs_relevant[['facet_temp2']]=degs_relevant[[facet2]]
  }
  
  graph_max=max(ceiling(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  graph_min=min(floor(degs_relevant[[fill_col]]/scale_temp)[!is.infinite(ceiling(degs_relevant[[fill_col]]))])*scale_temp
  if(graph_min>0){
    graph_min=0
  }
  
  if(is.null(graph_limits)){
    li <- c(graph_min, graph_max)
    la <- c(seq(graph_min,graph_max,scale_temp))
    br <- c(seq(graph_min,graph_max,scale_temp))
  }else{
    li <- c(graph_limits[1], graph_limits[2])
    la <- c(seq(graph_limits[1],graph_limits[2],scale_temp))
    br <- c(seq(graph_limits[1],graph_limits[2],scale_temp))
  }
  
  
  # Modify labels if fill_cap is not NULL
  if (!is.null(fill_cap)) {
    la[la == fill_cap] = paste0(">", fill_cap)
  }
  
  if(!is.null(log2_cutoff)){
    degs_relevant$padj.signif[abs(degs_relevant[[fill_col]])<log2_cutoff]='ns'
  }
  
  if(hide_ns){
    degs_relevant$padj.signif[degs_relevant$padj.signif=='ns']=''
    # degs_relevant$padj.signif=paste0('\n',degs_relevant$padj.signif)
    nudge_y=-0.2
  }else{
    nudge_y=0
  }
  
  if(!flip_axis){
    p_temp=degs_relevant%>%ggplot(aes(x=x,y=y,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), 
                           low = "blue", 
                           mid = "white", 
                           high = "red2", 
                           midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',
                                                  title.position = 'top',
                                                  title.hjust=0.5,
                                                  barwidth = unit(2, "inches")), # increase this for a longer legend
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(gene_lab)+
      ylab(condition_lab)
  }else{
    p_temp=degs_relevant%>%ggplot(aes(x=y,y=x,fill=fill))+
      geom_tile(colour='black')+
      geom_text(aes(label=padj.signif),nudge_y=nudge_y)+
      scale_fill_gradient2(expression('log'[2]*'(FC)'), 
                           low = "blue", 
                           mid = "white", 
                           high = "#ff0000", 
                           midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black",
                                                  ticks = TRUE, 
                                                  ticks.colour = 'black',
                                                  title.position = 'top',
                                                  title.hjust=0.5,
                                                  barwidth = unit(2, "inches")), # increase this for a longer legend
                           breaks=br,
                           labels=la,
                           limits=li)+
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0))+
      xlab(condition_lab)+
      ylab(gene_lab)
  }
  if (tilt_x==45) {
    p_temp = p_temp + theme(axis.text.x=element_text(angle=45,hjust=1, size=text_size))
  }else if(tilt_x==90){
    p_temp = p_temp + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=text_size))
  }else if(tilt_x==0){
    p_temp + theme(axis.text.x=element_text(size=text_size))
  }
  
  # add size control for axis labels, axis titles, legend and geom_text
  p_temp = p_temp + theme(
    text = element_text(size=text_size),
    axis.title = element_text(size=text_size),
    axis.text = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    strip.text = element_text(size=facet_text_size),
    # modify geom_text size
    plot.subtitle = element_text(size=text_size),
    plot.caption = element_text(size=text_size)
  )
  
  # move the legend to the bottom if legend_bottom is TRUE
  if (legend_bottom) {
    p_temp = p_temp + theme(legend.position = "bottom")
  }
  
  if(!is.null(facet)&is.null(facet2)){
    if(!flip_axis){
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_x')
    }else{
      p_temp=p_temp+
        facet_wrap(facet_temp~.,scales = 'free_y')
    }
  }else if(!is.null(facet)&!is.null(facet2)){
    p_temp=p_temp+
      facet_grid(facet_temp~facet_temp2,scales = 'free',space='free')
  }
  return(p_temp)
}

time_degs_plot$group_1=factor(time_degs_plot$group_1,levels = c('20 days',
                                                                '10 days',
                                                                '4 days'))
group_1=unique(time_degs_plot$group)[1:4]
group_1_genes=genes_plot[genes_plot%in%time_degs_plot$ensembl[time_degs_plot$group%in%group_1]]
time_degs_plot_1=time_degs_plot[time_degs_plot$group%in%group_1,]
group_2=unique(time_degs_plot$group)[5:7]
group_2_genes=genes_plot[genes_plot%in%time_degs_plot$ensembl[time_degs_plot$group%in%group_2]]
time_degs_plot_2=time_degs_plot[time_degs_plot$group%in%group_2,]
compare_expression_temporal_1=compare_expression_temporal(degs=time_degs_plot_1,
                            gene_col = 'ensembl',
                            group = 'group_1',
                            genes=unique(group_1_genes),
                            facet = 'cell_type',
                            facet2 = 'group',
                            flip_axis = FALSE,
                            scale_temp = 2,fill_cap = 6,
                            tilt_x = 90,
                            text_size = 13,
                            legend_bottom = TRUE,
                            condition_lab='Days Post-Irradiation',
                            graph_limits = c(-3,4))

compare_expression_temporal_2=compare_expression_temporal(degs=time_degs_plot_2,
                                                          gene_col = 'ensembl',
                                                          group = 'group_1',
                                                          genes=unique(group_2_genes),
                                                          facet = 'cell_type',
                                                          facet2 = 'group',
                                                          flip_axis = FALSE,
                                                          scale_temp = 2,fill_cap = 6,
                                                          tilt_x = 90,
                                                          text_size = 13,
                                                          legend_bottom = TRUE,
                                                          condition_lab='Days Post-Irradiation',
                                                          graph_limits = c(-3,4))

save_p(compare_expression_temporal_1,
       file_name = 'temporal_metabolism_genes_p1',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_figures/',
       p_width = 20)

save_p(compare_expression_temporal_2,
       file_name = 'temporal_metabolism_genes_p2',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_figures/',
       p_width = 20)

pathway_genes=read.csv('/Users/ravelarvargas/Downloads/marian/pathways_genes.csv')
time_degs$group_1=as.character(time_degs$group_1)
time_degs_sig=time_degs[time_degs$sig=='y',]
time_degs_sig$dir_accession=paste0(time_degs_sig$cell_type,'_',time_degs_sig$dir_accession)
time_degs_sig$dir_accession=paste0(time_degs_sig$cell_type,'_',
                                   time_degs_sig$dir_accession)
metabolism_arrest_overlaps_temporal=overlap_function(df_1 = time_degs_sig,
                                            df_2 = pathway_genes,
                                            gene_col_1 = 'ensembl',
                                            gene_col_2 = 'gene',
                                            group_col_1 = 'dir_accession',
                                            carry_col_1 = c('group_1','direction_1','cell_type'),
                                            group_col_2 = 'class',
                                            background = unique(time_degs$ensembl))

save_csv(metabolism_arrest_overlaps_temporal,file_name = 'metabolism_temporal_fishers',
         path = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_tables/')

metabolism_arrest_overlaps_temporal$direction_1[metabolism_arrest_overlaps_temporal$direction_1=='up']='Up\nin Arrest'
metabolism_arrest_overlaps_temporal$direction_1[metabolism_arrest_overlaps_temporal$direction_1=='down']='Down\nin Arrest'
metabolism_arrest_overlaps_temporal$group_1=factor(metabolism_arrest_overlaps_temporal$group_1,
                                                   levels = c('4 days','10 days','20 days'))
overlaps_metabolism_temporal=create_overlap_plot(metabolism_arrest_overlaps_temporal,
                                        x = 'group_1',
                                        y='direction_1',
                                        facet_col = 'cell_type',
                                        facet_2 ='class', x_tilt = 90)

save_p(overlaps_metabolism_temporal,
       file_name = 'metabolism_temporal',
       save_dir = '/Volumes/GoogleDrive/My Drive/PhD_to_publish/systems_analysis_arrest/Final/SI_figures/',
       p_height = 10)
