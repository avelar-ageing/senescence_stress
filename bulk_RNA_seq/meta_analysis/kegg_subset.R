kegg_no_disease=read.csv('/Users/ravelarvargas/Downloads/temp_kegg.csv')

kegg_no_disease[1:5,]

filter_dataframe <- function(dataframe, column_name, terms) {
  # Check if the column exists in the dataframe
  if(!column_name %in% names(dataframe)) {
    stop("Column name not found in the dataframe")
  }
  
  # Initialize a logical vector to keep track of rows to keep
  keep_rows <- rep(TRUE, nrow(dataframe))
  
  # Loop through each term and update the logical vector
  for(term in terms) {
    # Update keep_rows to FALSE for rows where the term is found
    keep_rows <- keep_rows & !grepl(term, dataframe[[column_name]], ignore.case = TRUE)
  }
  
  # Return the filtered dataframe
  return(dataframe[keep_rows, ])
}

kegg_test=filter_dataframe(dataframe=kegg_no_disease,
                           column_name='Description',
                           terms=c('cancer','virus',
                                   'hepatitis','leukemia',
                                   'Measles','diabetes','diabetic','viral',
                                   'atherosclerosis','disease',
                                   'glioma','carcinoma','Melanoma',
                                   'shigellosis','Toxoplasmosis',
                                   'carcinogenesis','infection',
                                   'Influenza','Endocrine resistance',
                                   'EGFR tyrosine kinase inhibitor resistance',
                                   'Leishmaniasis','Cushing',
                                   'Tuberculosis'))

kegg_filter=kegg_test[1:50,]

enrich_kegg_subset=enrichment_dotplot(kegg_filter)

save_p(enrich_kegg_subset,file_name = 'kegg_subset',
       save_dir =  RERUN_DIR,p_width = 10,p_height = 11)
