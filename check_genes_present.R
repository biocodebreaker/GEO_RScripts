
# Define the genes to check
genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

# Check if genes are present in the data frame
genes_present <- genes_to_check %in% GSE84433_expression_data_GeneSymbol$Symbol

# Subset the DataFrame to include only rows where genes are present
genes_present_df <- subset(GSE62254_estimate_input_Gene_Symbols_Average_Expression, external_gene_symbol %in% genes_to_check)

# Print the resulting DataFrame
print(genes_present_df[,1:5])




# Genes to check
genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

# Initialize an empty data frame to store the results
result_df <- data.frame()

# Loop through each gene and check if it's present
for (gene_name_to_check in genes_to_check) {
  is_present <- gene_name_to_check %in% GSE62254_estimate_input_probeids_GeneSymbols$external_gene_name
  
  if (is_present) {
    rows_with_gene <- GSE62254_estimate_input_probeids_GeneSymbols[GSE62254_estimate_input_probeids_GeneSymbols$external_gene_name == gene_name_to_check, ]
    result_df <- rbind(result_df, rows_with_gene)
  }
}

# Print the resulting data frame
print(result_df)
