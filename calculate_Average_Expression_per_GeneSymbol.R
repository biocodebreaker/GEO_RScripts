# Remove rows with blank or NA values in external_gene_symbol column
filtered_df <- subset(GSE62254_series_matrix_simple_estimate_input_df_with_gene_symbols, 
                      !is.na(external_gene_symbol) & external_gene_symbol != "")

# Check the dimensions of the filtered dataframe
dim(filtered_df)

# Sort filtered_df by external_gene_symbol in A-Z order
filtered_df_sorted <- filtered_df[order(filtered_df$external_gene_symbol), ]



# Remove the first column (external_gene_symbol) for aggregation
df_for_aggregation <- filtered_df_sorted[, -1]

# Calculate the mean for each gene symbol
gene_symbol_means <- aggregate(. ~ filtered_df_sorted$external_gene_symbol, data = df_for_aggregation, mean)

# Rename the columns for clarity
colnames(gene_symbol_means) <- c("external_gene_symbol", colnames(df_for_aggregation))

# Print the resulting DataFrame
gene_symbol_means
