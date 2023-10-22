library("tidyverse")
library("dplyr")
#GSE84433 ANALYSIS 

GSE84433 <- getGEO("GSE84433", GSEMatrix = TRUE)

#Get Clinical/phenotypic data

#GSE84433_phenoData <- pData(GSE84433[[1]])


GSE84433_phenoData_subset <- pData(GSE84433[[1]]) %>%
  as_tibble() %>%
  select(sample_GSM = geo_accession,
         age = `age:ch1`,
         event = `death:ch1`,
         time = `duration overall survival:ch1`,
         Nstage = `pnstage:ch1`,
         Tstage = `ptstage:ch1`,
         Gender = `Sex:ch1`)


#Get ExpressionSet data from GSE84433$GSE84433_series_matrix.txt.gz

#GSE84433_expression_data <- exprs(GSE84433[[1]])


#Join the Feature Data with the Expression data
#First transform the expression matrix to a data.frame 
#library("tidyverse")

exprs(GSE84433[[1]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") -> GSE84433_expression_data

#Get Annotation/ Feature data
GSE84433_featureData <- fData(GSE84433[[1]])

#subset the GSE84433_featureData data frame to include only the columns "ID", "Entrez_Gene_ID" and "Symbol" 
GSE84433_featureData_subset <- GSE84433_featureData %>%
  select(ID, Entrez_Gene_ID, Symbol)


GSE84433_expression_data %>%
  inner_join(GSE84433_featureData_subset) -> GSE84433_expression_data_GeneSymbol_original

#join phenoData with estimate scores 
GSE84433_estimateScores %>%
  inner_join(GSE84433_phenoData_subset) -> GSE84433_estimateScores_phenoData


GSE84433_expression_data_GeneSymbol %>%
  inner_join(GSE84433_phenoData_subset ) -> GSE84433_expression_surv_info


#Check if genes to be validated are present 

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

genes_to_check %in% GSE84433_expression_data_GeneSymbol$Symbol
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# Rearrange the columns
GSE84433_expression_data_GeneSymbol_arranged <- GSE84433_expression_data_GeneSymbol_original %>%
  select(Symbol, Entrez_Gene_ID, everything())

# Sort the data frame by the "Symbol" column
GSE84433_expression_data_GeneSymbol_sorted <- GSE84433_expression_data_GeneSymbol_arranged %>%
  arrange(Symbol)

# Remove rows with blank or NA values in Symbol column
GSE84433_expression_data_GeneSymbol_filtered <- GSE84433_expression_data_GeneSymbol_sorted %>%
  filter(!is.na(Symbol) & Symbol != "")


# First remove the Entrez_Gene_ID column
GSE84433_expression_data_GeneSymbol_filtered <- GSE84433_expression_data_GeneSymbol_filtered %>%
  select(-Entrez_Gene_ID)

GSE84433_expression_data_GeneSymbol_filtered <- GSE84433_expression_data_GeneSymbol_filtered %>%
  select(-ID)

#library("tidyverse")

# Remove the first column (external_gene_symbol) for aggregation
df_for_aggregation <- GSE84433_expression_data_GeneSymbol_filtered[, -1]

# Calculate the mean for each gene symbol
gene_symbol_means <- aggregate(. ~ GSE84433_expression_data_GeneSymbol_filtered$Symbol, data = df_for_aggregation, mean)

# Rename the columns for clarity
colnames(gene_symbol_means) <- c("Symbol", colnames(df_for_aggregation))

# Print the resulting DataFrame
#gene_symbol_means

GSE84433_Gene_Symbols_Average_Expression_estimate_input <- gene_symbol_means

#GSE84433_Gene_Symbols_Average_Expression_estimate_input <- GSE84433_Gene_Symbols_Average_Expression_estimate_input %>%
  #rename("Gene Symbol" = Symbol)

#ESTIMATE ANALYSIS

filterCommonGenes("GSE84433_Gene_Symbols_Average_Expression_estimate_input.txt", "GSE84433_Gene_Symbols_Average_Expression_estimate_input.gct")
#[1] "Merged dataset includes 10280 genes (132 mismatched)."

estimateScore("GSE84433_Gene_Symbols_Average_Expression_estimate_input.gct", "GSE84433_Gene_Symbols_Average_Expression_estimate_scores.gct")
#"1 gene set: StromalSignature  overlap= 140"
#"2 gene set: ImmuneSignature  overlap= 141"



