#Mapping Affymetrix Probesets to Gene Symbols

#NOTE entrezgene_id attribute is not available with affy_hg_u133_plus_2 filter

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
GSE62254_estimate_input_probeids_GeneSymbols <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = GSE62254_series_matrix_simple_estimate_input_df$X, uniqueRows=TRUE)


# NOT A GOOD IDEA produces error. Remove rows with empty external_gene_name values
GSE62254_estimate_input_probeids_GeneSymbols_filtered <- GSE62254_estimate_input_probeids_GeneSymbols %>%
  filter(external_gene_name != "")


matching_indices <- match(rownames(GSE62254_series_matrix_simple_estimate_input_df), GSE62254_estimate_input_probeids_GeneSymbols$affy_hg_u133_plus_2)

rownames(GSE62254_series_matrix_simple_estimate_input_df) <- paste(GSE62254_estimate_input_probeids_GeneSymbols[matching_indices, "external_gene_name"])

merged_df <- data.frame(rownames(GSE62254_series_matrix_simple_estimate_input_df), GSE62254_estimate_input_probeids_GeneSymbols[matching_indices, c("affy_hg_u133_plus_2", "external_gene_name")])


rownames(GSE62254_series_matrix_simple_estimate_input_df) <- paste(GSE62254_estimate_input_probeids_GeneSymbols[matching_indices, "external_gene_name"], c(1:length(matching_indices)), sep="_")


#Get average expression per gene symbol using avereps in limma

GSE62254_estimate_input_Gene_Symbols_Average_Expression <- data.frame (avereps(GSE62254_estimate_input_Gene_Symbols_notUnique, ID=GSE62254_estimate_input_Gene_Symbols_notUnique[1:22112,1]))

