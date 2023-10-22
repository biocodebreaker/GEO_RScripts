library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Create a vector of unique ensembl_gene_id values
ensembl_gene_ids <- unique(result_df$ensembl_gene_id)

# Get corresponding entrezgene_id values
genes <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = ensembl_gene_ids,
  mart = mart
)

# Merge the gene information with result_df based on ensembl_gene_id
result_df <- merge(result_df, genes, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
