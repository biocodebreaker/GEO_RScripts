
# Combine Immune DEGs (OverExpressed and UnderExpressed)
GSE62254_Immune_DEGs_combined <- rbind(
  GSE62254_Immune_OverExpressed_positive_logFC_0.05,
  GSE62254_Immune_underexpressed_negative_logFC_0.05
)

# Combine Stromal DEGs (OverExpressed and UnderExpressed)

GSE62254_Stromal_DEGs_combined <- rbind(
  GSE62254_Stromal_OverExpressed_positive_logFC_0.05,
  GSE62254_Stromal_underexpressed_negative_logFC_0.05
)

# Find the intersection with Stromal DEGs
intersection_genes <- intersect(
  rownames(GSE62254_Immune_DEGs_combined),
  rownames(GSE62254_Stromal_DEGs_combined)
)

# Convert to a data frame with one column
GSE62254_stromal_immune_intersection_genes_624 <- data.frame(Intersection_Genes = intersection_genes)

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

# Check if genes are present in the intersection data frame
genes_present <- genes_to_check %in% GSE62254_stromal_immune_intersection_genes_624$Intersection_Genes

# Print the result
genes_present
