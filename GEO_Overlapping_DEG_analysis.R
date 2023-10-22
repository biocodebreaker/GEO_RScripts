#GEO Limma DEG analysis

##create the design matrix for both stromal and immune groups

Stromal_design <- data.frame(
  Stromallow = as.numeric(GSE62254_Average_Expression_estimateScores_categorized$StromalScore_group == "low"),
  Stromalhigh = as.numeric(GSE62254_Average_Expression_estimateScores_categorized$StromalScore_group == "high")
)
rownames(Stromal_design) <- GSE62254_Average_Expression_estimateScores_categorized$SampleName

Immune_design <- data.frame(
  Immunelow = as.numeric(GSE62254_Average_Expression_estimateScores_categorized$ImmuneScore_group == "low"),
  Immunehigh = as.numeric(GSE62254_Average_Expression_estimateScores_categorized$ImmuneScore_group == "high")
)
rownames(Immune_design) <- GSE62254_Average_Expression_estimateScores_categorized$SampleName


# Create a contrast matrix for the Stromal_group
Stromal_contrasts <- makeContrasts(high_vs_low = Stromalhigh - Stromallow, levels = Stromal_design)

# Create a contrast matrix for the Immune_group
Immune_contrasts <- makeContrasts(high_vs_low = Immunehigh - Immunelow, levels = Immune_design)

#set the rownames of the design matrices to be the same as the submitter_id column
#in the merged_survival_info_estimate_scores_gene_expression_limma

rownames(Stromal_design) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id
rownames(Immune_design) <- merged_survival_info_estimate_scores_gene_expression_limma$submitter_id

#make sure that the gene expression data also has the SampleName as its rownames
GSE62254_expression_data <- GSE62254_Average_Expression_estimateScores_categorized[, 6:16388]

rownames(GSE62254_expression_data) <- GSE62254_Average_Expression_estimateScores_categorized$SampleName



#The lmFit function expects the number of columns in the expresion data to match the number of rows in the design matrix.

# Transpose the gene expression data
GSE62254_expression_data_transposed_for_lmfit <- t(GSE62254_expression_data)

# Run lmFit for the Stromal_group
Stromal_fit <- lmFit(GSE62254_expression_data_transposed_for_lmfit, Stromal_design)

# Apply contrasts and compute empirical Bayes statistics
Stromal_fit <- contrasts.fit(Stromal_fit, Stromal_contrasts)
Stromal_fit <- eBayes(Stromal_fit)

# Linear modeling and empirical Bayes moderation for Immune_group
Immune_fit <- lmFit(GSE62254_expression_data_transposed_for_lmfit, Immune_design)
Immune_fit <- contrasts.fit(Immune_fit, Immune_contrasts)
Immune_fit <- eBayes(Immune_fit)

# Find top DEGs for Stromal_group
Stromal_topTable <- topTable(Stromal_fit, n = Inf, adjust.method = "fdr")

# Find top DEGs for Immune_group
Immune_topTable <- topTable(Immune_fit, n = Inf, adjust.method = "fdr")

# Filter DEGs based on FDR-adjusted p-value and log fold-change thresholds
GSE62254_Stromal_DEGs <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.0001 & abs(Stromal_topTable$logFC) > 4, ]




#There were no DEGs with abs logFC > 1
GSE62254_Immune_DEGs <- Immune_topTable[Immune_topTable$adj.P.Val < 0.05 & abs(Immune_topTable$logFC) > 1, ]

Stromal_DEGs_p_value_0.05_logFC_1 <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.05 & abs(Stromal_topTable$logFC) > 1, ]

#4122 Overexpressed Stromal DEGs had abs logFC > 0 (positive logFC)and adj.P.Val < 0.05

GSE62254_Stromal_overexpressed <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.05 & abs(Stromal_topTable$logFC) > 0, ]

#683 Overexpressed Stromal DEGs had abs logFC > 0 (positive logFC)and adj.P.Val < 0.0001

GSE62254_Stromal_overexpressed_0.0001 <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.0001 & abs(Stromal_topTable$logFC) > 0, ]

#2285 Overexpressed Immune DEGs had abs logFC > 0 (positive logFC)and adj.P.Val < 0.05
GSE62254_Immune_overexpressed <- Immune_topTable[Immune_topTable$adj.P.Val < 0.05 & abs(Immune_topTable$logFC) > 0, ]

#64 Overexpressed Immune DEGs adj.P.Val < 0.0001
GSE62254_Immune_overexpressed_0.0001 <- Immune_topTable[Immune_topTable$adj.P.Val < 0.0001 & abs(Immune_topTable$logFC) > 0, ]

#There were no Underexpressed Stromal DEGs (negative logFC)

GSE62254_Stromal_underexpressed_0.0001 <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.0001 & abs(Stromal_topTable$logFC) < 0, ]

GSE62254_Stromal_underexpressed_0.05 <- Stromal_topTable[Stromal_topTable$adj.P.Val < 0.05 & abs(Stromal_topTable$logFC) < 0, ]


#There were no Underexpressed Immune DEGs (negative logFC)

GSE62254_Immune_underexpressed_0.0001 <- Immune_topTable[Immune_topTable$adj.P.Val < 0.0001 & abs(Immune_topTable$logFC) < 0, ]

GSE62254_Immune_underexpressed_0.05 <- Immune_topTable[Immune_topTable$adj.P.Val < 0.05 & abs(Immune_topTable$logFC) < 0, ]

#Find overlapping DEGs between Stromal and Immune Groups 
#There were no overlapping DEGs abs.logFC > 0,adj.P.Val < 0.05 

GSE62254_intersection_genes <- intersect(rownames(GSE62254_Immune_overexpressed), rownames(GSE62254_Stromal_overexpressed))

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

dim(GSE62254_Stromal_overexpressed)
#[1] 4122    6

present_genes <- genes_to_check[genes_to_check %in% rownames(GSE62254_Stromal_overexpressed)]
present_genes
#[1] "MAGEA11"

dim(GSE62254_Immune_overexpressed)
#[1] 2285    6
#> present_genes <- genes_to_check[genes_to_check %in% rownames(GSE62254_Immune_overexpressed)]
#> present_genes
#[1] "OMD"

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")
present_genes <- genes_to_check[genes_to_check %in% rownames(GSE62254_Stromal_overexpressed)]
rows_with_genes <- GSE62254_Stromal_overexpressed[present_genes, ]


#          logFC  AveExpr        t     P.Value  adj.P.Val         B
#MAGEA11 -0.04432287 2.739737 -2.79154 0.005581871 0.02708761 -3.249327