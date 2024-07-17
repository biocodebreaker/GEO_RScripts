BiocManager::install("GEOquery")

library(GEOquery)

library("dplyr")

library("tidyverse")

# GSE147163 ANALYSIS 

GSE147163_matrix <- getGEO("GSE147163", GSEMatrix = TRUE)


GSE54129_matrix <- getGEO("GSE54129", GSEMatrix = TRUE)

GSE15081_matrix <- getGEO("GSE15081", GSEMatrix = TRUE)

GSE84426_matrix <- getGEO("GSE84426", GSEMatrix = TRUE)

GSE15459_matrix <- getGEO("GSE15459", GSEMatrix = TRUE)

GSE34942_matrix <- getGEO("GSE34942", GSEMatrix = TRUE)

GSE57303_matrix <- getGEO("GSE57303", GSEMatrix = TRUE)

GSE51105_matrix <- getGEO("GSE51105", GSEMatrix = TRUE)

# GSE33335 has clinicopathologic x-tics (age, gender, TNM) but not surv info (time, event)

GSE33335_matrix <- getGEO("GSE33335", GSEMatrix = TRUE)

GSE33335_SOFT_format <- getGEO("GSE33335", GSEMatrix = FALSE )

GSE103236_matrix <- getGEO("GSE103236", GSEMatrix = TRUE)

GSE15456_matrix <- getGEO("GSE15456", GSEMatrix = TRUE)
GSE15456_matrix_phenoData <- pData(GSE15456_matrix[[1]])

# GSE26899 has clinicopathologic info but not surv info 
GSE26899_matrix <- getGEO("GSE26899", GSEMatrix = TRUE)
GSE26899_matrix_phenoData <- pData(GSE26899_matrix[[1]])

# GSE27342 has clinicopathologic but not survival info 
GSE27342_matrix <- getGEO("GSE27342", GSEMatrix = TRUE)
GSE27342_matrix_phenoData <- pData(GSE27342_matrix[[1]])

# GSE26942 has no surv info but has clinicopathologic 
GSE26942_matrix <- getGEO("GSE26942", GSEMatrix = TRUE)
GSE26942_matrix_phenoData <- pData(GSE26942_matrix[[1]])

# Does not have clinicopathologic and surv info 
GSE66229_matrix <- getGEO("GSE66229", GSEMatrix = TRUE)
GSE66229_matrix_phenoData <- pData(GSE66229_matrix[[1]])



#Get Clinical/phenotypic data

GSE147163_phenoData <- pData(GSE147163_matrix[[1]])

GSE54129_phenoData <- pData(GSE54129_matrix[[1]])

GSE15081_phenoData <- pData(GSE15081_matrix[[1]])

# GSE84426_matrix sample size 76, surv info yes, clinical info yes
GSE84426_phenoData <- pData(GSE84426_matrix[[1]])

GSE15459_phenoData <- pData(GSE15459_matrix[[1]])

GSE34942_phenoData <- pData(GSE34942_matrix[[1]])


GSE51105_phenoData <- pData(GSE51105_matrix[[1]])


GSE33335_phenoData <- pData(GSE33335_matrix[[1]])

GSE103236_phenoData <- pData(GSE103236_matrix[[1]])


# GSE147163 does not have survival information (time, event) and no clinicopatholigical characteristics

# GSE62254 doesn't have clinical data and no surv info


GSE147163_phenoData_subset <- pData(GSE147163_matrix[[1]]) %>%
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

exprs(GSE84433_matrix[[1]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") -> GSE84433_expression_data

#Get Annotation/ Feature data
GSE84433_featureData <- fData(GSE84433_matrix[[1]])

#subset the GSE84433_featureData data frame to include only the columns "ID", "Entrez_Gene_ID" and "Symbol" 
GSE84433_featureData_subset <- GSE84433_featureData %>%
  select(ID, Entrez_Gene_ID, Symbol)


GSE84433_expression_data %>%
  inner_join(GSE84433_featureData_subset) -> GSE84433_expression_data_GeneSymbol_original



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

# Calculate Average Expression per Gene Symbol using the aggregate function

# Remove the first column (GeneSymbol) for aggregation
df_for_aggregation <- GSE84433_expression_data_GeneSymbol_filtered[, -1]

# Calculate the mean for each gene symbol
gene_symbol_means <- aggregate(. ~ GSE84433_expression_data_GeneSymbol_filtered$Symbol, data = df_for_aggregation, mean)

# Rename the columns for clarity
colnames(gene_symbol_means) <- c("Symbol", colnames(df_for_aggregation))

# Print the resulting DataFrame
#gene_symbol_means

# Calculate Median Expression per Gene Symbol using the aggregate function

# Remove the first column (GeneSymbol) for aggregation
df_for_aggregation <- GSE84433_expression_data_GeneSymbol_filtered[, -1]

# Calculate the median for each gene symbol

gene_symbol_medians <- aggregate(. ~ GSE84433_expression_data_GeneSymbol_filtered$Symbol, data = df_for_aggregation, median)

# Rename the columns for clarity
colnames(gene_symbol_medians) <- c("Symbol", colnames(df_for_aggregation))

GSE84433_Gene_Symbols_Median_Expression_estimate_input <- gene_symbol_medians

# Send to Excel to edit for ESTIMATE input 

write.table(GSE84433_Gene_Symbols_Median_Expression_estimate_input, "GSE84433_Gene_Symbols_Median_Expression_estimate_input.txt", sep = "\t", row.names = TRUE)

# Read in the edited file for ESTIMATE input

GSE84433_Gene_Symbols_Estimate_input_Median_Expression <- read.table("GSE84433_Gene_Symbols_Estimate_input_Median_Expression.txt", header = TRUE, sep = "\t")

library(estimate)

# ESTIMATE ANALYSIS BASED ON MEDIAN EXPRESSION 
filterCommonGenes("GSE84433_Gene_Symbols_Estimate_input_Median_Expression.txt", "GSE84433_Gene_Symbols_Estimate_input_Median_Expression.gct")

#[1] "Merged dataset includes 10280 genes (132 mismatched)."

estimateScore("GSE84433_Gene_Symbols_Estimate_input_Median_Expression.gct", "GSE84433_Gene_Symbols_Median_Expression_estimate_scores.gct")

#[1] "1 gene set: StromalSignature  overlap= 140"
#[1] "2 gene set: ImmuneSignature  overlap= 141"

# ESTIMATE ANALYSIS BASED ON AVERAGE EXPRESSION

filterCommonGenes("GSE84433_Gene_Symbols_Average_Expression_estimate_input.txt", "GSE84433_Gene_Symbols_Average_Expression_estimate_input.gct")
#[1] "Merged dataset includes 10280 genes (132 mismatched)."

estimateScore("GSE84433_Gene_Symbols_Average_Expression_estimate_input.gct", "GSE84433_Gene_Symbols_Average_Expression_estimate_scores.gct")
#"1 gene set: StromalSignature  overlap= 140"
#"2 gene set: ImmuneSignature  overlap= 141"


# Transpose the estimate scores .gct file such that the samples are rows and Gene Symbols are columns
# This is necessary befofe merging with phenoData subset (survival info)
# Transpose using excel and rename sample column to "sample_GSM" to match phenoData subset


#join phenoData with estimate scores 
GSE84433_Gene_Symbols_Median_Expression_estimateScores %>%
  inner_join(GSE84433_phenoData_subset) -> GSE84433_Median_Expression_estimateScores_phenoData

# This doesn't work because you have to transpose the expression data first before merge with phenoData (surv info)
#GSE84433_expression_data_GeneSymbol %>%
# inner_join(GSE84433_phenoData_subset ) -> GSE84433_expression_surv_info

# GSE84433_Gene_Symbols_Estimate_input_Median_Expression data.frame has Gene Symbols as rows and samples as columns so it needs to be transposed
# such that the samples become rows and Gene Symbols become columns
# However transposing converts data.frame to array matrix so there is need to remove the V1, V2, V3 column names added


library(tidyverse)

# Transpose the data frame and convert it back to a data frame
GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed <- as.data.frame(t(GSE84433_Gene_Symbols_Estimate_input_Median_Expression))


# Check dimensions of transposed data.frame 
dim(GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed)
#[1]   358 25124

# Excel has a limit of 16,384 columns and 1,048,576 rows
# Split the columns so that you get two equal files 
# While in Excel, remove row with V1, V2, V3.... Rename "Gene.Symbol" column to "sample_GSM"

GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1 <- GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed[, 1:12562]

GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2 <- GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed[, 12563:25124]

# Write the files to Excel for proper editing 
write.table(GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1, "GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1.txt", sep = "\t", row.names = TRUE)

write.table(GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2, "GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2.txt", sep = "\t", row.names = TRUE)

GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1 <- read.table("GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1.txt", header = TRUE, sep = "\t")

library(dplyr)

# First inner join between GSE84433_Median_Expression_estimateScores_phenoData and part1
result_part1 <- inner_join(GSE84433_Median_Expression_estimateScores_phenoData, GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part1, by = "sample_GSM")

# Second inner join between the result from part1 and part2
GSE84433_Median_Expression_surv_info_estimateScores <- inner_join(result_part1, GSE84433_Gene_Symbols_Estimate_input_Median_Expression_transposed_part2, by = "sample_GSM")

dim(GSE84433_Median_Expression_surv_info_estimateScores)
#[1]   357 25135

library(survival)
library(maxstat)

# NOTE: Find out if survival time is in days or months. Presumption: It is in days. 

GSE84433_Median_Expression_surv_info_estimateScores_surv_obj <- with(GSE84433_Median_Expression_surv_info_estimateScores, Surv(as.numeric(time), as.numeric(event)))


GSE84433_Median_Expr_estimateScores_maxstat_stromal_optimalcutoff <- maxstat.test(GSE84433_Median_Expression_surv_info_estimateScores_surv_obj  ~ GSE84433_Median_Expression_surv_info_estimateScores$StromalScore, data = GSE84433_Median_Expression_surv_info_estimateScores, smethod = "LogRank")

#estimated cutpoint 419.5145

GSE84433_Median_Expr_estimateScores_maxstat_IMMUNE_optimalcutoff <- maxstat.test(GSE84433_Median_Expression_surv_info_estimateScores_surv_obj  ~ GSE84433_Median_Expression_surv_info_estimateScores$ImmuneScore, data = GSE84433_Median_Expression_surv_info_estimateScores, smethod = "LogRank")
#estimated cutpoint 1900.368 

# Optimal cutoff values
cutoff_stromal <- 419.5145
cutoff_immune <- 1900.368

# Create new columns for categorizing StromalScore and ImmuneScore
GSE84433_Median_Expression_surv_info_estimateScores$StromalScore_group <- ifelse(
  GSE84433_Median_Expression_surv_info_estimateScores$StromalScore >= cutoff_stromal,
  "high",
  "low"
)

GSE84433_Median_Expression_surv_info_estimateScores$ImmuneScore_group <- ifelse(
  GSE84433_Median_Expression_surv_info_estimateScores$ImmuneScore >= cutoff_immune,
  "high",
  "low"
)

library(dplyr)

# Rearrange the data.frame
GSE84433_Median_Expression_surv_info_estimateScores_categorized <- GSE84433_Median_Expression_surv_info_estimateScores %>%
  select(1:4, StromalScore_group, ImmuneScore_group, 5:ncol(.))


# Convert time and event to appropriate data types
GSE84433_Median_Expression_surv_info_estimateScores_categorized$time <- as.numeric(GSE84433_Median_Expression_surv_info_estimateScores_categorized$time)
GSE84433_Median_Expression_surv_info_estimateScores_categorized$event <- as.numeric(GSE84433_Median_Expression_surv_info_estimateScores_categorized$event)


# Create a Surv object
surv_obj <- Surv(time = GSE84433_Median_Expression_surv_info_estimateScores_categorized$time, event = GSE84433_Median_Expression_surv_info_estimateScores_categorized$event)


#High vs Low ImmuneScore_group survival analysis 

immune_km_fit <- survfit(surv_obj ~ ImmuneScore_group, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized)

library(survminer)

survminer::ggsurvplot(immune_km_fit, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized, risk.table = T, pval = TRUE)


#High vs Low StromalScore_group survival analysis 
stromal_km_fit <- survfit(surv_obj ~ StromalScore_group, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized)

survminer::ggsurvplot(stromal_km_fit, data = GSE84433_Median_Expression_surv_info_estimateScores_categorized, risk.table = T, pval = TRUE)

#count number of high vs low
GSE84433_immune_counts <- table(GSE84433_Median_Expression_surv_info_estimateScores_categorized$ImmuneScore_group)

GSE84433_stromal_counts <- table(GSE84433_Median_Expression_surv_info_estimateScores_categorized$StromalScore_group)


# Load the limma library
library(limma)

# Create a factor variable for high and low ImmuneScore

ImmuneScore_group <- factor(GSE84433_Median_Expression_surv_info_estimateScores_categorized$ImmuneScore_group, levels = c("high", "low"))

# Create the design matrix for ImmuneScore_group 
immune_design_matrix <- model.matrix(~0 + ImmuneScore_group)

# Create contrasts using the immune_design_matrix
immune_contrast_matrix <- makeContrasts(high_vs_low = ImmuneScore_grouphigh - ImmuneScore_grouplow, levels = colnames(immune_design_matrix))


# Extract gene expression data from the data frame # Exclude non-gene columns

GSE84433_Median_Expression_lmFit <- GSE84433_Median_Expression_surv_info_estimateScores_categorized[, -(1:13)] 

# Transpose it such that the row dimension of design matrix matches column dimension of expression data
GSE84433_Median_Expression_lmFit <- t(GSE84433_Median_Expression_lmFit)

# Perform differential expression analysis
immune_lmFit <- lmFit(GSE84433_Median_Expression_lmFit, immune_design_matrix)

fit_contrast_immune <- contrasts.fit(immune_lmFit, immune_contrast_matrix)

immune_eBayes <- eBayes(fit_contrast_immune)

# Extract the differentially expressed genes

GSE84433_Immune_DEGs <- topTable(immune_eBayes, number = Inf)

# Get all significant Immune DEGs

GSE84433_Immune_DEGs_adj.P.Val_0.05 <- GSE84433_Immune_DEGs[GSE84433_Immune_DEGs$adj.P.Val < 0.05, ]

# Define the genes you want to check

genes_to_check <- c("MAGEA11", "CNKSR2", "SELP", "CYP1B1", "OMD", "CNR1", "CCR4", "FMO2")

# Check if the genes are present in the row names of GSE84433_Stromal_DEGS

matching_rows <- rownames(GSE84433_Immune_DEGs) %in% genes_to_check

# Print the entire rows where genes are present
print(GSE84433_Immune_DEGs[matching_rows, ])

# StromalScore_group ANALYSIS 

StromalScore_group <- factor(GSE84433_Median_Expression_surv_info_estimateScores_categorized$StromalScore_group, levels = c("high", "low"))

stromal_design_matrix <- model.matrix(~0 + StromalScore_group)

stromal_contrast_matrix <- makeContrasts(high_vs_low = StromalScore_grouphigh - StromalScore_grouplow, levels = colnames(stromal_design_matrix))

stromal_lmFit <- lmFit(GSE84433_Median_Expression_lmFit, stromal_design_matrix)

fit_contrast_stromal <- contrasts.fit(stromal_lmFit, stromal_contrast_matrix)

stromal_eBayes <- eBayes(fit_contrast_stromal)

GSE84433_Stromal_DEGs <- topTable(stromal_eBayes, number = Inf)

GSE84433_Stromal_DEGs_adj.P.Val_0.05 <- GSE84433_Stromal_DEGs[GSE84433_Stromal_DEGs$adj.P.Val < 0.05, ]



GSE84433_intersection_genes <- intersect(rownames(GSE84433_Immune_DEGs_adj.P.Val_0.05), rownames(GSE84433_Stromal_DEGs_adj.P.Val_0.05))
length(GSE84433_intersection_genes)
# [1] 2263

# UNIVARIATE ANALYSIS 

# Prepare the data
GSE84433_univariate_analysis_input_clean_data <- GSE84433_Median_Expression_surv_info_estimateScores_categorized[!is.na(GSE84433_Median_Expression_surv_info_estimateScores_categorized$time), ]
GSE84433_univariate_analysis_input_clean_data <- GSE84433_univariate_analysis_input_clean_data[GSE84433_univariate_analysis_input_clean_data$time > 0, ]

# Extract survival data
GSE84433_survival_data <- GSE84433_univariate_analysis_input_clean_data[, c("time", "event")]


# Create the Surv() object
surv_obj <- Surv(GSE84433_survival_data$time, GSE84433_survival_data$event)

# Initialize an empty list to store results
GSE84433_Median_Expression_univariate_results <- list()

# Perform univariate Cox regression for each gene in intersection_genes
for (gene in GSE84433_intersection_genes) {
  cox_model <- coxph(surv_obj ~ get(gene), data = GSE84433_univariate_analysis_input_clean_data)
  result <- summary(cox_model)
  GSE84433_Median_Expression_univariate_results[[gene]] <- result
}

# Can't use all the columns, predictors too many (2263)
#GSE84433_Median_Expression_multivariate_results <- coxph(surv_obj ~ ., data = GSE84433_univariate_analysis_input_clean_data)  # Only include predictor variables


library(dplyr)

#rearranges the columns so that CCR4, FMO2, and SELP come first
GSE84433_univariate_analysis_input_clean_data <- GSE84433_univariate_analysis_input_clean_data %>%
  select(sample_GSM, StromalScore, ImmuneScore, ESTIMATEScore, 
         StromalScore_group, ImmuneScore_group, TumorPurity, age, event, time,
         Nstage, Tstage, Gender, CCR4, FMO2, SELP, everything())

# Perform multivariate Cox regression using columns 14:100

GSE84433_Median_Expression_multivariate_results <- coxph(surv_obj ~ ., data = GSE84433_univariate_analysis_input_clean_data[, 14:100])


# Perform multivariate Cox regression with SELP, CCR4, and FMO2 as predictors
multivariate_cox_model <- coxph(
  surv_obj ~ SELP + CCR4 + FMO2,
  data = GSE84433_univariate_analysis_input_clean_data
)

# Get the summary of the multivariate Cox regression
#GSE84433_Median_Expression_multivariate_results <- summary(multivariate_cox_model)
