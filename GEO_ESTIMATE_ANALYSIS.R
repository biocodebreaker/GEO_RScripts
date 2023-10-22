#ESTIMATE ANALYSIS OF GEO CEL FILES

#RMA 

#Download the raw CEL data from GEO 

#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62254/suppl/GSE62254_RAW.tar

#Create an AffyBatch object from the CEL files

GSE62254_celFiles_AffyBatch <- ReadAffy()

#Create an ExpressionSet object using RMA (Background correcting, Normalizing and Calculating Expression)

GSE62254_ExpressionSet_RMA <- rma(GSE62254_celFiles_AffyBatch, background=TRUE, normalize=TRUE, target="core")

#ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data) 
#Installation

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)


filterCommonGenes("GSE62254_estimate_input_Gene_Symbols_Average_Expression.txt", "GSE62254_estimate_input_Gene_Symbols_Average_Expression.gct")
#[1] "Merged dataset includes 9693 genes (719 mismatched)."

estimateScore("GSE62254_estimate_input_Gene_Symbols_Average_Expression.gct", "GSE62254_Gene_Symbols_Average_Expression_estimate_scores.gct")
#[1] "1 gene set: StromalSignature  overlap= 135"
#[1] "2 gene set: ImmuneSignature  overlap= 136"


filterCommonGenes("GSE62254_Gene_Symbols_Average_Expression_estimate_input_FINAL.txt", "GSE62254_Gene_Symbols_Average_Expression_estimate_input_FINAL.gct")
#[1] "Merged dataset includes 5746 genes (4666 mismatched)."

estimateScore("GSE62254_Gene_Symbols_Average_Expression_estimate_input_FINAL.gct", "GSE62254_Gene_Symbols_Average_Expression_estimate_scores.gct")
#[1] "1 gene set: StromalSignature  overlap= 92"
#[1] "2 gene set: ImmuneSignature  overlap= 99"

