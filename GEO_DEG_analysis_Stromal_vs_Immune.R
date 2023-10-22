
# Divide the StromalScore and ImmuneScore into the high and low risk groups using optimal cutoff
#Use maxstat cutoffs like in the TCGA dataset
GSE62254_Average_Expression_estimateScores_categorized <- GSE62254_Average_Expression_estimateScores %>%
  mutate(StromalScore_group = ifelse(StromalScore > 1323.487, "high", "low"),
         ImmuneScore_group = ifelse(ImmuneScore > 1957.582, "high", "low"))


# Calculate the median of StromalScore
median_StromalScore <- median(GSE62254_Average_Expression_estimateScores$StromalScore)
# Calculate the median of ImmuneScore
median_ImmuneScore <- median(GSE62254_Average_Expression_estimateScores$ImmuneScore)

# Divide the scores into the high and low groups using median
GSE62254_Average_Expression_estimateScores_categorized <- GSE62254_Average_Expression_estimateScores %>%
  mutate(StromalScore_group = ifelse(StromalScore > median_StromalScore, "high", "low"),
         ImmuneScore_group = ifelse(ImmuneScore > median_ImmuneScore, "high", "low"))

#count number of high vs low
GSE62254_immune_counts <- table(GSE62254_Average_Expression_estimateScores_categorized$ImmuneScore_group)

GSE62254_stromal_counts <- table(GSE62254_Average_Expression_estimateScores_categorized$StromalScore_group)
