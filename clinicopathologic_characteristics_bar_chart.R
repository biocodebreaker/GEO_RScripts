# Load ggplot2
#library(ggplot2)

# Create a data frame with subgroup names and totals
subgroups <- c("Age <60", "Age >60", "Female", "Male", "T1-T2", "T3-T4", "N0", "N1-N3", "M0", "M1", "Stage I-II", "Stage III-IV")
totals <- c(131, 274, 142, 263, 106, 291, 122, 264, 360, 25, 177, 204)
data <- data.frame(Subgroup = factor(subgroups, levels = subgroups), Total = totals)

# Define custom colors for subgroups within a group
group_colors <- c("#1f77b4", "#1f77b4", "#ff7f0e", "#ff7f0e", "#2ca02c", "#2ca02c",
                  "#d62728", "#d62728", "#9467bd", "#9467bd", "#8c564b", "#8c564b")

# Create the bar chart with custom colors and without title and grid lines
ggplot(data, aes(x = Subgroup, y = Total, fill = Subgroup)) +
  geom_bar(stat = "identity") +
  xlab("Clinicopathologic characteristics") +
  ylab("subgroup sample size") +
  scale_fill_manual(values = group_colors) +  # Assign custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_blank(),  # Remove title
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
