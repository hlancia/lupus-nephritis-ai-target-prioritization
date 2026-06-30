############################################################
# 04_project_summary_outputs.R
#
# Purpose:
# Export concise summary tables for reporting and GitHub.
############################################################

library(dplyr)

feature_matrix <- read.csv(
  "results/Final_AI_Target_Feature_Matrix.csv"
)

top_targets <- feature_matrix %>%
  select(
    Gene,
    Weighted_AI_Score,
    AI_Target_Score,
    avg_log2FC,
    Spatial_FC,
    Spearman_rho,
    Mouse_Validated,
    Communication_Count,
    Pathway_Count
  ) %>%
  arrange(desc(Weighted_AI_Score))

write.csv(
  top_targets,
  "results/Top_AI_Prioritized_Targets.csv",
  row.names = FALSE
)

top_targets

list.files("results")