############################################################
# 03_interpretable_ranking.R
#
# Purpose:
# Generate interpretable AI target rankings and feature
# importance estimates from the multi-modal feature matrix.
############################################################

library(dplyr)
library(ggplot2)

feature_matrix <- read.csv(
  "results/Final_AI_Target_Feature_Matrix.csv"
)

############################################################
# Final target ranking plot
############################################################

png(
  "figures/Figure4_Weighted_AI_Target_Ranking.png",
  width = 2200,
  height = 1600,
  res = 300
)

ggplot(
  feature_matrix,
  aes(
    x = reorder(Gene, Weighted_AI_Score),
    y = Weighted_AI_Score
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "AI-driven therapeutic target ranking in lupus nephritis",
    x = "Target",
    y = "Weighted AI score"
  )

dev.off()

list.files("figures")

############################################################
# Feature contribution weights
############################################################

feature_weights <- data.frame(
  Feature = c(
    "avg_log2FC",
    "Specificity_Score",
    "Spatial_FC",
    "Spatial_Score",
    "Spearman_rho",
    "Correlation_Score",
    "Mouse_Validated",
    "Communication_Count",
    "Pathway_Count"
  ),
  Weight = c(
    1.0,
    1.5,
    2.0,
    2.0,
    1.5,
    1.5,
    2.0,
    1.0,
    1.0
  )
)

write.csv(
  feature_weights,
  "results/AI_Feature_Weights.csv",
  row.names = FALSE
)

png(
  "figures/Figure5_AI_Feature_Weights.png",
  width = 2200,
  height = 1600,
  res = 300
)

ggplot(
  feature_weights,
  aes(
    x = reorder(Feature, Weight),
    y = Weight
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Feature weights used for interpretable AI target ranking",
    x = "Feature",
    y = "Weight"
  )

dev.off()

list.files("figures")
list.files("results")