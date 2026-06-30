############################################################
# 02_unsupervised_learning.R
#
# Project:
# AI-driven therapeutic target prioritization in lupus nephritis
#
# Purpose:
# Use unsupervised learning to visualize and cluster
# therapeutic target candidates based on multi-modal
# features.
############################################################

library(dplyr)
library(ggplot2)
library(pheatmap)

############################################################
# Load feature matrix
############################################################

feature_matrix <- read.csv(
  "results/Final_AI_Target_Feature_Matrix.csv"
)

rownames(feature_matrix) <- feature_matrix$Gene

############################################################
# Select ML features
############################################################

ml_features <- feature_matrix %>%
  select(
    avg_log2FC,
    Specificity_Score,
    Spatial_FC,
    Spatial_Score,
    Spearman_rho,
    Correlation_Score,
    Mouse_Expression,
    Mouse_Validated,
    Communication_Count,
    Pathway_Count
  )

ml_scaled <- scale(ml_features)

############################################################
# Feature heatmap
############################################################

png(
  "figures/Figure1_Target_Feature_Heatmap.png",
  width = 2400,
  height = 1800,
  res = 300
)

pheatmap(
  ml_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Multi-modal target prioritization features"
)

dev.off()

############################################################
# PCA
############################################################

pca_res <- prcomp(
  ml_scaled,
  center = TRUE,
  scale. = FALSE
)

pca_df <- data.frame(
  Gene = rownames(ml_scaled),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Weighted_AI_Score = feature_matrix$Weighted_AI_Score
)

png(
  "figures/Figure2_Target_PCA.png",
  width = 2200,
  height = 1600,
  res = 300
)

ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    label = Gene,
    color = Weighted_AI_Score
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 4) +
  theme_classic() +
  labs(
    title = "PCA of lupus nephritis target candidates",
    x = "PC1",
    y = "PC2",
    color = "Weighted AI score"
  )

dev.off()

list.files("figures")

############################################################
# Hierarchical clustering
############################################################

dist_targets <- dist(ml_scaled)

hc_targets <- hclust(
  dist_targets,
  method = "ward.D2"
)

png(
  "figures/Figure3_Target_Hierarchical_Clustering.png",
  width = 2200,
  height = 1600,
  res = 300
)

plot(
  hc_targets,
  main = "Hierarchical clustering of candidate targets",
  xlab = "",
  sub = ""
)

dev.off()


write.csv(
  pca_df,
  "results/Target_PCA_Coordinates.csv",
  row.names = FALSE
)


list.files("figures")
list.files("results")