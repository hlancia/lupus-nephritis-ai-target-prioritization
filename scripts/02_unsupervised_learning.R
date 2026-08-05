############################################################
# 02_unsupervised_learning.R
#
# Project:
# Interpretable multi-evidence therapeutic target
# prioritization in lupus nephritis
#
# Purpose:
# Perform exploratory unsupervised analysis of the integrated
# target feature matrix using:
#   1. Feature heatmap
#   2. Principal Component Analysis
#   3. Hierarchical clustering
#
# Important:
# These analyses visualize similarities among ten preselected
# candidate targets. They are exploratory and are not used to
# train or validate the target-ranking model.
#
# Inputs:
#   results/Final_AI_Target_Feature_Matrix.csv
#
# Outputs:
#   results/Target_PCA_Coordinates.csv
#   results/Target_Hierarchical_Clustering.csv
#
# Figures:
#   figures/Figure1_Target_Feature_Heatmap.png
#   figures/Figure2_Target_PCA.png
#   figures/Figure3_Target_Hierarchical_Clustering.png
############################################################


############################################################
# Load packages
############################################################

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)


############################################################
# Create output directories if missing
############################################################

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# Load final integrated feature matrix
############################################################

feature_matrix <- read.csv(
  "results/Final_AI_Target_Feature_Matrix.csv",
  check.names = FALSE
)


############################################################
# Define features used for exploratory analysis
#
# The features below correspond to the evidence variables used
# in the final weighted prioritization model.
############################################################

analysis_features <- c(
  "avg_log2FC",
  "Specificity_Score",
  "Spatial_FC",
  "Spatial_Score",
  "Spearman_rho",
  "Correlation_Score",
  "Mouse_Validated",
  "Communication_Count",
  "Pathway_Count"
)


############################################################
# Validate input columns
############################################################

missing_features <- setdiff(
  analysis_features,
  colnames(feature_matrix)
)

if (length(missing_features) > 0) {
  
  stop(
    paste(
      "Missing required features:",
      paste(
        missing_features,
        collapse = ", "
      )
    )
  )
}


if (anyDuplicated(feature_matrix$Gene) > 0) {
  
  stop(
    "Duplicate target names were detected in the feature matrix."
  )
}


############################################################
# Construct target-by-feature matrix
############################################################

target_matrix <- feature_matrix %>%
  
  select(
    Gene,
    all_of(analysis_features)
  ) %>%
  
  column_to_rownames(
    "Gene"
  ) %>%
  
  as.matrix()


############################################################
# Verify complete numeric matrix
############################################################

storage.mode(target_matrix) <- "numeric"


if (anyNA(target_matrix)) {
  
  stop(
    "Missing values were detected in the exploratory feature matrix."
  )
}


############################################################
# Standardize features
#
# Features are standardized by column so variables measured on
# different numerical scales contribute comparably.
############################################################

scaled_target_matrix <- scale(
  target_matrix
)


if (anyNA(scaled_target_matrix)) {
  
  stop(
    paste(
      "Standardization produced missing values.",
      "Check for features with zero variance."
    )
  )
}


############################################################
# Figure 1: target feature heatmap
############################################################

png(
  "figures/Figure1_Target_Feature_Heatmap.png",
  width = 2800,
  height = 2200,
  res = 300
)

pheatmap(
  scaled_target_matrix,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  border_color = NA,
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 10,
  angle_col = 45,
  main = paste(
    "Integrated evidence profiles of",
    "lupus nephritis target candidates"
  )
)

dev.off()


############################################################
# Principal Component Analysis
############################################################

target_pca <- prcomp(
  scaled_target_matrix,
  center = FALSE,
  scale. = FALSE
)


############################################################
# Calculate explained variance
############################################################

explained_variance <-
  100 *
  target_pca$sdev^2 /
  sum(target_pca$sdev^2)


############################################################
# Create and save PCA coordinate table
############################################################

pca_coordinates <- as.data.frame(
  target_pca$x
) %>%
  
  rownames_to_column(
    "Gene"
  ) %>%
  
  left_join(
    feature_matrix %>%
      select(
        Gene,
        Weighted_AI_Rank,
        Weighted_AI_Score
      ),
    by = "Gene"
  ) %>%
  
  arrange(
    Weighted_AI_Rank
  )


write.csv(
  pca_coordinates,
  "results/Target_PCA_Coordinates.csv",
  row.names = FALSE
)


############################################################
# Figure 2: target PCA
############################################################

pca_plot <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    size = Weighted_AI_Score,
    label = Gene
  )
) +
  
  geom_point(
    alpha = 0.8
  ) +
  
  geom_text_repel(
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.4,
    min.segment.length = 0
  ) +
  
  scale_size_continuous(
    range = c(3, 9)
  ) +
  
  labs(
    title = paste(
      "Exploratory PCA of integrated",
      "target evidence profiles"
    ),
    subtitle = paste(
      "Point size represents the weighted",
      "target-prioritization score"
    ),
    x = paste0(
      "PC1 (",
      round(explained_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(explained_variance[2], 1),
      "%)"
    ),
    size = "Weighted\nscore"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  "figures/Figure2_Target_PCA.png",
  pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)


############################################################
# Hierarchical clustering
############################################################

target_distance <- dist(
  scaled_target_matrix,
  method = "euclidean"
)

target_hclust <- hclust(
  target_distance,
  method = "ward.D2"
)


############################################################
# Save hierarchical cluster order
############################################################

clustering_table <- data.frame(
  Dendrogram_Order =
    seq_along(target_hclust$order),
  
  Gene =
    target_hclust$labels[
      target_hclust$order
    ],
  
  stringsAsFactors = FALSE
) %>%
  
  left_join(
    feature_matrix %>%
      select(
        Gene,
        Weighted_AI_Rank,
        Weighted_AI_Score
      ),
    by = "Gene"
  )


write.csv(
  clustering_table,
  "results/Target_Hierarchical_Clustering.csv",
  row.names = FALSE
)


############################################################
# Figure 3: hierarchical clustering
############################################################

png(
  "figures/Figure3_Target_Hierarchical_Clustering.png",
  width = 2600,
  height = 1800,
  res = 300
)

plot(
  target_hclust,
  main = paste(
    "Hierarchical clustering of integrated",
    "target evidence profiles"
  ),
  xlab = "",
  sub = "",
  ylab = "Euclidean distance",
  cex = 0.9,
  hang = -1
)

dev.off()


############################################################
# Save PCA explained variance
############################################################

pca_variance_table <- data.frame(
  Principal_Component =
    paste0(
      "PC",
      seq_along(explained_variance)
    ),
  
  Explained_Variance_Percent =
    explained_variance
)


write.csv(
  pca_variance_table,
  "results/Target_PCA_Explained_Variance.csv",
  row.names = FALSE
)


############################################################
# Console summary
############################################################

cat(
  "\nScript 02 completed successfully.\n"
)

cat(
  "\nPCA variance explained:\n"
)

print(
  head(
    pca_variance_table,
    5
  )
)


cat(
  "\nFiles generated:\n"
)

print(
  c(
    "figures/Figure1_Target_Feature_Heatmap.png",
    "figures/Figure2_Target_PCA.png",
    "figures/Figure3_Target_Hierarchical_Clustering.png",
    "results/Target_PCA_Coordinates.csv",
    "results/Target_PCA_Explained_Variance.csv",
    "results/Target_Hierarchical_Clustering.csv"
  )
)