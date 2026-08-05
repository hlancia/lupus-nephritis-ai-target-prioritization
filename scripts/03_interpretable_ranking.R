############################################################
# 03_interpretable_ranking.R
#
# Project:
# Interpretable multi-evidence therapeutic target
# prioritization in lupus nephritis
#
# Purpose:
# Generate an interpretable description of the weighted target
# prioritization model.
#
# This script:
#   1. Validates the weighted target scores
#   2. Calculates feature-level score contributions
#   3. Generates the final target-ranking figure
#   4. Visualizes the expert-defined weighting scheme
#   5. Generates a target-by-feature contribution heatmap
#
# Important:
# The feature weights used in this project are expert-defined.
# They are not learned from a supervised machine-learning model.
#
# Inputs:
#   results/Final_AI_Target_Feature_Matrix.csv
#   results/AI_Feature_Weights.csv
#
# Outputs:
#   results/Target_Feature_Contributions.csv
#   results/Target_Feature_Contributions_Wide.csv
#   results/Top_AI_Prioritized_Targets.csv
#
# Figures:
#   figures/Figure4_Weighted_Target_Ranking.png
#   figures/Figure5_Expert_Defined_Feature_Weights.png
#   figures/Figure6_Target_Feature_Contributions.png
############################################################


############################################################
# Load packages
############################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
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
# Load finalized feature matrix and weighting scheme
############################################################

feature_matrix <- read.csv(
  "results/Final_AI_Target_Feature_Matrix.csv",
  check.names = FALSE
)

feature_weights <- read.csv(
  "results/AI_Feature_Weights.csv",
  check.names = FALSE
)


############################################################
# Validate required columns
############################################################

required_matrix_columns <- c(
  "Gene",
  "Weighted_AI_Score",
  "Weighted_AI_Rank"
)

missing_matrix_columns <- setdiff(
  required_matrix_columns,
  colnames(feature_matrix)
)

if (length(missing_matrix_columns) > 0) {
  
  stop(
    paste(
      "The final feature matrix is missing required columns:",
      paste(
        missing_matrix_columns,
        collapse = ", "
      )
    )
  )
}


required_weight_columns <- c(
  "Feature",
  "Evidence_Domain",
  "Weight"
)

missing_weight_columns <- setdiff(
  required_weight_columns,
  colnames(feature_weights)
)

if (length(missing_weight_columns) > 0) {
  
  stop(
    paste(
      "The feature-weight table is missing required columns:",
      paste(
        missing_weight_columns,
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


if (anyDuplicated(feature_weights$Feature) > 0) {
  
  stop(
    "Duplicate feature names were detected in the feature-weight table."
  )
}


############################################################
# Verify that standardized features are available
############################################################

feature_weights <- feature_weights %>%
  
  mutate(
    Standardized_Column =
      paste0(
        "z_",
        Feature
      )
  )


missing_standardized_columns <- setdiff(
  feature_weights$Standardized_Column,
  colnames(feature_matrix)
)

if (length(missing_standardized_columns) > 0) {
  
  stop(
    paste(
      "Missing standardized feature columns:",
      paste(
        missing_standardized_columns,
        collapse = ", "
      )
    )
  )
}


############################################################
# Build target-by-feature contribution table
#
# Feature contribution:
#
# standardized feature value × expert-defined feature weight
############################################################

contribution_list <- lapply(
  
  seq_len(
    nrow(feature_weights)
  ),
  
  function(i) {
    
    current_feature <-
      feature_weights$Feature[i]
    
    current_domain <-
      feature_weights$Evidence_Domain[i]
    
    current_weight <-
      feature_weights$Weight[i]
    
    standardized_column <-
      feature_weights$Standardized_Column[i]
    
    data.frame(
      Gene =
        feature_matrix$Gene,
      
      Feature =
        current_feature,
      
      Evidence_Domain =
        current_domain,
      
      Weight =
        current_weight,
      
      Standardized_Value =
        feature_matrix[[standardized_column]],
      
      Weighted_Contribution =
        feature_matrix[[standardized_column]] *
        current_weight,
      
      stringsAsFactors = FALSE
    )
  }
)


feature_contributions <- bind_rows(
  contribution_list
)


############################################################
# Add target rank and final score
############################################################

feature_contributions <- feature_contributions %>%
  
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
    Weighted_AI_Rank,
    Feature
  )


############################################################
# Verify score reconstruction
############################################################

reconstructed_scores <- feature_contributions %>%
  
  group_by(Gene) %>%
  
  summarise(
    Reconstructed_Weighted_Score =
      sum(
        Weighted_Contribution,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    feature_matrix %>%
      select(
        Gene,
        Weighted_AI_Score
      ),
    by = "Gene"
  ) %>%
  
  mutate(
    Score_Difference =
      Reconstructed_Weighted_Score -
      Weighted_AI_Score
  )


maximum_score_difference <- max(
  abs(
    reconstructed_scores$Score_Difference
  ),
  na.rm = TRUE
)


if (maximum_score_difference > 1e-8) {
  
  stop(
    paste(
      "The reconstructed weighted scores do not match",
      "the scores in the final feature matrix.",
      "Maximum difference:",
      maximum_score_difference
    )
  )
}


cat(
  "\nWeighted-score reconstruction validated.\n"
)

cat(
  "Maximum numerical difference:",
  maximum_score_difference,
  "\n"
)


############################################################
# Save long-format contribution table
############################################################

write.csv(
  feature_contributions,
  "results/Target_Feature_Contributions.csv",
  row.names = FALSE
)


############################################################
# Create and save wide-format contribution matrix
############################################################

feature_contributions_wide <- feature_contributions %>%
  
  select(
    Gene,
    Feature,
    Weighted_Contribution
  ) %>%
  
  pivot_wider(
    names_from = Feature,
    values_from = Weighted_Contribution
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
  feature_contributions_wide,
  "results/Target_Feature_Contributions_Wide.csv",
  row.names = FALSE
)


############################################################
# Create final ranked target table
############################################################

ranked_targets <- feature_matrix %>%
  
  arrange(
    Weighted_AI_Rank
  ) %>%
  
  select(
    Weighted_AI_Rank,
    Gene,
    Weighted_AI_Score,
    AI_Target_Rank,
    AI_Target_Score,
    Mouse_Ortholog_Status,
    Mouse_Validated,
    Communication_Count,
    Pathway_Count
  )


write.csv(
  ranked_targets,
  "results/Top_AI_Prioritized_Targets.csv",
  row.names = FALSE
)


############################################################
# Figure 4: final weighted target ranking
#
# Negative scores indicate evidence profiles below the average
# of the ten-candidate set. They do not represent biological
# evidence against a target.
############################################################

ranking_plot_data <- ranked_targets %>%
  
  mutate(
    Gene =
      factor(
        Gene,
        levels = rev(Gene)
      ),
    
    Priority_Group = case_when(
      Weighted_AI_Rank <= 4 ~
        "Higher-priority candidate group",
      
      TRUE ~
        "Remaining candidate group"
    )
  )


ranking_plot <- ggplot(
  ranking_plot_data,
  aes(
    x = Gene,
    y = Weighted_AI_Score,
    fill = Priority_Group
  )
) +
  
  geom_col(
    width = 0.72
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label =
        sprintf(
          "%.2f",
          Weighted_AI_Score
        )
    ),
    
    hjust = ifelse(
      ranking_plot_data$Weighted_AI_Score >= 0,
      -0.12,
      1.12
    ),
    
    size = 3.7
  ) +
  
  coord_flip(
    clip = "off"
  ) +
  
  scale_fill_manual(
    values = c(
      "Higher-priority candidate group" =
        "#3B82F6",
      
      "Remaining candidate group" =
        "#B8C2CC"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.16,
        0.15
      )
    )
  ) +
  
  labs(
    title =
      "Interpretable multi-evidence ranking of therapeutic target candidates",
    
    subtitle =
      paste(
        "Scores are expert-weighted sums of standardized",
        "evidence features"
      ),
    
    x =
      "Candidate target",
    
    y =
      "Weighted integrated evidence score",
    
    fill =
      NULL,
    
    caption =
      paste(
        "Negative values indicate below-average evidence",
        "within the ten-candidate comparison set."
      )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "bottom",
    
    plot.caption =
      element_text(
        hjust = 0
      ),
    
    plot.margin =
      margin(
        10,
        35,
        10,
        10
      )
  )


ggsave(
  "figures/Figure4_Weighted_Target_Ranking.png",
  ranking_plot,
  width = 10,
  height = 7,
  dpi = 300
)


############################################################
# Figure 5: expert-defined feature weights
############################################################

weight_plot_data <- feature_weights %>%
  
  mutate(
    Feature_Label = recode(
      Feature,
      
      avg_log2FC =
        "Human macrophage log2 fold change",
      
      Specificity_Score =
        "Human macrophage specificity",
      
      Spatial_FC =
        "Spatial fold change",
      
      Spatial_Score =
        "Spatial statistical support",
      
      Spearman_rho =
        "Spatial macrophage correlation",
      
      Correlation_Score =
        "Correlation statistical support",
      
      Mouse_Validated =
        "Mouse cross-species validation",
      
      Communication_Count =
        "CellChat interaction count",
      
      Pathway_Count =
        "CellChat pathway count"
    ),
    
    Feature_Label =
      factor(
        Feature_Label,
        levels =
          rev(
            Feature_Label[
              order(Weight)
            ]
          )
      ),
    
    Evidence_Domain =
      recode(
        Evidence_Domain,
        
        Human_scRNA =
          "Human scRNA-seq",
        
        Spatial =
          "Human spatial",
        
        Cross_species =
          "Cross-species",
        
        Communication =
          "Cell communication"
      )
  )


weight_plot <- ggplot(
  weight_plot_data,
  aes(
    x = Feature_Label,
    y = Weight,
    fill = Evidence_Domain
  )
) +
  
  geom_col(
    width = 0.72
  ) +
  
  geom_text(
    aes(
      label = Weight
    ),
    
    hjust = -0.25,
    size = 3.8
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(
      0,
      max(
        weight_plot_data$Weight
      ) +
        0.5
    ),
    
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  labs(
    title =
      "Expert-defined feature weighting scheme",
    
    subtitle =
      paste(
        "Weights reflect the predefined contribution",
        "of each evidence feature"
      ),
    
    x =
      "Evidence feature",
    
    y =
      "Assigned feature weight",
    
    fill =
      "Evidence domain"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "bottom"
  )


ggsave(
  "figures/Figure5_Expert_Defined_Feature_Weights.png",
  weight_plot,
  width = 10,
  height = 7,
  dpi = 300
)


############################################################
# Figure 6: target-by-feature contribution heatmap
############################################################

contribution_matrix <- feature_contributions %>%
  
  select(
    Gene,
    Feature,
    Weighted_Contribution
  ) %>%
  
  pivot_wider(
    names_from = Feature,
    values_from = Weighted_Contribution
  ) %>%
  
  column_to_rownames(
    "Gene"
  ) %>%
  
  as.matrix()


############################################################
# Reorder rows according to final target ranking
############################################################

ranked_gene_order <- ranked_targets$Gene

contribution_matrix <-
  contribution_matrix[
    ranked_gene_order,
    feature_weights$Feature,
    drop = FALSE
  ]


############################################################
# Improve feature labels
############################################################

colnames(
  contribution_matrix
) <- c(
  "Human\nlog2FC",
  "Human\nspecificity",
  "Spatial\nfold change",
  "Spatial\nsignificance",
  "Spatial\ncorrelation",
  "Correlation\nsignificance",
  "Mouse\nvalidation",
  "Communication\ncount",
  "Pathway\ncount"
)


############################################################
# Symmetric color scale centered on zero
############################################################

maximum_absolute_contribution <- max(
  abs(
    contribution_matrix
  ),
  na.rm = TRUE
)


heatmap_breaks <- seq(
  -maximum_absolute_contribution,
  maximum_absolute_contribution,
  length.out = 101
)


png(
  "figures/Figure6_Target_Feature_Contributions.png",
  width = 3000,
  height = 2200,
  res = 300
)

pheatmap(
  contribution_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  breaks = heatmap_breaks,
  color = colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(
    100
  ),
  border_color = "white",
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 10,
  angle_col = 45,
  main =
    paste(
      "Feature-level contributions to the",
      "weighted target-prioritization score"
    )
)

dev.off()


############################################################
# Summarize evidence-domain contributions
############################################################

domain_contributions <- feature_contributions %>%
  
  group_by(
    Gene,
    Evidence_Domain
  ) %>%
  
  summarise(
    Domain_Contribution =
      sum(
        Weighted_Contribution,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    ranked_targets %>%
      select(
        Gene,
        Weighted_AI_Rank,
        Weighted_AI_Score
      ),
    by = "Gene"
  ) %>%
  
  arrange(
    Weighted_AI_Rank,
    Evidence_Domain
  )


write.csv(
  domain_contributions,
  "results/Target_Evidence_Domain_Contributions.csv",
  row.names = FALSE
)


############################################################
# Console summary
############################################################

cat(
  "\nScript 03 completed successfully.\n"
)


cat(
  "\nFinal target ranking:\n"
)

print(
  ranked_targets %>%
    select(
      Weighted_AI_Rank,
      Gene,
      Weighted_AI_Score
    )
)


cat(
  "\nEvidence-domain contributions for the top four targets:\n"
)

print(
  domain_contributions %>%
    filter(
      Weighted_AI_Rank <= 4
    )
)


cat(
  "\nFiles generated:\n"
)

print(
  c(
    "figures/Figure4_Weighted_Target_Ranking.png",
    "figures/Figure5_Expert_Defined_Feature_Weights.png",
    "figures/Figure6_Target_Feature_Contributions.png",
    "results/Target_Feature_Contributions.csv",
    "results/Target_Feature_Contributions_Wide.csv",
    "results/Target_Evidence_Domain_Contributions.csv",
    "results/Top_AI_Prioritized_Targets.csv"
  )
)