############################################################
# 05_weight_sensitivity_analysis.R
#
# Project:
# Interpretable multi-evidence therapeutic target
# prioritization in lupus nephritis
#
# Purpose:
# Evaluate whether the target ranking remains stable when the
# expert-defined feature weights are varied within a plausible
# range.
#
# Sensitivity strategy:
#   - Independently perturb each feature weight by ±20%
#   - Recalculate target scores
#   - Recalculate target ranks
#   - Repeat for 1,000 simulations
#
# Important:
# This is a local sensitivity analysis around the predefined
# weighting model. It does not test every possible scoring
# architecture or eliminate the subjectivity of expert-defined
# feature weights.
#
# Inputs:
#   results/Final_AI_Target_Feature_Matrix.csv
#   results/AI_Feature_Weights.csv
#
# Outputs:
#   results/Weight_Sensitivity_All_Simulations.csv
#   results/Weight_Sensitivity_Summary.csv
#   results/Weight_Sensitivity_Perturbed_Weights.csv
#   results/Weight_Sensitivity_Pairwise_Probability.csv
#
# Figures:
#   figures/Figure8_Weight_Sensitivity_Rank_Distributions.png
#   figures/Figure9_Weight_Sensitivity_TopRank_Frequencies.png
#   figures/Figure10_Weight_Sensitivity_Pairwise_Probability.png
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
# Analysis parameters
############################################################

set.seed(20260804)

n_simulations <- 1000

perturbation_fraction <- 0.20

top_group_threshold <- 4


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
# Load validated feature matrix and default weights
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

required_feature_columns <- c(
  "Gene",
  "Weighted_AI_Rank",
  "Weighted_AI_Score"
)

missing_feature_columns <- setdiff(
  required_feature_columns,
  colnames(feature_matrix)
)

if (length(missing_feature_columns) > 0) {
  
  stop(
    paste(
      "The final feature matrix is missing required columns:",
      paste(
        missing_feature_columns,
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
    "Duplicate feature names were detected in the weight table."
  )
}


############################################################
# Identify standardized scoring columns
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
      "Missing standardized scoring columns:",
      paste(
        missing_standardized_columns,
        collapse = ", "
      )
    )
  )
}


############################################################
# Construct target-by-feature standardized matrix
############################################################

standardized_feature_matrix <- feature_matrix %>%
  
  select(
    Gene,
    all_of(
      feature_weights$Standardized_Column
    )
  ) %>%
  
  column_to_rownames(
    "Gene"
  ) %>%
  
  as.matrix()


storage.mode(
  standardized_feature_matrix
) <- "numeric"


if (anyNA(standardized_feature_matrix)) {
  
  stop(
    paste(
      "Missing values were detected in the standardized",
      "feature matrix."
    )
  )
}


############################################################
# Preserve target and feature order
############################################################

target_names <- rownames(
  standardized_feature_matrix
)

feature_names <- feature_weights$Feature

default_weights <- feature_weights$Weight

names(default_weights) <- feature_names


############################################################
# Verify the default score before simulations
############################################################

default_reconstructed_score <-
  as.numeric(
    standardized_feature_matrix %*%
      default_weights
  )


default_score_check <- data.frame(
  
  Gene =
    target_names,
  
  Reconstructed_Default_Score =
    default_reconstructed_score,
  
  stringsAsFactors = FALSE
  
) %>%
  
  left_join(
    feature_matrix %>%
      select(
        Gene,
        Weighted_AI_Score,
        Weighted_AI_Rank
      ),
    by = "Gene"
  ) %>%
  
  mutate(
    Score_Difference =
      Reconstructed_Default_Score -
      Weighted_AI_Score
  )


maximum_default_difference <- max(
  abs(
    default_score_check$Score_Difference
  ),
  na.rm = TRUE
)


if (maximum_default_difference > 1e-8) {
  
  stop(
    paste(
      "The default score could not be reconstructed.",
      "Maximum difference:",
      maximum_default_difference
    )
  )
}


cat(
  "\nDefault weighted scores successfully reconstructed.\n"
)

cat(
  "Maximum numerical difference:",
  maximum_default_difference,
  "\n"
)


############################################################
# Generate perturbed feature weights
#
# Each feature weight is sampled independently from:
#
# default weight × [0.80, 1.20]
############################################################

weight_multiplier_matrix <- matrix(
  
  runif(
    n =
      n_simulations *
      length(default_weights),
    
    min =
      1 -
      perturbation_fraction,
    
    max =
      1 +
      perturbation_fraction
  ),
  
  nrow =
    n_simulations,
  
  ncol =
    length(default_weights),
  
  byrow =
    FALSE
)


perturbed_weight_matrix <- sweep(
  
  weight_multiplier_matrix,
  
  MARGIN = 2,
  
  STATS = default_weights,
  
  FUN = "*"
)


colnames(
  perturbed_weight_matrix
) <- feature_names


############################################################
# Save perturbed weights
############################################################

perturbed_weights_table <-
  
  as.data.frame(
    perturbed_weight_matrix
  ) %>%
  
  rownames_to_column(
    "Simulation"
  ) %>%
  
  mutate(
    Simulation =
      as.integer(
        Simulation
      )
  )


write.csv(
  perturbed_weights_table,
  "results/Weight_Sensitivity_Perturbed_Weights.csv",
  row.names = FALSE
)


############################################################
# Calculate simulated target scores
#
# standardized_feature_matrix:
#   targets × features
#
# perturbed_weight_matrix:
#   simulations × features
#
# Result:
#   targets × simulations
############################################################

simulated_score_matrix <-
  
  standardized_feature_matrix %*%
  
  t(
    perturbed_weight_matrix
  )


rownames(
  simulated_score_matrix
) <- target_names


colnames(
  simulated_score_matrix
) <- paste0(
  "Simulation_",
  seq_len(
    n_simulations
  )
)


############################################################
# Convert scores to ranks
#
# Rank 1 represents the highest simulated score.
############################################################

simulated_rank_matrix <- apply(
  
  simulated_score_matrix,
  
  MARGIN = 2,
  
  FUN = function(current_scores) {
    
    rank(
      -current_scores,
      ties.method = "min"
    )
  }
)


rownames(
  simulated_rank_matrix
) <- target_names


colnames(
  simulated_rank_matrix
) <- paste0(
  "Simulation_",
  seq_len(
    n_simulations
  )
)


############################################################
# Convert simulation results to long format
############################################################

simulation_scores_long <-
  
  as.data.frame(
    simulated_score_matrix
  ) %>%
  
  rownames_to_column(
    "Gene"
  ) %>%
  
  pivot_longer(
    cols = -Gene,
    names_to = "Simulation",
    values_to = "Simulated_Score"
  )


simulation_ranks_long <-
  
  as.data.frame(
    simulated_rank_matrix
  ) %>%
  
  rownames_to_column(
    "Gene"
  ) %>%
  
  pivot_longer(
    cols = -Gene,
    names_to = "Simulation",
    values_to = "Simulated_Rank"
  )


simulation_results <- simulation_scores_long %>%
  
  left_join(
    simulation_ranks_long,
    by = c(
      "Gene",
      "Simulation"
    )
  ) %>%
  
  mutate(
    Simulation =
      as.integer(
        sub(
          "Simulation_",
          "",
          Simulation
        )
      )
  ) %>%
  
  left_join(
    feature_matrix %>%
      select(
        Gene,
        Default_Rank =
          Weighted_AI_Rank,
        Default_Score =
          Weighted_AI_Score
      ),
    by = "Gene"
  ) %>%
  
  arrange(
    Simulation,
    Simulated_Rank,
    Gene
  )


write.csv(
  simulation_results,
  "results/Weight_Sensitivity_All_Simulations.csv",
  row.names = FALSE
)


############################################################
# Summarize ranking robustness
############################################################

sensitivity_summary <- simulation_results %>%
  
  group_by(
    Gene
  ) %>%
  
  summarise(
    
    Default_Rank =
      first(
        Default_Rank
      ),
    
    Default_Score =
      first(
        Default_Score
      ),
    
    Mean_Simulated_Rank =
      mean(
        Simulated_Rank
      ),
    
    Median_Simulated_Rank =
      median(
        Simulated_Rank
      ),
    
    Rank_SD =
      sd(
        Simulated_Rank
      ),
    
    Best_Rank =
      min(
        Simulated_Rank
      ),
    
    Worst_Rank =
      max(
        Simulated_Rank
      ),
    
    Mean_Simulated_Score =
      mean(
        Simulated_Score
      ),
    
    Simulated_Score_SD =
      sd(
        Simulated_Score
      ),
    
    Top1_Frequency =
      mean(
        Simulated_Rank == 1
      ),
    
    Top3_Frequency =
      mean(
        Simulated_Rank <= 3
      ),
    
    Top4_Frequency =
      mean(
        Simulated_Rank <=
          top_group_threshold
      ),
    
    Default_Rank_Frequency =
      mean(
        Simulated_Rank ==
          Default_Rank
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    Top1_Percent =
      100 *
      Top1_Frequency,
    
    Top3_Percent =
      100 *
      Top3_Frequency,
    
    Top4_Percent =
      100 *
      Top4_Frequency,
    
    Default_Rank_Percent =
      100 *
      Default_Rank_Frequency
  ) %>%
  
  arrange(
    Default_Rank
  )


write.csv(
  sensitivity_summary,
  "results/Weight_Sensitivity_Summary.csv",
  row.names = FALSE
)


############################################################
# Pairwise outranking probabilities
#
# Cell [i, j] represents the percentage of simulations in
# which target i ranked above target j.
############################################################

pairwise_probability_matrix <- matrix(
  
  NA_real_,
  
  nrow =
    length(target_names),
  
  ncol =
    length(target_names),
  
  dimnames = list(
    target_names,
    target_names
  )
)


for (
  target_i in target_names
) {
  
  for (
    target_j in target_names
  ) {
    
    if (
      target_i == target_j
    ) {
      
      pairwise_probability_matrix[
        target_i,
        target_j
      ] <- 0.5
      
    } else {
      
      pairwise_probability_matrix[
        target_i,
        target_j
      ] <-
        
        mean(
          simulated_rank_matrix[
            target_i,
          ] <
            simulated_rank_matrix[
              target_j,
            ]
        )
    }
  }
}


############################################################
# Reorder pairwise matrix by default rank
############################################################

default_gene_order <- feature_matrix %>%
  
  arrange(
    Weighted_AI_Rank
  ) %>%
  
  pull(
    Gene
  )


pairwise_probability_matrix <-
  
  pairwise_probability_matrix[
    default_gene_order,
    default_gene_order,
    drop = FALSE
  ]


pairwise_probability_table <-
  
  as.data.frame(
    pairwise_probability_matrix
  ) %>%
  
  rownames_to_column(
    "Target"
  )


write.csv(
  pairwise_probability_table,
  "results/Weight_Sensitivity_Pairwise_Probability.csv",
  row.names = FALSE
)


############################################################
# Figure 8: simulated rank distributions
############################################################

rank_plot_data <- simulation_results %>%
  
  left_join(
    sensitivity_summary %>%
      select(
        Gene,
        Median_Simulated_Rank
      ),
    by = "Gene"
  ) %>%
  
  mutate(
    Gene =
      factor(
        Gene,
        levels =
          rev(
            default_gene_order
          )
      )
  )


rank_distribution_plot <- ggplot(
  rank_plot_data,
  aes(
    x = Gene,
    y = Simulated_Rank
  )
) +
  
  geom_violin(
    fill = "#D9E6F2",
    color = "#4C78A8",
    scale = "width",
    trim = TRUE
  ) +
  
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white"
  ) +
  
  geom_point(
    data =
      sensitivity_summary %>%
      mutate(
        Gene =
          factor(
            Gene,
            levels =
              rev(
                default_gene_order
              )
          )
      ),
    
    aes(
      x = Gene,
      y = Default_Rank
    ),
    
    inherit.aes = FALSE,
    shape = 21,
    size = 3.5,
    fill = "#B2182B"
  ) +
  
  coord_flip() +
  
  scale_y_reverse(
    breaks = seq_len(
      length(target_names)
    )
  ) +
  
  labs(
    title =
      "Target-rank stability under feature-weight perturbation",
    
    subtitle =
      paste0(
        n_simulations,
        " simulations with independent ±",
        perturbation_fraction * 100,
        "% variation in each feature weight"
      ),
    
    x =
      "Candidate target",
    
    y =
      "Simulated rank",
    
    caption =
      paste(
        "Red points indicate the default target rank.",
        "Rank 1 represents the highest integrated score."
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
    
    plot.caption =
      element_text(
        hjust = 0
      )
  )


ggsave(
  "figures/Figure8_Weight_Sensitivity_Rank_Distributions.png",
  rank_distribution_plot,
  width = 10,
  height = 7.5,
  dpi = 300
)


############################################################
# Figure 9: top-rank frequencies
############################################################

frequency_plot_data <- sensitivity_summary %>%
  
  select(
    Gene,
    Default_Rank,
    Top1_Percent,
    Top3_Percent,
    Top4_Percent
  ) %>%
  
  pivot_longer(
    cols = c(
      Top1_Percent,
      Top3_Percent,
      Top4_Percent
    ),
    names_to = "Rank_Category",
    values_to = "Frequency_Percent"
  ) %>%
  
  mutate(
    
    Rank_Category = recode(
      Rank_Category,
      
      Top1_Percent =
        "Ranked #1",
      
      Top3_Percent =
        "Ranked in top 3",
      
      Top4_Percent =
        "Ranked in top 4"
    ),
    
    Rank_Category =
      factor(
        Rank_Category,
        levels = c(
          "Ranked #1",
          "Ranked in top 3",
          "Ranked in top 4"
        )
      ),
    
    Gene =
      factor(
        Gene,
        levels =
          rev(
            default_gene_order
          )
      )
  )


top_frequency_plot <- ggplot(
  frequency_plot_data,
  aes(
    x = Gene,
    y = Frequency_Percent,
    fill = Rank_Category
  )
) +
  
  geom_col(
    position =
      position_dodge(
        width = 0.8
      ),
    
    width = 0.72
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    
    breaks = seq(
      0,
      100,
      by = 20
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
      "Frequency of high-priority rankings across weight simulations",
    
    subtitle =
      "Robust candidates retain high ranks despite perturbation of the weighting scheme",
    
    x =
      "Candidate target",
    
    y =
      "Simulation frequency (%)",
    
    fill =
      "Rank category"
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
  "figures/Figure9_Weight_Sensitivity_TopRank_Frequencies.png",
  top_frequency_plot,
  width = 10,
  height = 7.5,
  dpi = 300
)


############################################################
# Figure 9: top-rank frequencies
############################################################

frequency_plot_data <- sensitivity_summary %>%

  select(
    Gene,
    Default_Rank,
    Top1_Percent,
    Top3_Percent,
    Top4_Percent
  ) %>%

  pivot_longer(
    cols = c(
      Top1_Percent,
      Top3_Percent,
      Top4_Percent
    ),
    names_to = "Rank_Category",
    values_to = "Frequency_Percent"
  ) %>%

  mutate(

    Rank_Category = recode(
      Rank_Category,

      Top1_Percent =
        "Ranked #1",

      Top3_Percent =
        "Ranked in top 3",

      Top4_Percent =
        "Ranked in top 4"
    ),

    Rank_Category =
      factor(
        Rank_Category,
        levels = c(
          "Ranked #1",
          "Ranked in top 3",
          "Ranked in top 4"
        )
      ),

    Gene =
      factor(
        Gene,
        levels =
          rev(
            default_gene_order
          )
      )
  )


top_frequency_plot <- ggplot(
  frequency_plot_data,
  aes(
    x = Gene,
    y = Frequency_Percent,
    fill = Rank_Category
  )
) +

  geom_col(
    position =
      position_dodge(
        width = 0.8
      ),

    width = 0.72
  ) +

  coord_flip() +

  scale_y_continuous(
    limits = c(
      0,
      100
    ),

    breaks = seq(
      0,
      100,
      by = 20
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
      "Frequency of high-priority rankings across weight simulations",

    subtitle =
      "Robust candidates retain high ranks despite perturbation of the weighting scheme",

    x =
      "Candidate target",

    y =
      "Simulation frequency (%)",

    fill =
      "Rank category"
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
  "figures/Figure9_Weight_Sensitivity_TopRank_Frequencies.png",
  top_frequency_plot,
  width = 10,
  height = 7.5,
  dpi = 300
)


############################################################
# Figure 10: pairwise outranking probability
############################################################

pairwise_display_matrix <-
  100 *
  pairwise_probability_matrix


png(
  "figures/Figure10_Weight_Sensitivity_Pairwise_Probability.png",
  width = 2800,
  height = 2400,
  res = 300
)

pheatmap(
  pairwise_display_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(
    100
  ),
  breaks = seq(
    0,
    100,
    length.out = 101
  ),
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  border_color = "white",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  main =
    paste(
      "Probability that the row target ranks",
      "above the column target (%)"
    )
)

dev.off()


############################################################
# Create concise robustness interpretation
############################################################

robustness_interpretation <- sensitivity_summary %>%
  
  mutate(
    
    Rank_Stability = case_when(
      
      Rank_SD <= 0.10 ~
        "Highly stable",
      
      Rank_SD <= 0.50 ~
        "Stable",
      
      Rank_SD <= 1.00 ~
        "Moderately stable",
      
      TRUE ~
        "Sensitive to weight variation"
    ),
    
    Priority_Group_Stability = case_when(
      
      Top4_Percent >= 95 ~
        "Stable top-four member",
      
      Top4_Percent >= 50 ~
        "Occasional top-four member",
      
      TRUE ~
        "Outside top-four group"
    )
  ) %>%
  
  select(
    Default_Rank,
    Gene,
    Default_Score,
    Median_Simulated_Rank,
    Rank_SD,
    Best_Rank,
    Worst_Rank,
    Top1_Percent,
    Top3_Percent,
    Top4_Percent,
    Default_Rank_Percent,
    Rank_Stability,
    Priority_Group_Stability
  ) %>%
  
  arrange(
    Default_Rank
  )


############################################################
# Save robustness interpretation
############################################################

write.csv(
  robustness_interpretation,
  "results/Weight_Sensitivity_Interpretation.csv",
  row.names = FALSE
)


############################################################
# Console summary
############################################################

cat(
  "\nScript 05 completed successfully.\n"
)

cat(
  "\nSensitivity-analysis parameters:\n"
)

cat(
  "Number of simulations:",
  n_simulations,
  "\n"
)

cat(
  "Feature-weight perturbation: ±",
  perturbation_fraction * 100,
  "%\n",
  sep = ""
)

cat(
  "\nRanking robustness summary:\n"
)

print(
  robustness_interpretation
)

cat(
  "\nFiles generated:\n"
)

print(
  c(
    "figures/Figure8_Weight_Sensitivity_Rank_Distributions.png",
    "figures/Figure9_Weight_Sensitivity_TopRank_Frequencies.png",
    "figures/Figure10_Weight_Sensitivity_Pairwise_Probability.png",
    "results/Weight_Sensitivity_All_Simulations.csv",
    "results/Weight_Sensitivity_Summary.csv",
    "results/Weight_Sensitivity_Perturbed_Weights.csv",
    "results/Weight_Sensitivity_Pairwise_Probability.csv",
    "results/Weight_Sensitivity_Interpretation.csv"
  )
)