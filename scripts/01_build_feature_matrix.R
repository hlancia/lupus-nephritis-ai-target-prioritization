############################################################
# 01_build_feature_matrix.R
#
# Project:
# Interpretable multi-evidence therapeutic target
# prioritization in lupus nephritis
#
# Purpose:
# Integrate evidence from:
#   1. Human inflammatory macrophage expression
#   2. Human spatial transcriptomics
#   3. Spatial macrophage correlation
#   4. Mouse cross-species validation
#   5. CellChat communication analysis
#
# Outputs:
#   results/Target_Feature_Matrix_v1.csv
#   results/Final_AI_Target_Feature_Matrix.csv
#   results/AI_Feature_Weights.csv
############################################################


############################################################
# Load packages
############################################################

library(dplyr)
library(readr)
library(tidyr)


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
# Helper functions
############################################################

safe_neg_log10 <- function(x, minimum_p = 1e-300) {
  
  x <- as.numeric(x)
  
  x[is.na(x)] <- 1
  
  x <- pmax(
    x,
    minimum_p
  )
  
  -log10(x)
}


safe_zscore <- function(x) {
  
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  
  x[is.na(x)] <- median(
    x,
    na.rm = TRUE
  )
  
  current_sd <- sd(
    x,
    na.rm = TRUE
  )
  
  if (
    is.na(current_sd) ||
    current_sd == 0
  ) {
    return(rep(0, length(x)))
  }
  
  as.numeric(scale(x))
}


############################################################
# Load Project 1 evidence
############################################################

human_targets <- read.csv(
  "C:/Users/hlanc/Documents/AMP_LN_Atlas/results/Human_Target_Table.csv",
  check.names = FALSE
)

spatial_validation <- read.csv(
  "C:/Users/hlanc/Documents/AMP_LN_Atlas/results/Human_Spatial_Validation_GSE263909.csv",
  check.names = FALSE
)

spatial_correlation <- read.csv(
  "C:/Users/hlanc/Documents/AMP_LN_Atlas/results/Human_Spatial_Macrophage_Correlation.csv",
  check.names = FALSE
)

mouse_expression <- read.csv(
  "C:/Users/hlanc/Documents/AMP_LN_Atlas/results/Mouse_Target_Expression_Cluster9.csv",
  check.names = FALSE
)


############################################################
# Load Project 2 evidence
############################################################

macro_comm <- read.csv(
  "C:/Users/hlanc/Documents/LN_AI_Target_Prioritization/data/project2/Inflammatory_Macrophage_Communication_Table.csv",
  check.names = FALSE
)


############################################################
# Define candidate targets
#
# These are preselected therapeutic target candidates carried
# forward from Projects 1 and 2. This project prioritizes these
# candidates; it does not perform genome-wide target discovery.
############################################################

candidate_genes <- c(
  "C5AR1",
  "CSF1R",
  "LILRB2",
  "PILRA",
  "CLEC7A",
  "TLR4",
  "P2RX7",
  "C3AR1",
  "CD300E",
  "SIGLEC1"
)


############################################################
# Initialize feature matrix
############################################################

feature_matrix <- data.frame(
  Gene = candidate_genes,
  stringsAsFactors = FALSE
)


############################################################
# Human single-cell features
############################################################

human_features <- human_targets %>%
  filter(
    Gene %in% candidate_genes
  ) %>%
  select(
    Gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  )

feature_matrix <- feature_matrix %>%
  left_join(
    human_features,
    by = "Gene"
  )


############################################################
# Human spatial transcriptomic features
############################################################

spatial_features <- spatial_validation %>%
  filter(
    Gene %in% candidate_genes
  ) %>%
  transmute(
    Gene,
    Spatial_FC = Spatial_FoldChange,
    Spatial_P = Spatial_Pvalue
  )

feature_matrix <- feature_matrix %>%
  left_join(
    spatial_features,
    by = "Gene"
  )


############################################################
# Spatial macrophage correlation features
############################################################

correlation_features <- spatial_correlation %>%
  filter(
    Gene %in% candidate_genes
  ) %>%
  select(
    Gene,
    Spearman_rho,
    Correlation_Pvalue
  )

feature_matrix <- feature_matrix %>%
  left_join(
    correlation_features,
    by = "Gene"
  )


############################################################
# Mouse cross-species validation
############################################################

mouse_features <- mouse_expression %>%
  transmute(
    Mouse_Gene = Gene,
    Mouse_Expression =
      Mouse_Macrophage_Cluster9_Expression
  )


############################################################
# Human-to-mouse mapping
#
# LILRB2 does not have a strict one-to-one mouse ortholog.
# Lilrb4a is retained as a related family proxy and should
# be interpreted cautiously.
############################################################

human_to_mouse <- data.frame(
  Gene = c(
    "C5AR1",
    "CSF1R",
    "LILRB2",
    "PILRA",
    "CLEC7A",
    "TLR4",
    "P2RX7",
    "C3AR1",
    "CD300E",
    "SIGLEC1"
  ),
  
  Mouse_Gene = c(
    "C5ar1",
    "Csf1r",
    "Lilrb4a",
    "Pilra",
    "Clec7a",
    "Tlr4",
    "P2rx7",
    "C3ar1",
    "Cd300e",
    "Siglec1"
  ),
  
  Mouse_Ortholog_Status = c(
    "Direct ortholog",
    "Direct ortholog",
    "Family proxy",
    "Direct ortholog",
    "Direct ortholog",
    "Direct ortholog",
    "Direct ortholog",
    "Direct ortholog",
    "Direct ortholog",
    "Direct ortholog"
  ),
  
  stringsAsFactors = FALSE
)


feature_matrix <- feature_matrix %>%
  left_join(
    human_to_mouse,
    by = "Gene"
  ) %>%
  left_join(
    mouse_features,
    by = "Mouse_Gene"
  ) %>%
  mutate(
    
    Mouse_Validated = case_when(
      
      Mouse_Ortholog_Status == "Family proxy" ~ 0,
      
      !is.na(Mouse_Expression) &
        Mouse_Expression > 0.25 ~ 1,
      
      TRUE ~ 0
    )
  )


############################################################
# CellChat communication features
#
# A gene is counted when it appears as either a ligand or a
# receptor in an interaction involving inflammatory
# macrophages.
############################################################

communication_features <- macro_comm %>%
  
  filter(
    ligand %in% candidate_genes |
      receptor %in% candidate_genes
  ) %>%
  
  mutate(
    Gene = case_when(
      ligand %in% candidate_genes ~ ligand,
      receptor %in% candidate_genes ~ receptor,
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(Gene)
  ) %>%
  
  group_by(Gene) %>%
  
  summarise(
    Communication_Count = n(),
    
    Pathway_Count =
      n_distinct(pathway_name),
    
    Mean_Communication_Prob =
      mean(prob, na.rm = TRUE),
    
    Max_Communication_Prob =
      max(prob, na.rm = TRUE),
    
    .groups = "drop"
  )


feature_matrix <- feature_matrix %>%
  
  left_join(
    communication_features,
    by = "Gene"
  ) %>%
  
  mutate(
    Communication_Count =
      replace_na(Communication_Count, 0),
    
    Pathway_Count =
      replace_na(Pathway_Count, 0),
    
    Mean_Communication_Prob =
      replace_na(Mean_Communication_Prob, 0),
    
    Max_Communication_Prob =
      replace_na(Max_Communication_Prob, 0)
  )


############################################################
# Derived interpretable features
############################################################

feature_matrix <- feature_matrix %>%
  
  mutate(
    
    Specificity_Score =
      pct.1 - pct.2,
    
    Human_DE_Score =
      safe_neg_log10(p_val_adj),
    
    Spatial_Score =
      safe_neg_log10(Spatial_P),
    
    Correlation_Score =
      safe_neg_log10(Correlation_Pvalue)
  )


############################################################
# Save raw integrated feature matrix
############################################################

write.csv(
  feature_matrix,
  "results/Target_Feature_Matrix_v1.csv",
  row.names = FALSE
)


############################################################
# Define final scoring features and expert-informed weights
#
# Mouse_Expression is retained in the feature matrix for
# exploration, but the final weighted score uses the binary
# Mouse_Validated feature. This avoids allowing differences
# in expression scale between datasets to dominate the score.
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
  
  Evidence_Domain = c(
    "Human_scRNA",
    "Human_scRNA",
    "Spatial",
    "Spatial",
    "Spatial",
    "Spatial",
    "Cross_species",
    "Communication",
    "Communication"
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
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  feature_weights,
  "results/AI_Feature_Weights.csv",
  row.names = FALSE
)


############################################################
# Standardize model features
############################################################

scoring_features <- feature_weights$Feature


for (
  current_feature in scoring_features
) {
  
  standardized_name <- paste0(
    "z_",
    current_feature
  )
  
  feature_matrix[[standardized_name]] <-
    safe_zscore(
      feature_matrix[[current_feature]]
    )
}


############################################################
# Unweighted integrated score
############################################################

standardized_columns <- paste0(
  "z_",
  scoring_features
)


feature_matrix$AI_Target_Score <-
  rowSums(
    feature_matrix[
      standardized_columns
    ],
    na.rm = TRUE
  )


############################################################
# Weighted integrated score
############################################################

weighted_components <- sapply(
  
  seq_len(
    nrow(feature_weights)
  ),
  
  function(i) {
    
    feature_name <-
      feature_weights$Feature[i]
    
    weight_value <-
      feature_weights$Weight[i]
    
    standardized_name <-
      paste0(
        "z_",
        feature_name
      )
    
    feature_matrix[[standardized_name]] * weight_value
  }
)

feature_matrix$Weighted_AI_Score <-
  rowSums(
    weighted_components,
    na.rm = TRUE
  )



############################################################
# Add rank variables
############################################################

feature_matrix <- feature_matrix %>%
  
  mutate(
    
    AI_Target_Rank =
      min_rank(
        desc(AI_Target_Score)
      ),
    
    Weighted_AI_Rank =
      min_rank(
        desc(Weighted_AI_Score)
      )
  ) %>%
  
  arrange(
    Weighted_AI_Rank,
    AI_Target_Rank,
    Gene
  )


############################################################
# Save final matrix used by Scripts 02–04
############################################################

write.csv(
  feature_matrix,
  "results/Final_AI_Target_Feature_Matrix.csv",
  row.names = FALSE
)


############################################################
# Console checks
############################################################

cat(
  "\nFinal candidate ranking:\n"
)

print(
  feature_matrix %>%
    select(
      Weighted_AI_Rank,
      Gene,
      AI_Target_Score,
      Weighted_AI_Score,
      Mouse_Validated,
      Communication_Count,
      Pathway_Count
    )
)


cat(
  "\nFiles written:\n"
)

print(
  c(
    "results/Target_Feature_Matrix_v1.csv",
    "results/Final_AI_Target_Feature_Matrix.csv",
    "results/AI_Feature_Weights.csv"
  )
)