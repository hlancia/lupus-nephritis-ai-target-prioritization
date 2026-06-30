############################################################
# 01_build_feature_matrix.R
#
# Project:
# AI-driven therapeutic target prioritization in lupus nephritis
#
# Purpose:
# Integrate outputs from target discovery, spatial validation,
# mouse validation, and CellChat communication analysis into
# a machine-learning-ready feature matrix.
############################################################

library(dplyr)
library(readr)
library(tidyr)


############################################################
# Load Project 1 target-discovery features
############################################################

human_targets <- read.csv(
  "data/project1/Human_Target_Table.csv"
)

spatial_validation <- read.csv(
  "data/project1/Human_Spatial_Validation_GSE263909.csv"
)

spatial_correlation <- read.csv(
  "data/project1/Human_Spatial_Macrophage_Correlation.csv"
)

mouse_expression <- read.csv(
  "data/project1/Mouse_Target_Expression_Cluster9.csv"
)

final_project1 <- read.csv(
  "data/project1/Final_Target_Ranking.csv"
)

human_targets

############################################################
# Load Project 2 cell-cell communication features
############################################################

macro_comm <- read.csv(
  "data/project2/Inflammatory_Macrophage_Communication_Table.csv"
)

macro_pathways <- read.csv(
  "data/project2/Macrophage_Pathway_Summary.csv"
)

pathway_strength <- read.csv(
  "data/project2/Pathway_Communication_Strength.csv"
)

dim(macro_comm)
head(macro_comm)

############################################################
# Define candidate targets
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

feature_matrix


############################################################
# Human single-cell features
############################################################

human_features <- human_targets %>%
  filter(Gene %in% candidate_genes) %>%
  select(
    Gene,
    avg_log2FC,
    pct.1,
    pct.2
  )

feature_matrix <- feature_matrix %>%
  left_join(
    human_features,
    by = "Gene"
  )

feature_matrix

############################################################
# Add spatial transcriptomic features
############################################################

spatial_features <- spatial_validation %>%
  rename(
    Spatial_FC = Spatial_FoldChange,
    Spatial_P = Spatial_Pvalue
  )

feature_matrix <- feature_matrix %>%
  left_join(
    spatial_features,
    by = "Gene"
  )

feature_matrix

############################################################
# Add macrophage correlation features
############################################################

correlation_features <- spatial_correlation %>%
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

feature_matrix


############################################################
# Add mouse validation features
############################################################

mouse_features <- mouse_expression %>%
  rename(
    Mouse_Gene = Gene,
    Mouse_Expression = Mouse_Macrophage_Cluster9_Expression
  )

human_to_mouse <- data.frame(
  Gene = c("C5AR1","CSF1R","LILRB2","PILRA","CLEC7A","TLR4","P2RX7","C3AR1","CD300E","SIGLEC1"),
  Mouse_Gene = c("C5ar1","Csf1r","Lilrb4a","Pilra","Clec7a","Tlr4","P2rx7","C3ar1","Cd300e","Siglec1")
)

feature_matrix <- feature_matrix %>%
  left_join(human_to_mouse, by = "Gene") %>%
  left_join(mouse_features, by = "Mouse_Gene") %>%
  mutate(
    Mouse_Validated = ifelse(Mouse_Expression > 0.25, 1, 0)
  )

feature_matrix

write.csv(
  feature_matrix,
  "results/Target_Feature_Matrix_v1.csv",
  row.names = FALSE
)

############################################################
# Add CellChat communication features
############################################################

communication_features <- macro_comm %>%
  dplyr::filter(
    ligand %in% candidate_genes | receptor %in% candidate_genes
  ) %>%
  dplyr::mutate(
    Gene = ifelse(
      ligand %in% candidate_genes,
      ligand,
      receptor
    )
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    Communication_Count = n(),
    Pathway_Count = dplyr::n_distinct(pathway_name),
    Mean_Communication_Prob = mean(prob),
    Max_Communication_Prob = max(prob),
    .groups = "drop"
  )

feature_matrix <- feature_matrix %>%
  dplyr::left_join(
    communication_features,
    by = "Gene"
  ) %>%
  dplyr::mutate(
    Communication_Count = ifelse(is.na(Communication_Count), 0, Communication_Count),
    Pathway_Count = ifelse(is.na(Pathway_Count), 0, Pathway_Count),
    Mean_Communication_Prob = ifelse(is.na(Mean_Communication_Prob), 0, Mean_Communication_Prob),
    Max_Communication_Prob = ifelse(is.na(Max_Communication_Prob), 0, Max_Communication_Prob)
  )

feature_matrix


write.csv(
  feature_matrix,
  "results/Target_Feature_Matrix_v1.csv",
  row.names = FALSE
)


############################################################
# Build interpretable AI target score
############################################################

feature_matrix <- feature_matrix %>%
  mutate(
    Spatial_Score = -log10(Spatial_P),
    Correlation_Score = -log10(Correlation_Pvalue),
    Specificity_Score = pct.1 - pct.2
  )

feature_matrix$AI_Target_Score <-
  scale(feature_matrix$avg_log2FC)[,1] +
  scale(feature_matrix$Specificity_Score)[,1] +
  scale(feature_matrix$Spatial_FC)[,1] +
  scale(feature_matrix$Spatial_Score)[,1] +
  scale(feature_matrix$Spearman_rho)[,1] +
  scale(feature_matrix$Correlation_Score)[,1] +
  scale(feature_matrix$Mouse_Expression)[,1] +
  scale(feature_matrix$Mouse_Validated)[,1] +
  scale(feature_matrix$Communication_Count)[,1] +
  scale(feature_matrix$Pathway_Count)[,1]

feature_matrix <- feature_matrix %>%
  arrange(desc(AI_Target_Score))

write.csv(
  feature_matrix,
  "results/Target_Feature_Matrix_v1.csv",
  row.names = FALSE
)

feature_matrix

############################################################
# Weighted AI target score
############################################################

feature_matrix$Weighted_AI_Score <-
  1.0 * scale(feature_matrix$avg_log2FC)[,1] +
  1.5 * scale(feature_matrix$Specificity_Score)[,1] +
  2.0 * scale(feature_matrix$Spatial_FC)[,1] +
  2.0 * scale(feature_matrix$Spatial_Score)[,1] +
  1.5 * scale(feature_matrix$Spearman_rho)[,1] +
  1.5 * scale(feature_matrix$Correlation_Score)[,1] +
  2.0 * scale(feature_matrix$Mouse_Validated)[,1] +
  1.0 * scale(feature_matrix$Communication_Count)[,1] +
  1.0 * scale(feature_matrix$Pathway_Count)[,1]

feature_matrix <- feature_matrix %>%
  arrange(desc(Weighted_AI_Score))

feature_matrix %>%
  select(
    Gene,
    AI_Target_Score,
    Weighted_AI_Score
  )

write.csv(
  feature_matrix,
  "results/Final_AI_Target_Feature_Matrix.csv",
  row.names = FALSE
)


############################################################
# Final output check
############################################################

feature_matrix %>%
  select(
    Gene,
    AI_Target_Score,
    Weighted_AI_Score
  )

list.files("results")