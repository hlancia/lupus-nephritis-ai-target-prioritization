############################################################
# 04_project_summary_outputs.R
#
# Project:
# Interpretable multi-evidence therapeutic target
# prioritization in lupus nephritis
#
# Purpose:
# Generate final project-level summary outputs from the
# validated results produced by Scripts 01–03.
#
# This script:
#   1. Creates a concise final target summary table
#   2. Defines the higher-priority candidate group
#   3. Summarizes evidence-domain contributions
#   4. Generates a domain-level evidence figure
#   5. Produces README- and manuscript-ready output tables
#
# Important:
# This script does not recalculate the weighted score.
# It uses the validated scoring results generated previously.
#
# Inputs:
#   results/Final_AI_Target_Feature_Matrix.csv
#   results/Target_Evidence_Domain_Contributions.csv
#   results/Top_AI_Prioritized_Targets.csv
#
# Outputs:
#   results/Final_Target_Summary.csv
#   results/High_Priority_Targets.csv
#   results/Target_Domain_Contribution_Matrix.csv
#   results/Project3_Key_Findings.csv
#
# Figure:
#   figures/Figure7_Evidence_Domain_Contributions.png
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
# Load validated outputs
############################################################

feature_matrix <- read.csv(
  "results/Final_AI_Target_Feature_Matrix.csv",
  check.names = FALSE
)

domain_contributions <- read.csv(
  "results/Target_Evidence_Domain_Contributions.csv",
  check.names = FALSE
)

ranked_targets <- read.csv(
  "results/Top_AI_Prioritized_Targets.csv",
  check.names = FALSE
)


############################################################
# Validate required columns
############################################################

required_feature_columns <- c(
  "Gene",
  "Weighted_AI_Rank",
  "Weighted_AI_Score",
  "avg_log2FC",
  "Specificity_Score",
  "Spatial_FC",
  "Spatial_Score",
  "Spearman_rho",
  "Correlation_Score",
  "Mouse_Ortholog_Status",
  "Mouse_Validated",
  "Communication_Count",
  "Pathway_Count"
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


required_domain_columns <- c(
  "Gene",
  "Evidence_Domain",
  "Domain_Contribution",
  "Weighted_AI_Rank",
  "Weighted_AI_Score"
)

missing_domain_columns <- setdiff(
  required_domain_columns,
  colnames(domain_contributions)
)

if (length(missing_domain_columns) > 0) {
  
  stop(
    paste(
      "The evidence-domain table is missing required columns:",
      paste(
        missing_domain_columns,
        collapse = ", "
      )
    )
  )
}


############################################################
# Confirm ranking consistency across files
############################################################

ranking_check <- feature_matrix %>%
  
  select(
    Gene,
    Weighted_AI_Rank,
    Weighted_AI_Score
  ) %>%
  
  inner_join(
    ranked_targets %>%
      select(
        Gene,
        Weighted_AI_Rank,
        Weighted_AI_Score
      ),
    by = "Gene",
    suffix = c(
      "_FeatureMatrix",
      "_RankedTable"
    )
  ) %>%
  
  mutate(
    Rank_Difference =
      Weighted_AI_Rank_FeatureMatrix -
      Weighted_AI_Rank_RankedTable,
    
    Score_Difference =
      Weighted_AI_Score_FeatureMatrix -
      Weighted_AI_Score_RankedTable
  )


if (
  any(
    ranking_check$Rank_Difference != 0
  ) ||
  any(
    abs(
      ranking_check$Score_Difference
    ) > 1e-8
  )
) {
  
  stop(
    paste(
      "The target ranking is inconsistent across",
      "the validated result files."
    )
  )
}


cat(
  "\nRanking consistency validated across input files.\n"
)


############################################################
# Define priority groups
#
# The top four targets form a clearly separated higher-
# priority group in the current ten-candidate comparison.
#
# This is a descriptive grouping based on the observed ranking
# and score distribution, not a clinically validated threshold.
############################################################

feature_matrix <- feature_matrix %>%
  
  mutate(
    Priority_Group = case_when(
      
      Weighted_AI_Rank <= 4 ~
        "Higher-priority candidate group",
      
      TRUE ~
        "Remaining candidate group"
    )
  )


############################################################
# Create concise evidence labels
############################################################

feature_matrix <- feature_matrix %>%
  
  mutate(
    
    Human_scRNA_Evidence = case_when(
      
      avg_log2FC >= 3 &
        Specificity_Score >= 0.5 ~
        "Strong",
      
      avg_log2FC >= 1.5 &
        Specificity_Score >= 0.25 ~
        "Moderate",
      
      TRUE ~
        "Limited"
    ),
    
    
    Spatial_Evidence = case_when(
      
      Spatial_FC >= 2 &
        Spatial_Score >= 2 &
        Spearman_rho >= 0.4 ~
        "Strong",
      
      Spatial_FC > 1 &
        (
          Spatial_Score >= 1 |
            Spearman_rho >= 0.25
        ) ~
        "Moderate",
      
      TRUE ~
        "Limited"
    ),
    
    
    Cross_Species_Evidence = case_when(
      
      Mouse_Validated == 1 ~
        "Supported",
      
      Mouse_Ortholog_Status == "Family proxy" ~
        "Not directly validated",
      
      TRUE ~
        "Not supported"
    ),
    
    
    Communication_Evidence = case_when(
      
      Communication_Count >= 5 ~
        "Strong",
      
      Communication_Count > 0 ~
        "Limited",
      
      TRUE ~
        "Not detected"
    )
  )


############################################################
# Create final target summary table
############################################################

final_target_summary <- feature_matrix %>%
  
  arrange(
    Weighted_AI_Rank
  ) %>%
  
  transmute(
    
    Rank =
      Weighted_AI_Rank,
    
    Target =
      Gene,
    
    Weighted_Integrated_Score =
      round(
        Weighted_AI_Score,
        3
      ),
    
    Priority_Group =
      Priority_Group,
    
    Human_scRNA_Evidence =
      Human_scRNA_Evidence,
    
    Human_Spatial_Evidence =
      Spatial_Evidence,
    
    Cross_Species_Evidence =
      Cross_Species_Evidence,
    
    CellChat_Evidence =
      Communication_Evidence,
    
    Mouse_Ortholog_Status =
      Mouse_Ortholog_Status,
    
    Communication_Count =
      Communication_Count,
    
    Pathway_Count =
      Pathway_Count
  )


write.csv(
  final_target_summary,
  "results/Final_Target_Summary.csv",
  row.names = FALSE
)


############################################################
# Export higher-priority targets
############################################################

high_priority_targets <- final_target_summary %>%
  
  filter(
    Priority_Group ==
      "Higher-priority candidate group"
  )


write.csv(
  high_priority_targets,
  "results/High_Priority_Targets.csv",
  row.names = FALSE
)


############################################################
# Build target-by-domain contribution matrix
############################################################

domain_contribution_wide <- domain_contributions %>%
  
  select(
    Gene,
    Evidence_Domain,
    Domain_Contribution
  ) %>%
  
  pivot_wider(
    names_from = Evidence_Domain,
    values_from = Domain_Contribution,
    values_fill = 0
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
  domain_contribution_wide,
  "results/Target_Domain_Contribution_Matrix.csv",
  row.names = FALSE
)


############################################################
# Prepare matrix for domain-contribution heatmap
############################################################

expected_domain_order <- c(
  "Human_scRNA",
  "Spatial",
  "Cross_species",
  "Communication"
)


missing_domains <- setdiff(
  expected_domain_order,
  colnames(domain_contribution_wide)
)

if (length(missing_domains) > 0) {
  
  stop(
    paste(
      "Missing expected evidence domains:",
      paste(
        missing_domains,
        collapse = ", "
      )
    )
  )
}


domain_matrix <- domain_contribution_wide %>%
  
  select(
    Gene,
    all_of(
      expected_domain_order
    )
  ) %>%
  
  column_to_rownames(
    "Gene"
  ) %>%
  
  as.matrix()


############################################################
# Preserve final rank order
############################################################

ranked_gene_order <- feature_matrix %>%
  
  arrange(
    Weighted_AI_Rank
  ) %>%
  
  pull(
    Gene
  )


domain_matrix <- domain_matrix[
  ranked_gene_order,
  expected_domain_order,
  drop = FALSE
]


colnames(
  domain_matrix
) <- c(
  "Human\nscRNA-seq",
  "Human\nspatial",
  "Cross-\nspecies",
  "Cell\ncommunication"
)


############################################################
# Figure 7: evidence-domain contributions
############################################################

maximum_absolute_domain_value <- max(
  abs(
    domain_matrix
  ),
  na.rm = TRUE
)


domain_breaks <- seq(
  -maximum_absolute_domain_value,
  maximum_absolute_domain_value,
  length.out = 101
)


png(
  "figures/Figure7_Evidence_Domain_Contributions.png",
  width = 2600,
  height = 2200,
  res = 300
)

pheatmap(
  domain_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  breaks = domain_breaks,
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
  fontsize_col = 11,
  angle_col = 0,
  main =
    paste(
      "Evidence-domain contributions to",
      "therapeutic target prioritization"
    )
)

dev.off()


############################################################
# Create concise key-findings table
############################################################

top_target <- final_target_summary %>%
  
  filter(
    Rank == 1
  )


second_target <- final_target_summary %>%
  
  filter(
    Rank == 2
  )


third_target <- final_target_summary %>%
  
  filter(
    Rank == 3
  )


fourth_target <- final_target_summary %>%
  
  filter(
    Rank == 4
  )


project_key_findings <- data.frame(
  
  Finding = c(
    "Lead prioritized target",
    "Higher-priority candidate group",
    "Lead-target score",
    "Targets with direct mouse validation",
    "Targets with detected CellChat evidence",
    "Interpretation of negative scores"
  ),
  
  Result = c(
    top_target$Target,
    
    paste(
      high_priority_targets$Target,
      collapse = ", "
    ),
    
    as.character(
      top_target$Weighted_Integrated_Score
    ),
    
    paste(
      feature_matrix %>%
        filter(
          Mouse_Validated == 1
        ) %>%
        arrange(
          Weighted_AI_Rank
        ) %>%
        pull(
          Gene
        ),
      collapse = ", "
    ),
    
    paste(
      feature_matrix %>%
        filter(
          Communication_Count > 0
        ) %>%
        arrange(
          Weighted_AI_Rank
        ) %>%
        pull(
          Gene
        ),
      collapse = ", "
    ),
    
    paste(
      "Below-average integrated evidence",
      "within the ten-candidate comparison set;",
      "not evidence against biological relevance."
    )
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  project_key_findings,
  "results/Project3_Key_Findings.csv",
  row.names = FALSE
)


############################################################
# Create a manuscript-ready top-four table
############################################################

top_four_manuscript_table <- feature_matrix %>%
  
  filter(
    Weighted_AI_Rank <= 4
  ) %>%
  
  arrange(
    Weighted_AI_Rank
  ) %>%
  
  transmute(
    
    Rank =
      Weighted_AI_Rank,
    
    Target =
      Gene,
    
    Weighted_Score =
      round(
        Weighted_AI_Score,
        2
      ),
    
    Human_Macrophage_Enrichment =
      Human_scRNA_Evidence,
    
    Human_Spatial_Support =
      Spatial_Evidence,
    
    Mouse_Validation =
      Cross_Species_Evidence,
    
    CellChat_Support =
      Communication_Evidence,
    
    Interpretation = case_when(
      
      Gene == "C5AR1" ~
        paste(
          "Strong convergent support from human macrophage",
          "enrichment, spatial evidence, macrophage correlation,",
          "and direct mouse validation."
        ),
      
      Gene == "LILRB2" ~
        paste(
          "Strong human and spatial evidence; no strict",
          "one-to-one mouse validation."
        ),
      
      Gene == "PILRA" ~
        paste(
          "Multi-domain support including direct mouse",
          "validation and detected CellChat interactions."
        ),
      
      Gene == "CSF1R" ~
        paste(
          "Human macrophage enrichment, spatial support,",
          "and direct mouse validation."
        ),
      
      TRUE ~
        "Higher-priority multi-evidence candidate."
    )
  )


write.csv(
  top_four_manuscript_table,
  "results/Top4_Targets_Manuscript_Table.csv",
  row.names = FALSE
)


############################################################
# Console summary
############################################################

cat(
  "\nScript 04 completed successfully.\n"
)


cat(
  "\nFinal target summary:\n"
)

print(
  final_target_summary
)


cat(
  "\nHigher-priority candidate group:\n"
)

print(
  high_priority_targets
)


cat(
  "\nProject key findings:\n"
)

print(
  project_key_findings
)


cat(
  "\nFiles generated:\n"
)

print(
  c(
    "figures/Figure7_Evidence_Domain_Contributions.png",
    "results/Final_Target_Summary.csv",
    "results/High_Priority_Targets.csv",
    "results/Target_Domain_Contribution_Matrix.csv",
    "results/Project3_Key_Findings.csv",
    "results/Top4_Targets_Manuscript_Table.csv"
  )
)