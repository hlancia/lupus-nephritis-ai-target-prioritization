# Scripts

This directory contains the complete computational workflow used to develop and evaluate an interpretable multi-evidence framework for therapeutic target prioritization in lupus nephritis.

The analysis integrates transcriptomic, spatial, cross-species, and cell–cell communication evidence into an expert-informed weighted scoring model, followed by unsupervised analyses, interpretability assessment, and sensitivity analysis.

---

## Workflow

```text
01_build_feature_matrix.R
        │
        ▼
Integrated multi-evidence feature matrix
        │
        ▼
02_unsupervised_learning.R
        │
        ▼
Heatmap • PCA • Hierarchical clustering
        │
        ▼
03_interpretable_ranking.R
        │
        ▼
Weighted target prioritization
Feature contributions
Evidence-domain contributions
        │
        ▼
04_project_summary_outputs.R
        │
        ▼
Final project summary tables
Manuscript-ready outputs
Publication-quality figures
        │
        ▼
05_weight_sensitivity_analysis.R
        │
        ▼
Robustness analysis
Rank stability
Weight perturbation simulations
```

---

## Script Descriptions

### 01_build_feature_matrix.R

Builds the integrated therapeutic target feature matrix by combining evidence from multiple independent analyses.

Major tasks include:

- Importing candidate therapeutic targets
- Integrating human single-cell RNA-seq features
- Incorporating human spatial transcriptomic evidence
- Adding cross-species validation results
- Integrating CellChat communication metrics
- Standardizing quantitative features
- Applying expert-defined feature weights
- Calculating the weighted integrated AI target score

Outputs:

- Final_AI_Target_Feature_Matrix.csv
- AI_Feature_Weights.csv
- Top_AI_Prioritized_Targets.csv

---

### 02_unsupervised_learning.R

Performs exploratory analyses of the integrated feature matrix.

Major tasks include:

- Feature heatmap generation
- Principal component analysis (PCA)
- Hierarchical clustering
- Visualization of target similarity

Outputs:

- Target_PCA_Coordinates.csv
- Target_PCA_Explained_Variance.csv
- Target_Hierarchical_Clustering.csv

Figures:

- Figure1_Target_Feature_Heatmap.png
- Figure2_Target_PCA.png
- Figure3_Target_Hierarchical_Clustering.png

---

### 03_interpretable_ranking.R

Generates the interpretable weighted prioritization model.

Major tasks include:

- Reconstruction of weighted AI scores
- Validation of score calculations
- Feature-level contribution analysis
- Evidence-domain contribution analysis
- Visualization of weighted target rankings

Outputs:

- Target_Feature_Contributions.csv
- Target_Evidence_Domain_Contributions.csv

Figures:

- Figure4_Weighted_Target_Ranking.png
- Figure5_Expert_Defined_Feature_Weights.png
- Figure6_Target_Feature_Contributions.png

---

### 04_project_summary_outputs.R

Creates publication-ready summary outputs from the validated prioritization results.

Major tasks include:

- Final therapeutic target summary
- Higher-priority candidate identification
- Evidence-domain summary tables
- Manuscript-ready target tables
- Project-level summary outputs

Outputs:

- Final_Target_Summary.csv
- High_Priority_Targets.csv
- Project3_Key_Findings.csv
- Top4_Targets_Manuscript_Table.csv

Figure:

- Figure7_Evidence_Domain_Contributions.png

---

### 05_weight_sensitivity_analysis.R

Evaluates the robustness of the prioritization framework by assessing the impact of uncertainty in expert-defined feature weights.

Major tasks include:

- Random perturbation of feature weights (±20%)
- 1,000 independent simulations
- Target rank recalculation
- Rank stability assessment
- Pairwise outranking probability analysis
- Robustness interpretation

Outputs:

- Weight_Sensitivity_All_Simulations.csv
- Weight_Sensitivity_Summary.csv
- Weight_Sensitivity_Interpretation.csv
- Weight_Sensitivity_Pairwise_Probability.csv

Figures:

- Figure8_Weight_Sensitivity_Rank_Distributions.png
- Figure9_Weight_Sensitivity_TopRank_Frequencies.png
- Figure10_Weight_Sensitivity_Pairwise_Probability.png

---

## Reproducibility

Run the scripts sequentially:

```text
01_build_feature_matrix.R

↓

02_unsupervised_learning.R

↓

03_interpretable_ranking.R

↓

04_project_summary_outputs.R

↓

05_weight_sensitivity_analysis.R
```

All intermediate tables, publication-quality figures, and summary outputs are generated automatically.

---

## Computational Framework

This repository implements an interpretable computational framework for therapeutic target prioritization by integrating four complementary sources of biological evidence:

- Human single-cell transcriptomics
- Human spatial transcriptomics
- Cross-species validation
- Cell–cell communication analysis

The framework combines expert-informed feature weighting with transparent score decomposition and robustness analysis, enabling reproducible prioritization of therapeutic targets while preserving interpretability.
