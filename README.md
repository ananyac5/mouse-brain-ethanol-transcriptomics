# Transcriptomic Analysis of Mouse Brain Regions Following Ethanol Exposure

### Differential Gene Expression and Machine Learning | University of Texas at Austin

## Project Overview

This undergraduate computational biology research project investigates how acute and chronic ethanol exposure affects gene expression across five brain regions in mice.

Using RNA sequencing (RNA-seq) gene-count data, our research team performed differential gene expression analysis, examined region-specific transcriptional responses, and developed a K-Nearest Neighbors (KNN) classifier to distinguish acute and chronic ethanol exposure conditions.

The project was conducted as part of the Big Data in Biology (BDiB) research program at the University of Texas at Austin.

## Research Objectives

The primary research question was:

**How do acute and chronic ethanol exposure influence gene expression across different regions of the mouse brain?**

The study focused on five brain regions:

1. Cortex
2. Dorsal hippocampus
3. Ventral hippocampus
4. Dorsal striatum
5. Ventral striatum

Our objectives included:

- Identifying differentially expressed genes (DEGs) associated with ethanol exposure.
- Comparing transcriptional differences across brain regions.
- Examining differences between acute and chronic exposure.
- Investigating biologically relevant genes and pathways.
- Developing a machine learning classifier to distinguish exposure conditions using gene-expression profiles.

## Technologies and Tools

**Programming**
- Python
- Jupyter Notebook

**Data Analysis**
- pandas
- NumPy
- SciPy
- Statistical hypothesis testing
- Counts per million (CPM) normalization
- Log2 transformation

**Machine Learning**
- scikit-learn
- K-Nearest Neighbors (KNN)
- StandardScaler
- Train-test splitting

**Data Visualization**
- Matplotlib
- seaborn
- Hierarchical clustering
- Gene expression heatmaps
- Volcano plots

## Methodology

### 1. RNA-seq Data Preprocessing

Raw gene-count data were imported and processed using Python.

The preprocessing workflow included separating gene metadata from expression counts, organizing samples by treatment condition, and normalizing sequencing counts.

Counts per million (CPM) normalization and log2 transformation were applied to facilitate exploratory gene-expression comparisons.

### 2. Differential Gene Expression Analysis

Differential gene expression was explored by comparing ethanol-treated and control samples under acute and chronic conditions.

The analysis calculated:

- Expression differences between experimental groups
- Log2-scale expression differences
- Gene-wise statistical test p-values
- Genes meeting exploratory significance and expression-change thresholds

Additional analyses examined ethanol-associated expression differences across the five brain regions.

Volcano plots were generated to visualize expression changes and statistical significance.

**Note:** The original exploratory workflow used gene-wise t-tests and unadjusted p-values. The reported DEG counts should be interpreted cautiously because the analysis did not implement a multiple-testing correction.

### 3. Gene Expression Heatmap

Hierarchical clustering and heatmap visualization were used to examine expression patterns among selected genes across experimental samples.

The visualization explored similarities and differences in expression profiles associated with ethanol exposure and brain region.

### 4. K-Nearest Neighbors Classification

A KNN classifier was developed using scikit-learn to distinguish samples associated with acute and chronic ethanol exposure.

The original workflow included:

- Gene-expression feature selection
- An 80/20 training and testing split
- StandardScaler feature normalization
- KNN classification
- Evaluation across multiple neighborhood sizes

The original analysis reported a highest test accuracy of **85.71% at K = 1**.

This is an exploratory result rather than a validated estimate of generalization performance. Because feature selection occurred before the train-test split and multiple K values were compared using the test data, the workflow requires further validation to address potential data leakage and model-selection bias.

### 5. Biological Pathway Exploration

Selected genes were examined in the context of biological pathways potentially associated with ethanol exposure and neurological function.

The project presentation discussed candidate genes including *Nr4a2*, *Tbr1*, and *Plin4*, as well as pathways related to neuronal development and signaling.

These interpretations are exploratory and do not establish causal mechanisms.

## Key Findings

### Differential Gene Expression Across Brain Regions

The original project presentation reported the following numbers of differentially expressed genes:

| Brain Region | Reported DEGs |
|---|---:|
| Cortex | 0 |
| Dorsal Hippocampus | 16 |
| Ventral Hippocampus | 11 |
| Dorsal Striatum | 118 |
| Ventral Striatum | 74 |

The dorsal and ventral striatum showed the largest numbers of identified DEGs under the project's exploratory thresholds.

### Acute and Chronic Conditions

The presentation reported:

- Acute ethanol versus control: 3 DEGs
- Chronic ethanol versus control: 2 DEGs

### Machine Learning

The original KNN classifier achieved a reported maximum test accuracy of **85.71%**, using one nearest neighbor.

Further evaluation with training-only feature selection, cross-validation, and an independent test set would be needed to establish a more reliable performance estimate.

## Research Visualizations

Selected visualizations are provided in the `images/` directory.

These include:

- Differential gene expression heatmap
- KNN classifier accuracy graph
- Selected brain-region differential expression plots

The complete research presentation is available in the `presentation/` directory.

## Repository Structure

```text
mouse-brain-ethanol-transcriptomics/
├── README.md
├── Code/
│   ├── 01_preprocessing_differential_expression.py
│   └── 02_knn_classifier.py
├── Figures/
└── Presentation/
    └── AUD_in_Mice_Research_Presentation.pdf
```

### Code Descriptions

| Script | Description |
|---|---|
| `01_preprocessing_differential_expression.py` | RNA-seq preprocessing, CPM normalization, differential expression analysis across five brain regions, and visualization. |
| `02_knn_classifier.py` | Gene-expression feature preparation, train-test splitting, feature scaling, KNN classification, and accuracy evaluation. |

Both scripts depend on the original RNA-seq dataset and should be reviewed before being treated as fully reproducible analytical pipelines.

## Research Team

**University of Texas at Austin — Big Data in Biology**

Research team:
- Ananya Chintalapati
- Kaye Chen
- Deepika Chitirala

### Individual Contributions

*Individual research and programming contributions to be documented here.*

This repository represents a collaborative undergraduate research project. All team members should receive appropriate credit.

## Dataset and References

The project analyzed mouse brain RNA-seq data associated with:

Periandri, Erica M., et al. (2026). *Epigenetic and Transcriptomic Impacts of Ethanol Vary by Brain Region and Extent of Exposure*. eNeuro.

DOI: https://doi.org/10.1523/ENEURO.0484-25.2026

The original gene-count dataset is not included in this repository.

## Limitations and Future Improvements

- The study used relatively small numbers of biological replicates per experimental group.
- The exploratory DEG workflow used unadjusted p-values and requires further statistical validation.
- Count-based differential expression methods and false discovery rate correction would improve the analysis.
- The KNN classifier requires a leakage-free feature-selection and model-validation workflow.
- Additional biological validation is needed to assess the significance of identified genes and pathways.

## Acknowledgments

We acknowledge the Big Data in Biology research program at the University of Texas at Austin for providing research guidance and the opportunity to apply computational and statistical methods to biological datasets.

---

**Research Areas:** Computational Biology | Bioinformatics | RNA Sequencing | Machine Learning | Statistical Analysis
