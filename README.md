# HCC_B4GALT3_612-AYR_Code

Code repository for the bioinformatics analyses accompanying the
manuscript:

**Enzyme-responsive 612-AYR self-assembled nanoparticles inhibit B4GALT3
activity, remodel ITGB1 glycan structures, and suppress hepatocellular
carcinoma growth and metastasis**

This repository contains all R and Perl scripts used for the
transcriptomic, survival, enrichment, and immune-infiltration analyses
based on TCGA-LIHC. These scripts support the identification of B4GALT3
as a prognostic glycosyltransferase in hepatocellular carcinoma (HCC),
reveal its regulatory effects on cell-cycle progression and ITGB1
glycosylation, and provide the computational foundation for the
development of enzyme-responsive 612-AYR nanoparticles.

## 1. Overview

Hepatocellular carcinoma (HCC) is a highly lethal malignancy
characterized by frequent recurrence and treatment resistance. Through
differential expression, survival screening, Cox regression, ROC
analysis and functional annotation using TCGA-LIHC data, we identified
**β-1,4-Galactosyltransferase 3 (B4GALT3)** as a prognostic biomarker.
B4GALT3 expression correlates with poor survival and enrichment of
cell-cycle pathways.

Functional assays demonstrate that B4GALT3 enhances tumor proliferation,
migration and invasion by stabilizing and glycosylating **Integrin
Beta-1 (ITGB1)**. An **enzyme-responsive 612-AYR nanoparticle delivery
system** was designed to target this axis, showing tumor-specific
accumulation, prolonged retention, and potent anti-HCC activity.

This repository provides the full set of scripts used to reproduce the
TCGA-based analyses in the study.

## 2. Contents of the Repository

All files are located in the root directory and are listed below in
approximate workflow order:

-   **1.diff.R** --- Differential expression analysis (limma)
-   **2.mergeTime.R** --- Merge expression data with survival/clinical
    data
-   **3.SurvivalFilter.R** --- KM survival screening
-   **4.geneSur.R** --- Survival curves
-   **5.indepFilter.R** --- Clinical independence preprocessing
-   **6.cox.R** --- Univariate Cox regression
-   **7.ROC.R** --- Time-dependent ROC analysis
-   **8.clincialCorFilter.R** --- Clinical data cleanup
-   **9.cliCor.R** --- Expression vs. clinical subgroup analysis
-   **10.preGSEA.pl** --- GSEA preprocessing
-   **11.multiGSEA.R** --- GSEA visualization
-   **12.CIBERSORT.R** --- Immune cell fraction estimation
-   **13.vioplot.R** --- Immune infiltration violin plots
-   **14.immunePlot.R** --- Immune heatmaps & composition plots

## 3. Data Sources

-   TCGA-LIHC RNA-Seq & clinical data (UCSC Xena)
-   Glycosyltransferase gene lists from GeneCards
-   CIBERSORT (LM22 matrix required separately)

## 4. Software Requirements

-   R ≥ 4.0\
-   Key packages: limma, survival, survminer, pheatmap, ggplot2, ggpubr,
    corrplot, vioplot, etc.
-   Perl for GSEA preprocessing

## 5. Typical Workflow

``` r
source("1.diff.R")
source("2.mergeTime.R")
source("3.SurvivalFilter.R")
source("6.cox.R")
source("7.ROC.R")
source("8.clincialCorFilter.R")
source("9.cliCor.R")
perl 10.preGSEA.pl
source("11.multiGSEA.R")
source("12.CIBERSORT.R")
source("13.vioplot.R")
source("14.immunePlot.R")
```

## 6. Citation

If you use this code, please cite:

**Enzyme-responsive 612-AYR nanoparticles targeting B4GALT3-mediated
ITGB1 glycosylation for precision treatment of hepatocellular
carcinoma.**

## 7. Contact

For questions regarding the code or workflow, please contact the
corresponding author.
