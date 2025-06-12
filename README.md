# MicroBiotR: A Comprehensive R Toolkit for SOM-Based Microbiota Flow Cytometry and Downstream Machine Learning

## Overview

MicroBiotR is a comprehensive R pipeline developed for analysis between microbiota And host phenotype such as HC vs disease. The package delivers a robust framework for examining microbiota flowcytometry, aiming to pinpoint difference of bacterial community composition potentially linked to environmental factors and host phenotypes.

## Key Features

MicrobiotR provide an All-in-One framework that supports a wide array of microbiota flow cytometry clustering and analysis, including:
* **MBR_som**: SOM calculation and to obtain count-table.
* **MBR_stat**: performe statistical analyses to identify differentially abundant features between groups.
* **MBR_circle**: visualize distributions of features within imported data by presenting them in a circular form.
* **MBR_violin**: visualize abundance or expression levels for selected features.
* **MBR_beta**: visualize and assess beta diversity between groups (PCOA plot).
* **MBR_fs**: identifying and selecting a subset of the most predictive features.
* **MBR_heatmap**: visualize patterns of channels used in flowcytometry using a color-coded matrix.
* **MBR_mantel**: quantify and visualize the relationship between features in data and clinical or demographic metadata.
* **MBR_ml**: apply machine learning algorithms and visualize using ROC curve.
* **MBR_conf**: visualize confusion matrix.
* **MBR_reclustering**: recluster the clustering results obtained from an initial Self-Organizing Map (SOM) analysis.
* **MBR_read**: load fcs file.
* **MBR_process**: prepare data for downstream analytical steps.
* **MBR_prepare**: prepare data for downstream analytical steps.
* **MBR_plot**: visualize dot plot of selected features.

## Installation

You can install the released version of `MicroBiotR` from [GitHub](https://github.com/9cGU/MicroBiotR) with the following R commands:

```R
# install.packages("devtools")
devtools::install_github("9cGU/MicroBiotR", force = TRUE)
packageVersion("MicroBiotR")

## SOM analysis

You can install the released version of `MicroBiotR` from [GitHub](https://github.com/9cGU/MicroBiotR) with the following R commands:

```R
# install.packages("devtools")
devtools::install_github("9cGU/MicroBiotR", force = TRUE)
packageVersion("MicroBiotR")
