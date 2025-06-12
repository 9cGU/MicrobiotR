# MicroBiotR: A Comprehensive R Toolkit for SOM-Based Microbiota Flow Cytometry and Downstream Machine Learning

## Overview

MicroBiotR is a comprehensive R pipeline developed for analysis between microbiota And host phenotype such as HC vs disease. The package delivers a robust framework for examining microbiota flowcytometry, aiming to pinpoint difference of bacterial community composition potentially linked to environmental factors and host phenotypes.

## Key Features

MicrobiotR provide an All-in-One framework that supports a wide array of microbiota flow cytometry clustering and analysis, including:
* **MBR_som**: SOM calculation and to obtain count-table.
* **MBR_stat**: performing statistical analyses to identify differentially abundant features between groups.
* **MBR_circle**: visualize distributions of features within imported data by presenting them in a circular form.
* **MBR_violin**: visualize abundance or expression levels for selected features.
* **MBR_beta**: visualize and assess beta diversity between groups (PCOA plot)
* **MBR_fs**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_heatmap**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_mantel**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_ml**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_conf**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_reclustering**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_process**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_prepare**: A toolkit for extracting relevant biomarkers from your microbiome data.
* **MBR_plot**: A toolkit for extracting relevant biomarkers from your microbiome data.

## Installation

You can install the released version of `MicroBiotR` from [GitHub](https://github.com/9cGU/MicroBiotR) with the following R commands:

```R
# install.packages("devtools")
devtools::install_github("9cGU/MicroBiotR", force = TRUE)
packageVersion("MicroBiotR")
