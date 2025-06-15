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

```markdown
# install.packages("devtools")
devtools::install_github("9cGU/MicroBiotR", force = TRUE)
packageVersion("MicroBiotR")
```

## SOM analysis
```markdown
# load matadata
meta<-read.delim('meta.txt',header = T, row.names = 1)
# som analysis
MBR_som(fl_data_ig)
```

## Statistics
```markdown
# test type should be one of 'wilcoxon', 'kruskal', 'anova'
# correction could be none or fdr
MBR_stat(data = count_table, group_col = 'Group', meta_data = meta, 
       test_type = 'wilcoxon', out_path = './', 
       correction = 'none', cutoff = 0.008)
```

## Analysis
```markdown
MBR_circle(data = significant_data, group_col = 'Group', meta_data = meta, width = 8, 
         height = 8, out_path = './')

MBR_violin(data = significant_data, meta_data = meta, pvalue_data = pvalue_data, 
         cluster = 709, group_col = 'Group', colors = c('#FF7F00', '#4DAF4A'), out_path = './')

MBR_beta(original, out_path = './', test = 't.test', 
       meta_data = meta, group_name = 'Group')

MBR_heatmap(data = MBR_selected_features[,1:6], 
          cohonen_information = cohonen_information, 
          out_path = './', scale = 'row', cluster_cols = F, 
          cluster_rows = T, display_numbers = F)

MBR_mantel(
  data = MBR_selected_features,
  meta_data = meta,
  clinical_cols = c("Clinical1", "Clinical2"),
  demographic_cols = c("Demographic1", "Demographic2"),
  spec_select_names = list(A = "Clinical", B = "Demographic"),
  out_path = './',  
  width = 6,
  height = 6
)

```



