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

You can firstly install these dependencies prior to installing `MicroBiotR`

```markdown
# install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)
packageVersion("linkET")

install.packages("viridis")

devtools::install_github("Bioconductor/Biobase", force = TRUE)
install.packages("BiocManager", repos = "https://cran.R-project.org")
BiocManager::install("flowCore")

```


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
# test type should be one of 'wilcoxon', 'kruskal', 'anova' ,'ttest'
# correction could be one of 'none', 'fdr', 'bonferroni' ,'BH'
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

![MBR Circle Plot](images/circle.png)

## Machine learning
```markdown
MBR_ml(
  data = MBR_selected_features,        # your feature table (samples as rows, taxa as columns)
  meta_data = meta,          # metadata dataframe (samples in same order as `data`)
  group_name = 'Group',             # column name in metadata containing group labels
  out_path = './',          # directory to save output PDF and group file
  reference_level = 'IBD',            # reference group (e.g., "A" or "Control")
  width = 6,                        # width of the output plot PDF
  height = 6,                       # height of the output plot PDF
  method = 'repeatedcv',           # resampling method (e.g., repeated cross-validation)
  number = 2,                       # number of folds
  repeats = 2                       # number of repeats
)

MBR_conf(
  data = MBR_selected_features,
  meta_data = meta,
  group_name = 'Group',
  reference_level = 'IBD',
  out_path = './'
)
```

## Reclustering
```markdown
MBR_reclustering(data = cohonen_information, num_clusters = 1000)
```

## flowcytometry dotplot
```markdown
rawdata_path <- "/MappedFCS"

fcs_files <- MBR_read(rawdata_path)

selected_rows <- c('V1269','V1544','V1020','V1252','V1118','V1117', 'V1295')


column_mapping <- c(
  "FSC PAR"        = "FSC.PAR",
  "SSC"            = "SSC",
  "Hoechst.Red.DNA"= "Hoechst.Red",   
  "Hoechst.Red"    = "Hoechst.Red",
  "Hoechst Red.DNA"    = "Hoechst.Red",
  "FITC.hIgA2"     = "FITC",
  "FITC"           = "FITC",
  "APC.hIgA1"      = "APC",
  "APC"            = "APC",
  "Pe-TR.hIgG"     = "Pe.TR",
  "Pe.TR.hIgG"     = "Pe.TR",
  "Pe.TR"          = "Pe.TR",
  "BV650.hIgM"     = "BV650",
  "BV650"          = "BV650"
)



plot_params <- list(
  list(x = "FSC.PAR", y = "SSC"),
  list(x = "FSC.PAR", y = "Hoechst.Red"),
  list(x = "FITC", y = "APC"),
  list(x = "Pe.TR", y = "BV650")
)

dat <- MBR_process(
  fcs_files,
  file_index = 1,
  transformation = function(x) 10^((4 * x) / 65000),
  column_mapping = column_mapping
)

bins <- MBR_prepare(
  cohonen_information,
  selected_rows = selected_rows,
  transformation = function(x) 10^((4 * x) / 65000),
  column_mapping = column_mapping
)


plots <- MBR_plot(
  dat = dat,
  bins = bins,
  plot_params = plot_params,
  selected_rows = selected_rows
)

print(plots)
```



