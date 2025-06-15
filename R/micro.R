#' @import pbapply
#' @import kohonen
#' @import reshape2
#' @import circlize
#' @import ggplot2
#' @import patchwork
#' @import ggsci
#' @import vegan
#' @import ggpubr
#' @import permute
#' @import lattice
#' @import dplyr
#' @import rstatix
#' @import forcats
#' @importFrom cowplot plot_grid
#' @importFrom caret rfe rfeControl rfFuncs trainControl twoClassSummary train
#' @importFrom ggh4x force_panelsizes
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom linkET mantel_test qcorrplot geom_square geom_mark geom_couple correlate nice_curvature color_pal
#' @importFrom pROC roc ci.auc ggroc ci.se
#' @import tibble
#' @import viridis
#' @import scales

#' @export
MBR_som <- function(gated_fcs = NULL, fcs_data = NULL, m = 2e3, n_hclust = 300,
              n_cells_sub = 3e5, cluster_to_show = NULL, out_path = './') {
  
  # ================================================================
  # Stage 1: Environment Setup and Validation
  # ================================================================
  cat("==========================================\n")
  cat("Starting Analysis Pipeline\n")
  cat("==========================================\n")
  
  # Create output directory
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
    cat("✓ Created output directory:", out_path, "\n")
  }
  
  # ================================================================
  # Stage 2: Random Sampling
  # ================================================================
  cat("\nStep 1: Performing Random Sampling\n")
  cat("----------------------------------------\n")
  cat("Target sampling size:", format(n_cells_sub, scientific = FALSE), "cells\n")
  
  # Execute sampling
  sub_sample <- downsample(
    gated_fcs, 
    n = n_cells_sub, 
    samples = names(gated_fcs) 
  )
  sub_sample <- list(name = "random sample", data = sub_sample)
  
  cat("✓ Sampling completed, obtained", nrow(sub_sample$data), "cells\n")
  
  # ================================================================
  # Stage 3: SOM Training
  # ================================================================
  cat("\nStep 2: Training Self-Organizing Map (SOM)\n")
  cat("----------------------------------------\n")
  cat("Calculating SOM grid size...\n")
  
  # Calculate SOM parameters
  cells_per_node <- nrow(sub_sample$data) / m
  cat("Expected cells per node:", round(cells_per_node, 2), "\n")
  
  # Train SOM
  sub_sample <- compute_som(sub_sample, n_cells = cells_per_node)
  
  cat("✓ SOM training completed\n")
  cat("SOM grid size:", dim(sub_sample$som$grid$pts)[1], "nodes\n")
  
  # ================================================================
  # Stage 4: Save SOM Information
  # ================================================================
  cat("\nStep 3: Saving SOM Codebook Information\n")
  cat("----------------------------------------\n")
  
  # Save SOM codebook
  cohonen_information <- as.data.frame(sub_sample[["som"]][["codes"]][[1]]) 
  assign("cohonen_information", cohonen_information, envir = .GlobalEnv) 
  
  som_file <- file.path(out_path, "cluster_information.csv")
  write.csv(cohonen_information, som_file) 
  
  cat("✓ SOM codebook saved to:", som_file, "\n")
  cat("Codebook dimensions:", dim(cohonen_information), "\n")
  
  # ================================================================
  # Stage 5: Sample Quality Control
  # ================================================================
  cat("\nStep 4: Sample Quality Control\n")
  cat("----------------------------------------\n")
  
  min_cells <- 2e5 
  original_sample_count <- length(gated_fcs)
  
  # Check sample cell counts
  cell_counts <- sapply(gated_fcs, function(x) nrow(x$data))
  low_count_samples <- cell_counts < min_cells
  
  if (any(low_count_samples)) {
    dropped_samples <- names(gated_fcs)[low_count_samples]
    warning(paste0("Dropping low cell count samples: ", 
                   paste(dropped_samples, collapse = ", "), 
                   " (cell count < ", format(min_cells, scientific = FALSE), ")"))
    cat("⚠ Dropped", length(dropped_samples), "samples\n")
  }
  
  # Filter samples
  gated_fcs <- gated_fcs[!low_count_samples]
  cat("✓ Retained", length(gated_fcs), "/", original_sample_count, "samples\n")
  
  # ================================================================
  # Stage 6: Mapping to SOM
  # ================================================================
  cat("\nStep 5: Mapping All Samples to Trained SOM\n")
  cat("----------------------------------------\n")
  
  n_subset_large <- 1e10  # Use all cells
  
  cat("Mapping", length(gated_fcs), "samples to SOM...\n")
  
  # Map to trained SOM (with progress bar)
  SOM_fcs <- pbapply::pblapply(gated_fcs, map_som, 
                               trained = sub_sample$som, 
                               n_subset = n_subset_large) 
  
  cat("✓ SOM mapping completed\n")
  
  # ================================================================
  # Stage 7: Cluster Assignment
  # ================================================================
  cat("\nStep 6: Assigning Cluster Labels\n")
  cat("----------------------------------------\n")
  
  # Create cluster number matrix
  cluster_number <- as.matrix(as.list(1:2025)) 
  colnames(cluster_number) <- as.character(2025) 
  
  cat("Total clusters:", ncol(cluster_number), "\n")
  cat("Assigning cluster labels...\n")
  
  # Assign clusters (with progress bar)
  SOM_fcs <- pbapply::pblapply(SOM_fcs, assign_clusters, clusters = cluster_number) 
  
  cat("✓ Cluster assignment completed\n")
  
  # ================================================================
  # Stage 8: Counting and Statistics
  # ================================================================
  cat("\nStep 7: Computing Cluster Statistics\n")
  cat("----------------------------------------\n")
  
  cat("Calculating cell counts for each cluster...\n")
  
  # Count observations (with progress bar)
  SOM_fcs <- pbapply::pblapply(SOM_fcs, count_observations, 
                               clusters = colnames(cluster_number))
  
  # Get count tables
  count_tables <- get_counts(SOM_fcs) 
  
  cat("✓ Counting completed\n")
  
  # ================================================================
  # Stage 9: Save Count Results
  # ================================================================
  cat("\nStep 8: Saving Count Results\n")
  cat("----------------------------------------\n")
  
  # Save count tables
  tmp_path <- file.path(out_path, "count_tables")
  if (!dir.exists(tmp_path)) dir.create(tmp_path, recursive = TRUE)
  
  count_file <- file.path(tmp_path, "count_table_SOM.save")
  save(count_tables, file = count_file, compress = TRUE) 
  
  cat("✓ Raw count table saved to:", count_file, "\n")
  
  # Process and normalize counts
  original <- as.data.frame(t(count_tables[[1]])) 
  count_table_2025 <- original / rowSums(original) 
  
  # Check for rare clusters
  rare_clusters <- sum(sapply(count_table_2025, min) < 1e-4)
  if (rare_clusters > 0) { 
    warning(paste0(rare_clusters, ' clusters contain less than 0.01% of cells!'))
    cat("⚠ Found", rare_clusters, "rare clusters\n")
  }
  
  # Save normalized count table
  som_file <- file.path(out_path, "SOM.csv")
  write.csv(count_table_2025, som_file) 
  
  cat("✓ Normalized count table saved to:", som_file, "\n")
  
  # Save to global environment
  assign("original", original, envir = .GlobalEnv) 
  assign("count_table", count_table_2025, envir = .GlobalEnv) 
  
  # ================================================================
  # Stage 10: Export FCS Files
  # ================================================================
  cat("\nStep 9: Exporting Mapped FCS Files\n")
  cat("----------------------------------------\n")
  
  # Create output directory
  dir_name <- file.path(out_path, "MappedFCS") 
  if (!dir.exists(dir_name)) { 
    dir.create(dir_name, recursive = TRUE) 
    cat("✓ Created FCS output directory:", dir_name, "\n")
  }
  
  cat("Exporting", length(SOM_fcs), "FCS files...\n")
  
  # Export each sample's FCS file
  pb <- txtProgressBar(min = 0, max = length(SOM_fcs), style = 3)
  
  for (i in seq_along(SOM_fcs)) { 
    # Prepare data
    fcs_data <- as.data.frame(SOM_fcs[[i]]$data) 
    fcs_data$classes <- as.numeric(SOM_fcs[[i]]$classes) 
    
    # Convert to matrix and create flow cytometry frame
    fcs_matrix <- as.matrix(fcs_data) 
    colnames(fcs_matrix) <- names(fcs_data) 
    current_frame <- flowFrame(fcs_matrix) 
    
    # Generate filename and save
    file_name <- file.path(dir_name, SOM_fcs[[i]]$name) 
    write.FCS(current_frame, file_name) 
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("\n✓ FCS file export completed\n")
  
  # ================================================================
  # Completion Summary
  # ================================================================
  cat("\n==========================================\n")
  cat("Analysis Pipeline Completed!\n")
  cat("==========================================\n")
  cat("Processed samples:", length(SOM_fcs), "\n")
  cat("Total clusters:", ncol(cluster_number), "\n")
  cat("Output directory:", out_path, "\n")
  cat("\nMain output files:\n")
  cat("- SOM codebook:", file.path(out_path, "cluster_information.csv"), "\n")
  cat("- Normalized counts:", file.path(out_path, "SOM.csv"), "\n")
  cat("- Raw counts:", file.path(tmp_path, "count_table_SOM.save"), "\n")
  cat("- Mapped FCS files:", dir_name, "\n")
  cat("==========================================\n")
}

#' @export
MBR_stat <- function(data = NULL, meta_data = NULL, group_col = NULL,
                   test_type = 'wilcox', cutoff = 0.05,
                   correction = 'none', out_path = './') {

  cat(paste0(test_type, " test\n"))
  if (length(unlist(meta_data[group_col])) != nrow(data)) {
    stop("Error: The lengths of group_labels and number of samples are not equal.")
  }

  if (test_type == 'wilcoxon') {
    pvalue <- sapply(data, function(t) {
      tryCatch({wilcox.test(t ~ as.vector(unlist(meta_data[group_col])))$p.value}, warning = function(e){1})})
  } else if (test_type == 'kruskal') {
    pvalue <- sapply(data, function(t) {
      tryCatch({kruskal.test(t ~ as.vector(unlist(meta_data[group_col])))$p.value}, warning = function(e){1})})
  } else if (test_type == 'anova') {
    pvalue <- sapply(data, function(t) {
      tryCatch({anova(aov(t ~ as.vector(unlist(meta_data[group_col]))))$`Pr(>F)`[1]}, warning = function(e){1})})
  } else {
    stop("Error: test_type should be one of 'wilcoxon', 'kruskal', 'anova'.")
  }

  if (correction == 'fdr') {
    pvalue <- p.adjust(pvalue, method = 'fdr')
  } else if (correction == 'bonferroni') {
    pvalue <- p.adjust(pvalue, method = 'bonferroni')
  } else if (correction == 'BH') {
    pvalue <- p.adjust(pvalue, method = 'BH')
  } else if (correction != 'none') {
    stop('Error: correction should be one of "none", "fdr", "bonferroni", "BH".')
  }
  pvalue_data = as.data.frame(pvalue) %>% replace(is.na(.), 1)
  significant_data = data[, pvalue < cutoff]

  assign("pvalue_data", pvalue_data, envir = .GlobalEnv)
  assign("significant_data", significant_data, envir = .GlobalEnv)

  write.csv(pvalue_data, paste0(out_path, "/Ig_SOM_1024_pvalue.csv"))
  write.csv(significant_data, paste0(out_path, "/Ig_SOM_1024_significant.csv"))
  cat("\tDONE.\n")
}

#' @export
MBR_circle <- function(data = NULL, meta_data = NULL, group_col = NULL,
                     out_path = './', width = 8, height = 8,
                     point_colors = c("red", "blue"),
                     cell_colors = c("blue", "white", "red")) {
  cat("Circled heatmap\n")
  group_labels <- as.vector(unlist(meta_data[group_col]))
  mean_table <- t(sapply(data, function(t) tapply(t, group_labels, mean)))
  rownames(mean_table) <- gsub('V', 'Bin', rownames(mean_table))
  min_val <- min(mean_table)
  max_val <- max(mean_table)
  mean_val <- (min_val + max_val) / 2
  col_fun = colorRamp2(c(min_val, mean_val, max_val), cell_colors)

  circos.clear()
  circos.par(gap.degree = 10)
  circos.heatmap(mean_table, col = col_fun,
                 rownames.side = 'outside')

  mean_diff <- as.vector(mean_table[,1] - mean_table[,2])
  suppressMessages({
    circos.track(ylim = range(mean_diff),
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   y = mean_diff[CELL_META$row_order]
                   circos.points(seq_along(y) - 0.5, y, cex = 1,
                                 col = ifelse(y > 0, point_colors[1], point_colors[2]))
                   circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                   # add column label
                   circos.text(c(-0.8, -0.8),
                               c(1.8 * max(mean_diff) - 1.8 * min(mean_diff),
                                 3.2 * max(mean_diff) - 3.2 * min(mean_diff)),
                               rev(unique(group_labels)),
                               facing = 'downward',
                               cex = 0.5)
                 },
                 cell.padding = c(0.02, 0, 0.02, 0))
  })

  pdf(paste0(out_path, '/circle.pdf'), height = height, width = width)
  circos.clear()
  circos.par(gap.degree = 10)
  circos.heatmap(mean_table, col = col_fun,
                 rownames.side = 'outside')

  mean_diff <- as.vector(mean_table[,1] - mean_table[,2])
  suppressMessages({
    circos.track(ylim = range(mean_diff),
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   y = mean_diff[CELL_META$row_order]
                   circos.points(seq_along(y) - 0.5, y, cex = 1,
                                 col = ifelse(y > 0, point_colors[1], point_colors[2]))
                   circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                   # add column label
                   circos.text(c(-0.8, -0.8),
                               c(1.8 * max(mean_diff) - 1.8 * min(mean_diff),
                                 3.2 * max(mean_diff) - 3.2 * min(mean_diff)),
                               rev(unique(group_labels)),
                               facing = 'downward',
                               cex = 0.5)
                 },
                 cell.padding = c(0.02, 0, 0.02, 0))
  })
  invisible(dev.off())
  cat("\tDONE.\n")
}

#' @export
MBR_violin <- function(data = NULL, meta_data = NULL, group_col = NULL,
                       pvalue_data = NULL, cluster = 1, out_path = './',
                       colors = c('#E41A1C', '#377EB8'),
                       width = 4, height = 4) {
  cat("Creating violin plot...\n")
  
  # Extract group labels from metadata
  group_labels <- as.vector(unlist(meta_data[[group_col]]))
  
  # Identify clusters by removing 'v' prefix from column names
  clusters <- as.numeric(gsub('v', '', colnames(data), ignore.case = TRUE))
  
  # Extract data for the requested cluster
  cluster_data <- data[, clusters == cluster]
  
  # Get p-value for the selected cluster and round it
  p_value <- format(pvalue_data[clusters == cluster, ], digits = 3, scientific = TRUE)
  
  # Determine max value for positioning the p-value label
  max_val <- max(cluster_data, na.rm = TRUE)
  if (!is.finite(max_val)) {
    max_val <- 1
  }
  
  # Create the violin plot with ggplot2
  p <- ggplot(mapping = aes(x = group_labels, y = cluster_data)) +
    geom_violin(aes(fill = group_labels), alpha = 0.5) +
    geom_boxplot(aes(fill = group_labels), width = 0.1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = colors) +
    labs(title = paste0('Cluster ', cluster),
         x = 'Group',
         y = 'Relative Abundance',
         fill = 'Group') +
    geom_text(aes(label = ifelse(p_value < 0.05, paste0('P = ', p_value, '*'), paste0('P = ', p_value))),
              x = 1.5, y = max_val * 0.95, size = 4)
  
  # Print plot to console and save to PDF
  print(p)
  pdf(file.path(out_path, 'violin.pdf'), width = width, height = height)
  print(p)
  dev.off()
  
  cat("Violin plot complete.\n")
}


#' @export
MBR_beta <- function(data, out_path = './', test = 'wilcox.test', meta_data = NULL,
                   group_name = NULL, colors = c('#E41A1C', '#377EB8'),
                   width = 5, height = 5) {

  df <- cbind(meta_data[group_name], data)
  colnames(df)[1] <- 'Health.State'
  rownames(df) <- NULL
  groups<-as.data.frame(cbind(sample=paste('sample',rownames(df),sep = ''),group=df$Health.State))
  write.table(groups, paste0(out_path, '/.group.txt'), row.names = F,sep = '\t',quote = F)
  df$Health.State<-groups$sample
  rownames(df)<-df$Health.State
  dataT<-df[,-1]

  dist <- vegdist(dataT, method="bray")
  dist <- as.matrix(dist)
  adist<- as.dist(dist)

  options(stringsAsFactors=F)

  sd <- groups
  rownames(sd) <- as.character(sd[,1])
  sd[,2] <- as.character(sd[,2])

  dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]

  pc_num <-c(1,2)
  pc_x <- pc_num[1]
  pc_y <- pc_num[2]

  pcoa <- cmdscale(dist, k=3, eig=TRUE)
  pc12 <- pcoa$points[,pc_num]
  pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits = 2)
  pc12 <- as.data.frame(pc12)
  colnames(pc12) <- c("pc_x","pc_y")
  pc12['sample'] <- rownames(pc12)
  colnames(sd)[1:2] <- c("sample","group")
  sd$group<-factor(sd$group,levels=sd$group[!duplicated(sd$group)])
  pc12 <- merge(pc12,sd,by="sample")
  pc12$group<-factor(pc12$group,levels=levels(sd$group))

  p<- ggscatter(pc12, x = "pc_x", y = "pc_y",
               color = "group", shape = "group", linewidth=3,
               ellipse = TRUE, conf.int.level = 0.95,
               #palette = colors,
               mean.point = TRUE,
               star.plot = TRUE,star.plot.lty = 1) +
    #ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
    #xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
    geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "solid")+
    geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "solid")+
    scale_color_manual(values = colors) +

    theme(axis.title = element_blank(),
          legend.position = "top",legend.title = element_blank(),
          panel.border = element_rect(color = "black",linewidth  = 1.0, fill = NA),
          text = element_text(size=12)) +theme(plot.margin = unit(c(0,0,0,0),'cm'))

  mycols=pal_npg("nrc")(10)[1:length(levels(sd$group))]

  ADONIS<-suppressMessages(adonis2(dist~sd$group))

  TEST<-ADONIS$`Pr(>F)`[1]
  R2adonis<-round(ADONIS$R2[1],digits = 3)
  sink(paste0(out_path, '/adonis.txt'))
  print(ADONIS)
  sink()

  p<-p+ggtitle(label =expr(paste(bold(Bray)," ",bold(Curtis),bold(","),
                                 bold(Adnois:R^2),bold('='),bold(!!R2adonis),
                                 bold(","),bold(P),bold('='),!!TEST))) +
    theme(title = element_text(size = 10))

  cp <- combn(levels(pc12$group),2)
  comp <- list()
  for(i in 1:ncol(cp)){
    comp[[i]] <- cp[,i]
  }

  pl<-ggboxplot(pc12, x="group", y="pc_y", fill = "group", palette = colors) +
    stat_compare_means(comparisons = comp, label = "p.signif",method="wilcox.test")+
    ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
    theme(panel.border = element_rect(color = "black",linewidth = 1.0,fill = NA),
          #axis.text.y = element_blank(),
          #axis.ticks.y= element_blank(),
          #axis.title.y = element_blank(),
          legend.position = "none",
          axis.text = element_text(face='bold'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12,angle = 60,hjust = 1,face='bold'),
          axis.title.y = element_text(size = 15,face='bold'))+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm'))

  pt<-ggboxplot(pc12, x="group", y="pc_x", fill = "group", palette = colors) + coord_flip() +
    stat_compare_means(comparisons = comp, label = "p.signif",method = test) +
    scale_x_discrete(limits = rev(levels(pc12$group)))+
    ylab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
    theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
          #axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          #axis.title.x = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 12, angle = 0,face='bold'),
          axis.title.y=element_blank(),
          axis.title.x = element_text(size = 15,face='bold'))+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm'))

  p0 <- ggplot() + theme(panel.background = element_blank(),
                         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"))

  p_final<- pl+p+p0+pt+plot_layout(ncol = 2,nrow = 2,heights = c(4,1),widths = c(1,4))+
    plot_annotation(theme = theme(plot.margin = margin()))

  ggsave(paste0(out_path, '/pcoa.pdf') ,width = width, height = height)
  return(p_final)
  cat("\tDONE.\n")
}

#' @export
MBR_fs <- function(data = NULL, out_path = './',
                   meta_data = NULL,
                   group_name = NULL,
                   nfolds_cv = 5,
                   top_n_features = 20,
                   rfe_size = 10,
                   ref_group = NULL,
                   colors = c('#E41A1C', '#377EB8'),
                   width = 5, height = 5) {

  df <- cbind(meta_data[group_name], data)
  colnames(df)[1] <- 'Health.State'
  rownames(df) <- NULL

  groups <- data.frame(sample = paste0('sample', seq_len(nrow(df))), group = df$Health.State)
  write.table(groups, file = file.path(out_path, '.group.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

  df$Health.State <- groups$sample
  rownames(df) <- df$Health.State
  df <- df[, -1]
  df <- as.data.frame(t(df))
  df2 <- as.data.frame(t(df))
  df2$group <- factor(groups$group)

  control <- rfeControl(functions = rfFuncs, method = "cv", number = nfolds_cv)
  results <- rfe(df2[, 1:(ncol(df2) - 1)], df2[, ncol(df2)], sizes = 1:rfe_size, rfeControl = control)

  print(results)
  predictors(results)
  p1_out <- plot(results, type = c("g", "o"), main = 'Feature Selection', col = '#377EB8', lwd = 2)

  best_selection <- results$optVariables
  message("The number of selected features is ", length(best_selection))

  df2 <- df2[, c(best_selection, 'group')]
  write.table(df2, file = file.path(out_path, 'figure1.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

  # Save selected features ONLY (no group) from input data to global environment
  MBR_selected_features <<- data[, best_selection, drop = FALSE]
  rownames(MBR_selected_features) <<- rownames(data)

  varimp_all <- varImp(results)
  varimp_df <- varimp_all[best_selection, , drop = FALSE] 

  top_n <- min(top_n_features, nrow(varimp_df))
  varimp_data <- data.frame(
    feature = rownames(varimp_df)[1:top_n],
    importance = varimp_df[1:top_n, 1]
  )
  top_features <- varimp_data$feature
  df2 <- df2[, c(top_features, 'group')]

  p2_out <- ggplot(drop_na(varimp_data),
                   aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
    geom_bar(stat = "identity") +
    labs(x = "Features", y = "Variable Importance") +
    geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "Variable Importance")

  esize <- df2 %>%
    reshape2::melt() %>%
    mutate(value = log(value + 1, 10)) %>%
    group_by(variable) %>%
    suppressMessages() %>%
    cohens_d(value ~ group, conf.level = 0.95, ci = TRUE, ref.group = ref_group) %>%
    arrange(desc(effsize))

  fig2 <- ggplot(esize, aes(y = reorder(variable, effsize), x = effsize, fill = variable)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
    geom_point(aes(x = effsize), size = 4, color = 'black') +
    geom_point(aes(x = effsize, color = ifelse(effsize > 0, 'positive', 'negative')), size = 3) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_color_manual(values = c('positive' = colors[1], 'negative' = colors[2])) +
    labs(title = "Effect Size", x = "Cohen's d", y = '') +
    force_panelsizes(rows = 0.5, cols = 0.5)

  fig1 <- df2 %>%
    reshape2::melt() %>%
    mutate(value = log(value + 1, 10)) %>%
    ggplot(aes(x = value, y = fct_rev(factor(variable, esize$variable)), fill = group)) +
    geom_boxplot(position = 'dodge') +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = "log10(abundance)", y = "Variable") +
    scale_fill_manual(values = colors) +
    force_panelsizes(rows = 0.5, cols = 0.5)

  fig3 <- merge(varimp_data, esize, by.x = 'feature', by.y = 'variable') %>%
    ggplot(aes(x = importance, y = fct_rev(factor(feature, esize$variable)))) +
    geom_bar(stat = "identity", color = 'black', aes(fill = ifelse(effsize > 0, 'positive', 'negative'))) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = c('positive' = colors[1], 'negative' = colors[2])) +
    labs(title = "Variable Importance", y = '') +
    force_panelsizes(rows = 0.5, cols = 0.5)

  p3_out <- plot_grid(fig1, fig2, fig3, ncol = 3, align = 'h', axis = 't')

  pdf(file.path(out_path, 'feature_exploration.pdf'), width = width, height = height)
  print(p1_out)
  print(p2_out)
  print(p3_out)
  dev.off()

  show_plots <- function(plots) {
    for (plot in plots) {
      print(plot)
      readline(prompt = "Press [enter] to see the next plot")
    }
  }

  show_plots(list(p1_out, p2_out, p3_out))
  cat("\tDONE.\n")
}

#' @export
MBR_heatmap <- function(data = NULL, cohonen_information = NULL,
                      out_path = './',
                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                      width = 10, height = 10,
                      scale = 'row', cluster_rows = F, cluster_cols = F,
                      display_numbers = T, boarder_color = 'grey60',
                      legend = T) {
  cluster_ids <- paste0('V', gsub('^[vV]', '', colnames(data)))
  rownames(cohonen_information) <- as.character(rownames(cohonen_information))
  df <- cohonen_information[cluster_ids, , drop = FALSE]
  df_mat <- as.matrix(sapply(df, as.numeric))
  rownames(df_mat) <- rownames(df)                    	
  p <- pheatmap(df, scale = scale, cluster_rows = cluster_rows,
           cluster_cols = cluster_cols, display_numbers = display_numbers,
           border_color = boarder_color, color = color, legend = legend)
  pdf(paste0(out_path, '/heatmap.pdf'), width = width, height = height)
  print(p)
  dev.off()
  print(p)
  cat("\tDONE.\n")
}



#' @export
MBR_mantel <- function(data = NULL, meta_data = NULL,
                       clinical_cols = NULL, demographic_cols = NULL,
                       spec_select_names = list(A = "A", B = "B"),
                       colors = RColorBrewer::brewer.pal(11, "RdBu"),
                       out_path = './', width = 8, height = 8) {

  spec_select <- list()
  if (!is.null(spec_select_names$A)) {
    spec_select[[ spec_select_names$A ]] <- clinical_cols
  }
  if (!is.null(spec_select_names$B)) {
    spec_select[[ spec_select_names$B ]] <- demographic_cols
  }

  mantel <- suppressMessages(
    mantel_test(meta_data, data, spec_select = spec_select) %>%
    mutate(
      rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
               labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
      pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
    )
  )
  
    p <- qcorrplot(correlate(data), type = "lower", diag = FALSE) +
    geom_square() +
    geom_mark(sep='\n', size = 1.8, sig_level = c(0.05, 0.01, 0.001),
              sig_thres = 0.05, color="white")+
    geom_couple(aes(colour = pd, size = rd),
                data = mantel,
                curvature = nice_curvature()) +
    scale_fill_gradientn(colours = colors) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(size = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"),
                               order = 2),
           colour = guide_legend(title = "Mantel's p",
                                 override.aes = list(size = 3),
                                 order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3))

  #pdf
  pdf(paste0(out_path, '/mantel.pdf'), width = width, height = height)
  print(p)
  dev.off()

  print(p)
  cat("\tDONE.\n")
}

#' @export
MBR_ml <- function(data = NULL, meta_data = NULL,
                   group_name = 'Group', out_path = './',
                   reference_level = 'A',
                   width = 5, height = 5,
                   method = "repeatedcv", number = 5, repeats = 5) {
  
  rownames(data) <- NULL
  
  # Create groups dataframe with column named 'group' (fixed name)
  groups <- data.frame(
    sample = paste('sample', seq_len(nrow(data)), sep = ''),
    group = meta_data[[group_name]]
  )
  
  write.table(groups, paste0(out_path, 'ml_group.txt'), row.names = FALSE, sep = '\t', quote = FALSE)
  
  # Add sample names as a new column for matching, then transpose
  data <- cbind(data, Health.State = groups$sample)
  rownames(data) <- data$Health.State
  data <- data[, -ncol(data)]  # remove Health.State column for modeling
  df <- as.data.frame(t(data))
  train <- as.data.frame(t(df))
  
  # Fix here: assign factor groups from groups$group (not groups[group_name])
  train$group <- factor(groups$group)
  
  # Load pROC inside function
  library(pROC)
  library(caret)
  library(ggplot2)
  library(tibble)
  
  fitControl <- trainControl(
    method = method,
    number = number,
    repeats = repeats,
    returnResamp = "final",
    classProbs = TRUE,
    savePredictions = TRUE,
    summaryFunction = twoClassSummary
  )
  
  rf <- train(
    group ~ .,
    data = train,
    method = "rf",
    trControl = fitControl,
    metric = "ROC",
    verbose = FALSE
  )
  
  rocs_train <- roc(response = ifelse(train$group == reference_level, 0, 1), predictor = rf$finalModel$votes[, 2])
  ci_auc_train <- suppressWarnings(round(as.numeric(ci.auc(rocs_train)), 3))
  ci_tb_train <- suppressWarnings(as.data.frame(ci.se(rocs_train)))
  ci_tb_train <- suppressWarnings(rownames_to_column(ci_tb_train, var = 'x'))
  ci_tb_train <- as.data.frame(sapply(ci_tb_train, as.numeric))
  names(ci_tb_train) <- c('x', 'low', 'mid', 'high')
  
  metrics_train <- confusionMatrix(rf$finalModel$predicted, train$group,
                                   positive = reference_level, mode = 'everything')
  
  options(warn = -1)
  g1_out <- ggroc(rocs_train, legacy.axes = TRUE) +
    coord_equal() +
    geom_ribbon(aes(x = 1 - x, ymin = low, ymax = high), data = ci_tb_train, alpha = 0.5, fill = 'lightblue') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', alpha = 0.7) +
    geom_text(aes(0.5, 0.25, hjust = 0,
                  label = paste0('AUC: ', round(rocs_train$auc, 3), ' 95%CI: ',
                                 ci_auc_train[1], ' ~ ', ci_auc_train[3]))) +
    geom_text(aes(0.5, 0.2, hjust = 0,
                  label = paste0('Sensitivity: ', round(as.numeric(metrics_train$byClass[1]), 3)))) +
    geom_text(aes(0.5, 0.15, hjust = 0,
                  label = paste0('Specificity: ', round(as.numeric(metrics_train$byClass[2]), 3)))) +
    geom_text(aes(0.5, 0.1, hjust = 0,
                  label = paste0('F1: ', round(as.numeric(metrics_train$byClass[7]), 3)))) +
    theme_classic() +
    labs(x = '1 - Specificity (% false positive)',
         y = 'Sensitivity (% true positive)')
  
  pdf(paste0(out_path, 'roc_confusion.pdf'), width = width, height = height)
  print(g1_out)
  dev.off()
  
  print(g1_out)
  cat("\tDONE.\n")
}


#' @export
MBR_conf <- function(data = NULL, meta_data = NULL,                     group_name = 'Group', out_path = './',                     reference_level = 'A',                     colors = c('#00BFC4', '#F8766D'),                     width = 5, height = 5,                     method = "repeatedcv", number = 5, repeats = 5) {  # Generate sample names  rownames(data) <- NULL  groups <- data.frame(sample = paste0('sample', seq_len(nrow(data))),                       group = meta_data[[group_name]])  write.table(groups, paste0(out_path, 'ml_group.txt'), row.names = FALSE, sep = '\t', quote = FALSE)    # Prepare data  data <- cbind(data, Health.State = groups$sample)  rownames(data) <- data$Health.State  data <- data[, -ncol(data)]  df <- as.data.frame(t(data))  train <- as.data.frame(t(df))  train$group <- factor(groups$group)    # Train control  fitControl <- trainControl(    method = method,    number = number,    repeats = repeats,    returnResamp = "final",    classProbs = TRUE,    savePredictions = TRUE,    summaryFunction = twoClassSummary  )    # Train model  rf <- train(    group ~ .,    data = train,    method = "rf",    trControl = fitControl,    metric = "ROC",    verbose = FALSE  )    # Confusion matrix  pred <- rf$finalModel$predicted  truth <- train$group  cm <- confusionMatrix(pred, truth, positive = reference_level)  cm_table <- as.data.frame(cm$table)    # Plot  library(ggplot2)  library(reshape2)    g2_out <- ggplot(cm_table, aes(x = Reference, y = Prediction)) +    geom_tile(aes(fill = log(Freq + 1))) +    geom_text(aes(label = Freq)) +    scale_fill_gradient2(low = colors[1], high = colors[2],                         midpoint = mean(log(cm_table$Freq + 1))) +    coord_equal() +    theme_minimal() +    labs(title = 'Confusion Matrix',         x = 'True Label',         y = 'Predicted Label',         fill = 'Log(Count)') +    theme(legend.position = 'right')    pdf(file = paste0(out_path, 'confusion_matrix.pdf'), width = width, height = height)  print(g2_out)  dev.off()    print(g2_out)  cat("\tDONE.\n")}

#' @export
MBR_reclustering <- function(data, num_clusters) {  if (missing(data) || !is.data.frame(data)) {    stop("Please provide a valid data frame as 'data'.")  }    features <- data[, sapply(data, is.numeric)]  if (ncol(features) == 0) {    stop("No numeric columns found in the input data.")  }    features_scaled <- scale(features)  dist_matrix <- dist(features_scaled, method = "euclidean")  hc <- hclust(dist_matrix, method = "ward.D2")    cluster_assignments <- cutree(hc, k = num_clusters)    reclustered_information <- data  reclustered_information$Cluster <- cluster_assignments    assign("reclustered_information", reclustered_information, envir = .GlobalEnv)    cat("Reclustering complete. 'reclustered_information' created with dimensions:\n")  print(dim(reclustered_information))    invisible(reclustered_information)}


#' @export
MBR_read <- function(rawdata_path) {
  # Check if the path exists
  if (!dir.exists(rawdata_path)) {
    stop("Directory does not exist: ", rawdata_path)
  }

  # List FCS files
  fcs_path <- list.files(rawdata_path, pattern = "fcs$", full.names = TRUE)

  if (length(fcs_path) == 0) {
    stop("No FCS files found in directory: ", rawdata_path)
  }

  # Read FCS files as flowSet
  fcs_files <- read.flowSet(files = fcs_path,
                            column.pattern = "\\*|Bits|Drop",
                            invert.pattern = TRUE,
                            alter.names = FALSE)

  return(fcs_files)
}

#' @export
MBR_process <- function(fcs_files,
                              file_index = 1,
                              transformation = function(x) 10^((4 * x) / 65000),
                              column_mapping = NULL) {

  # Check if file_index is valid
  if (file_index > length(fcs_files) || file_index < 1) {
    stop("Invalid file index. Must be between 1 and ", length(fcs_files))
  }

  # Get file name from flowSet
  file_name <- sampleNames(fcs_files)[file_index]

  # Extract data from specified file
  dat <- exprs(fcs_files@frames[[file_name]]) %>%
    as.data.frame()

  # Apply transformation to all columns except 'classes' if it exists
  if ("classes" %in% colnames(dat)) {
    dat <- dat %>%
      mutate_at(vars(-'classes'), transformation)
  } else {
    dat <- dat %>%
      mutate_all(transformation)
  }

  # Rename columns if mapping is provided
  if (!is.null(column_mapping)) {
    new_names <- colnames(dat)
    for (old_name in names(column_mapping)) {
      if (old_name %in% colnames(dat)) {
        new_names[which(colnames(dat) == old_name)] <- column_mapping[old_name]
      }
    }
    colnames(dat) <- new_names
  }

  return(dat)
}

#' @export
MBR_prepare <- function(bins,
                         selected_rows,
                         transformation = function(x) 10^((4 * x) / 65000),
                         column_mapping = NULL) {

  # Select rows
  bins <- bins[rownames(bins) %in% selected_rows, , drop = FALSE]

  # Convert to dataframe if not already
  bins <- as.data.frame(bins)

  # Rename columns if mapping is provided
  if (!is.null(column_mapping)) {
    new_names <- colnames(bins)
    for (old_name in names(column_mapping)) {
      if (old_name %in% colnames(bins)) {
        new_names[which(colnames(bins) == old_name)] <- column_mapping[old_name]
      }
    }
    colnames(bins) <- new_names
  }

  # Apply transformation
  bins <- bins %>%
    mutate_all(transformation)

  return(bins)
}

#' @export
create_flow_plot <- function(dat,
                             bins = NULL,
                             x,
                             y,
                             selected_rows = NULL,
                             x_limits = c(0.9, 11000),
                             y_limits = c(0.9, 11000),
                             hex_bins = 100,
                             point_color = "grey20",
                             point_alpha = 0.3,
                             point_size = 0.5,
                             label_color = "white",
                             label_size = 3,
                             theme_family = "Times") {

  # Create the base plot
  p <- dat %>%
    ggplot(aes_string(x = x, y = y)) +
    geom_hex(bins = hex_bins) +
    scale_fill_viridis(discrete = FALSE, trans = 'log') +
    scale_x_log10(breaks = c(0, 10, 100, 1000, 10000),
                  labels = trans_format("log10", scales::math_format(10^.x)),
                  limits = x_limits) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                  labels = trans_format("log10", scales::math_format(10^.x)),
                  limits = y_limits) +
    theme_bw() +
    theme(text = element_text(family = theme_family))

  # Add highlighted points if selected_rows is provided
  if (!is.null(selected_rows) && "classes" %in% colnames(dat)) {
    p <- p +
      geom_point(data = dat[paste0('V', dat[,"classes"]) %in% selected_rows, ],
                 aes_string(x = x, y = y),
                 color = point_color,
                 alpha = point_alpha,
                 size = point_size)
  }

  # Add bin labels if bins is provided
  if (!is.null(bins) && !is.null(selected_rows)) {
    p <- p +
      annotate('text',
               x = bins[, x],
               y = bins[, y],
               label = selected_rows,
               color = label_color,
               size = label_size)
  }

  return(p)
}

#' @export
MBR_plot <- function(dat,
                                  bins = NULL,
                                  plot_params,
                                  selected_rows = NULL,
                                  ncol = 2,
                                  nrow = 2,
                                  common_legend = TRUE,
                                  ...) {

  # Create a list to store plots
  plots <- list()

  # Create each plot
  for (i in seq_along(plot_params)) {
    params <- plot_params[[i]]
    if (!("x" %in% names(params)) || !("y" %in% names(params))) {
      stop("Each plot_params item must contain 'x' and 'y' parameters")
    }

    plots[[i]] <- create_flow_plot(dat = dat,
                                   bins = bins,
                                   x = params$x,
                                   y = params$y,
                                   selected_rows = selected_rows,
                                   ...)
  }

  # Arrange plots in a grid
  combined_plot <- ggarrange(plotlist = plots,
                             ncol = ncol,
                             nrow = nrow,
                             common.legend = common_legend)

  return(combined_plot)
}