
##########################################################################################
#### Example Code for Pairwise-Differential Expression, ECDF Plots, and Density Plots ####
##########################################################################################

#The PairwiseFindMarkers function performs differential-expression testing for
#multiple batches (or whatever is metadata is specified in the split.by
#arguement), using Seurat's FindMarkers function. It is similar to Seurat's
#FindConservedMarkers, albeit with modified and additional functionality. If
#only_compare_same_split_idents = TRUE, then differential expression will only
#be performed between cells within the same batch, for each batches. If
#only_compare_same_split_idents = FALSE, then differential expression will be
#performed between cells for all permutations of batches (e.g., Batch 1 cells v.
#Batch 1 cells, Batch 1 cells v. Batch 2 cells, Batch 1 cells v. Batch 3 cells,
#etc.).

#In our Everetts, Worley, et al. 2021 paper, developmental time points are
#separate batched. Whenever biological signal is confounded by batch effects, we
#recommend setting only_compare_same_split_idents = FALSE to produce
#conservative results.

#When producing the ECDF and density plots in Figure 1 - figure supplement 3F-H,
#we set only_compare_same_split_idents = TRUE, as our biological signal was not
#confounded by batch effects.

PairwiseFindMarkers <- function(seurat_object, cells.1, cells.2, split.by, slot = "data", assay = "RNA", features = NULL,
                                only_compare_same_split_idents = FALSE, perform_de_direction_check = TRUE,
                                logfc.threshold = 0.15, min.pct = 0.15, min.diff.pct = -Inf, min.cells.group = 0,
                                p.adjust_n = NULL, ...){
  split_idents <- unique(seurat_object@meta.data[[split.by]])
  cells.1_split <- sapply(split_idents, FUN = function(x){cells.1[seurat_object[[split.by]][cells.1,] == x]})
  cells.2_split <- sapply(split_idents, FUN = function(x){cells.2[seurat_object[[split.by]][cells.2,] == x]})
  cells.1_split <- cells.1_split[sapply(cells.1_split, FUN = function(x){length(x) > 0})]
  cells.2_split <- cells.2_split[sapply(cells.2_split, FUN = function(x){length(x) > 0})]
  
  #only_compare_same_split_idents restricts DE testing such that cells are compared to each other only if they have the same split.by ident
  #if only_compare_same_split_idents is FALSE, then the function will determine all permutations of cells.1_split vs. cells.2_split for DE testing
  if (only_compare_same_split_idents){
    split_idents_intersect <- intersect(names(cells.1_split), names(cells.2_split))
    de_tests_permutations <- as.data.frame(matrix(split_idents_intersect, nrow = length(split_idents_intersect), ncol = 2))
  } else {
    de_tests_permutations <- expand.grid(names(cells.1_split), names(cells.2_split), stringsAsFactors = FALSE)
  }
  
  if (NROW(de_tests_permutations) == 0) {
    stop("Incompatible split.by: Separating either cells.1 or cells.2 by split.by results in zero cells.")
  }
  
  de_results <- apply(de_tests_permutations, MARGIN = 1, FUN = function(x){
    FindMarkers(GetAssayData(seurat_object, slot = slot, assay = assay), slot = slot,
                cells.1 = cells.1_split[[x[1]]],
                cells.2 = cells.2_split[[x[2]]],
                features = features,
                logfc.threshold = logfc.threshold, min.pct = min.pct, min.diff.pct = min.diff.pct, min.cells.group, ...)  
  })
  
  de_test_names <- apply(de_tests_permutations, MARGIN = 1, FUN = function(x){paste0(x[1], "_v_", x[2])})
  names(de_results) <- de_test_names
  
  #only proceed with genes that were DE in all test permutations
  #the following code will find the intersection of genes between all test permutations
  genes_from_de_tests <- lapply(de_results, FUN = rownames)
  intersect_de_genes <- Reduce(intersect, x = genes_from_de_tests)
  
  #change the name of the column containing the Bonferroni-corrected p-values
  de_results <- lapply(de_results, FUN = function(x){colnames(x)[colnames(x) == "p_val_adj"] <- "p_val_adj_Bonf"; return(x)})
  
  #calculate Benjamini-Hochberg- (a.k.a. FDR) and Bonferroni-adjusted p-values for DE genes in all test permutations
  #note that the number of comparisons (the "n" in p.adjust) can be set as an argument, but by default, is the number of DE genes (recommended)
  de_results <- lapply(de_results, FUN = function(x){x$p_val_adj_BH <- p.adjust(x$p_val, method = "BH",
                                                                                n = p.adjust_n %||% length(x$p_val)); return(x)})
  de_results <- lapply(de_results, FUN = function(x){x$p_val_adj_Bonf <- p.adjust(x$p_val, method = "bonferroni",
                                                                                  n = p.adjust_n %||% length(x$p_val)); return(x)})
  
  #shrink the results down to only include logFC and p-values, and only for DE genes that were common among all permutations
  de_results <- lapply(de_results, FUN = function(x){x[intersect_de_genes,c("avg_logFC", "p_val", "p_val_adj_BH", "p_val_adj_Bonf")]})
  
  #de_direction_check is to ensure that avg_logFC is the same for each gene between all tests
  #any gene which has a conflicting avg_logFC magnitudes will be thrown out
  if (perform_de_direction_check){
    de_direction_check <- do.call(cbind, lapply(de_results, FUN = function(x){x$avg_logFC}))
    de_direction_check <- apply(de_direction_check, MARGIN = 1, FUN = function(x){all(x > 0) | all(x < 0)})
    de_results <- lapply(de_results, FUN = function(x){x[de_direction_check,]})
  }
  
  #calculate statistics applied over the results of all tests (e.g., the maximum p-value among all DE tests conducted)
  MAX_p_val <- apply(do.call(cbind, lapply(de_results, FUN = function(x){x$p_val})), MARGIN = 1, FUN = max)
  MAX_p_val_adj_BH <- apply(do.call(cbind, lapply(de_results, FUN = function(x){x$p_val_adj_BH})), MARGIN = 1, FUN = max)
  MAX_p_val_adj_Bonf <- apply(do.call(cbind, lapply(de_results, FUN = function(x){x$p_val_adj_Bonf})), MARGIN = 1, FUN = max)
  MEAN_avg_logFC <- log(apply(exp(do.call(cbind, lapply(de_results, FUN = function(x){x$avg_logFC}))), MARGIN = 1, FUN = mean))
  MAX_avg_logFC <- apply(do.call(cbind, lapply(de_results, FUN = function(x){x$avg_logFC})), MARGIN = 1, FUN = max)
  MIN_avg_logFC <- apply(do.call(cbind, lapply(de_results, FUN = function(x){x$avg_logFC})), MARGIN = 1, FUN = min)
  
  de_results <- lapply(names(de_results), FUN = function(x){names(de_results[[x]]) <- paste0(names(de_results[[x]]), "_", x);
  return(de_results[[x]])})
  de_results <- cbind(MEAN_avg_logFC, MAX_avg_logFC, MIN_avg_logFC, MAX_p_val, MAX_p_val_adj_BH, MAX_p_val_adj_Bonf,
                      do.call(cbind, de_results))
  de_results <- de_results[order(de_results$MEAN_avg_logFC, decreasing = TRUE),]
  return(de_results)
}

#define cell barcodes for each of the wing disc domains (notum, hinge, pouch, anterior, posterior)
#anterior and posterior must be mutually exclusive with one another, and notum, hinge, and pouch must be mutually exclusive with each other
#e.g., the same barcode cannot be assigned to both notum and hinge (but a barcode can be assigned to both notum and anterior)

#anterior_cells <- vector_of_cell_barcodes_classified_as_anterior
#posterior_cells <- vector_of_cell_barcodes_classified_as_posterior
#notum_cells <- vector_of_cell_barcodes_classified_as_notum
#hinge_cells <- vector_of_cell_barcodes_classified_as_hinge
#pouch_cells <- vector_of_cell_barcodes_classified_as_pouch


#perform pairwise-differential expression for anterior v. posterior, notum v. hinge, notum v. pouch, and hinge v. pouch
#replace "orig.ident" with the meta.data name in the seurat_object that corresponds to batch identities 

de_ant_v_post <- PairwiseFindMarkers(seurat_object,
                                     cells.1 = anterior_cells,
                                     cells.2 = posterior_cells,
                                     split.by = "orig.ident", features = VariableFeatures(seurat_object),
                                     only_compare_same_split_idents = TRUE, perform_de_direction_check = TRUE,
                                     logfc.threshold = -Inf, min.diff.pct = -Inf, min.pct = -Inf, min.cells.group = 0)

de_notum_v_pouch <- PairwiseFindMarkers(seurat_object,
                                        cells.1 = names(seurat_object$Epi_Domains[seurat_object$Epi_Domains %in% c("NOTUM")]),
                                        cells.2 = names(seurat_object$Epi_Domains[seurat_object$Epi_Domains %in% c("POUCH")]),
                                        split.by = "orig.ident", features = VariableFeatures(seurat_object),
                                        only_compare_same_split_idents = TRUE, perform_de_direction_check = TRUE,
                                        logfc.threshold = -Inf, min.diff.pct = -Inf, min.pct = -Inf, min.cells.group = 0)

de_hinge_v_pouch <- PairwiseFindMarkers(seurat_object,
                                        cells.1 = names(seurat_object$Epi_Domains[seurat_object$Epi_Domains %in% c("HINGE")]),
                                        cells.2 = names(seurat_object$Epi_Domains[seurat_object$Epi_Domains %in% c("POUCH")]),
                                        split.by = "orig.ident", features = VariableFeatures(seurat_object),
                                        only_compare_same_split_idents = TRUE, perform_de_direction_check = TRUE,
                                        logfc.threshold = -Inf, min.diff.pct = -Inf, min.pct = -Inf, min.cells.group = 0)

de_notum_v_hinge <- PairwiseFindMarkers(seurat_object,
                                        cells.1 = names(seurat_object$Epi_Domains[seurat_object$Epi_Domains %in% c("NOTUM")]),
                                        cells.2 = names(seurat_object$Epi_Domains[seurat_object$Epi_Domains %in% c("HINGE")]),
                                        split.by = "orig.ident", features = VariableFeatures(seurat_object),
                                        only_compare_same_split_idents = TRUE, perform_de_direction_check = TRUE,
                                        logfc.threshold = -Inf, min.diff.pct = -Inf, min.pct = -Inf, min.cells.group = 0)


#combine all DE results into a single list
de_master_df <- list(ant_v_post = de_ant_v_post, notum_v_hinge = de_notum_v_hinge,
                     notum_v_pouch = de_notum_v_pouch, hinge_v_pouch = de_hinge_v_pouch)


#flatten this list to be format appropriate for ggplot
plot_dataframe <- data.frame(variable = factor(do.call(c, lapply(names(de_master_df), FUN = function(x){rep(x, times = NROW(de_master_df[[x]]))})),
                                               levels = rev(names(de_master_df))),
                             logFC = do.call(c, lapply(de_master_df, FUN = function(x){return(x[["MEAN_avg_logFC"]])})),
                             max_pval = do.call(c, lapply(de_master_df, FUN = function(x){return(x[["MAX_p_val"]])})),
                             max_pval_BH = do.call(c, lapply(de_master_df, FUN = function(x){return(x[["MAX_p_val_adj_BH"]])})))


#ggplot code for plotting the ECDF of log-fold changes from DE results
ggplot(plot_dataframe, aes(x = abs(logFC))) +
  stat_ecdf(aes(group = variable, colour = variable), size = 2.5) +
  labs(x = "Abs. Log(e) F.C.", y = "Empirical CDF", color = "Differential Expression") +
  xlim(0, 0.5) + theme_bw() + scale_color_manual(values=c("purple", "green", "blue", "red")) +
  theme(legend.title = element_text(size = 35, face = "bold"), legend.text = element_text(size = 28, face = "bold"), legend.key.size = unit(0.8, "cm"),
        axis.title = element_text(size = 35, face = "bold"), axis.text = element_text(size = 28), legend.position = "none",
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
# ggsave(filename = "./LogFC_ECDF.png", device = "png", height = 8, width = 11, dpi = 500)


#ggplot code for plotting the ECDF of FDR values from DE results
ggplot(plot_dataframe, aes(x = -log10(max_pval_BH))) +
  stat_ecdf(aes(group = variable, colour = variable), size = 2.5) +
  labs(x = "-Log10(FDR)", y = "Empirical CDF", color = "Differential Expression") +
  xlim(0, 30) + theme_bw() + scale_color_manual(values=c("purple", "green", "blue", "red")) +
  theme(legend.title = element_text(size = 35, face = "bold"), legend.text = element_text(size = 28, face = "bold"), legend.key.size = unit(0.8, "cm"),
        axis.title = element_text(size = 35, face = "bold"), axis.text = element_text(size = 28), legend.position = "none",
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
# ggsave(filename = "./pVal_ECDF.png", device = "png", height = 8, width = 11, dpi = 500)


#ggplot code for plotting the density of log-fold changes from DE results
ggplot(plot_dataframe, aes(x=logFC, color=variable)) +
  labs(x = "Log(e) F.C.", y = "Density", fill = "Differential Expression") +
  geom_line(size = 2.5, stat = "density") +
  coord_cartesian(ylim=c(0, 7), xlim = c(-0.5, 0.5)) + theme_classic() +
  theme(legend.title = element_text(size = 28, face = "bold"), legend.text = element_text(size = 35, face = "bold"), legend.key.size = unit(1.5, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 35, face = "bold"), axis.text = element_text(size = 28),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  scale_color_manual(values=c("purple", "green", "blue", "red"))
# ggsave(filename = "./LogFC_Density.png", device = "png", height = 8, width = 11, dpi = 500)