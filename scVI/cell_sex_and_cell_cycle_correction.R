#########################################################################
#### CLASSIFY CELLS AS ORIGINATING FROM EITHER MALE OR FEMALE LARVAE ####
#########################################################################

#define the genes used for identifying male and female cells
#we suggest using lncRNA:roX1 and lncRNA:roX2, as they are both expressed at relatively high levels
sex_specific_genes <- c("lncRNA:roX1", "lncRNA:roX2")

#extract the log-normalized data from the Seurat object
data_matrix <- as.matrix(GetAssayData(seurat_object, assay = "RNA", slot = "data"))

#initialize two empty lists that will store the expression cutoffs for the sex
#genes, and the classification of cells as male / female
sex_gene_cutoff <- vector(mode = "list", length = length(sex_specific_genes))
sex_threshold_cells <- vector(mode = "list", length = length(sex_specific_genes))
names(sex_gene_cutoff) <- sex_specific_genes
names(sex_threshold_cells) <- sex_specific_genes

for (gene in sex_specific_genes) {
  data_density <- density(data_matrix[gene,], adjust = 2)
  
  #NOTE: We expect the data for sex genes to be bimodal, and use the minimum
  #between peaks as a cutoff for determining male and female cells. This usually
  #corresponds to the second extrema in the density (from left to right).
  #However, noise and drop-outs in the data may make this difficult to
  #calculate, so a manual cutoff may be needed instead.
  
  #find the second extrema in the density, and use this as a cutoff for cells
  sex_gene_cutoff[[gene]] <- data_density$x[which(diff(diff(data_density$y)>0)!=0)][2]
  sex_threshold_cells[[gene]] <- colnames(data_matrix)[data_matrix[gene,] > sex_gene_cutoff[[gene]]]
  
  #plot a histogram of the gene expression data, with a vertical line indicating the cutoff
  hist(data_matrix[gene,], breaks = 50, probability = TRUE,
       main = paste("Density of", gene, "expression"), xlab = paste("Log Normalized Counts of", gene))
  lines(data_density, col = "blue", lwd = 2)
  abline(v = sex_gene_cutoff[[gene]], col = "red", lwd = 2)
}

#print out statistics for male and female cells
for (gene in names(sex_threshold_cells)) {
  cat("Number of cells with high", gene, "expression: ", length(sex_threshold_cells[[gene]]), "\n")
  cat("Number of cells with low", gene, "expression: ", NCOL(seurat_object) - length(sex_threshold_cells[[gene]]), "\n\n")
}

male_cells <- unique(do.call(what = c, sex_threshold_cells))
female_cells <- colnames(seurat_object)[!(colnames(seurat_object) %in% male_cells)]

cat("\nNumber of predicted MALE cells: ", length(male_cells))
cat("\nNumber of predicted FEMALE cells: ", length(female_cells))
cat("\n\nRatio of MALE to FEMALE cells: ", length(male_cells)/length(female_cells))
cat("\nFraction of cells that are MALE: ", length(male_cells)/(length(male_cells)+length(female_cells)))
cat("\n\nRatio of FEMALE to MALE cells: ", length(female_cells)/length(male_cells))
cat("\nFraction of cells that are FEMALE: ", length(female_cells)/(length(male_cells)+length(female_cells)))

#this classification can be corrected for similar to batch correction


############################################################################
#### EXAMINE THE CORRELATION OF LATENT DIMENSIONS WITH CELL CYCLE GENES ####
############################################################################

#create a heatmap showing the (abs.) correlation of cell cycle genes vs. latent dimensions, as shown in Figure 1 - figure supplement 5D

#load the cell-cycle genes into R (you can supply your own list as a .csv, if preferable)
cell_cycle_genes <- read.delim(file = "E:/scRNAseq/FlyBase_Batch_Downloads/cell_cycle_genes.csv",
                               stringsAsFactors = FALSE, header = FALSE)[,1]

#subset the cell-cycle genes to only those that are variably expressed
cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes %in% VariableFeatures(seurat_object, assay = paste0(reduction_key, "Imputed"))]
cell_cycle_genes <- sort(cell_cycle_genes)


corr_methods <- c("Pearson", "Spearman")
for (md in corr_methods) {
  #In this code, we have stored the scVI latent space as a reduction in the Seurat object (similar to storing PCA, tSNE, UMAP, etc.).
  
  #calculate the correlation between each scVI latent dimension and the cell-cycle genes provided
  latent_cell_cycle_correlation <- cor(seurat_object@reductions$scVI@cell.embeddings, 
                                       t(as.matrix(GetAssayData(seurat_object, assay = "RNA", slot = "data")[cell_cycle_genes,])),
                                       method = tolower(md))
  rownames(latent_cell_cycle_correlation) <- paste0("scVI_", 1:n_latent)
  
  #OPTIONAL: reorder the rows of the correlation matrix via hierarchical clustering
  #this will put genes with similar latent space correlations next to each other, may look more visually appealing
  hclust_results <- hclust(dist.matrix(latent_cell_cycle_correlation, byrow = FALSE, method = "euclidean", as.dist = TRUE), method = "complete")
  latent_cell_cycle_correlation <- latent_cell_cycle_correlation[, hclust_results$order]
  
  #OPTIONAL: for better visualization, cap the correlation between 0 and 0.5
  corr_limits <- c(0, 0.5)
  
  #flatten the correlation matrix data
  #we also take the absolute value of correlation values, as we only care about
  #the correlation magnitudes (a negative correlation with cell cycle genes is
  #just as troublesome as a positive correlation with cell cycle genes)
  latent_cell_cycle_corr_melted <- melt(abs(latent_cell_cycle_correlation))
  
  #for correlation values that fall outside the correlation limits, provide text values in the appropriate cells
  text_values <- round(latent_cell_cycle_corr_melted$value, digits = 2)
  text_values[text_values >= corr_limits[1] & text_values <= corr_limits[2]] <- ""
  
  #ggplot code for the correlation heatmap with text values in the appropriate cells
  p <- ggplot(data = latent_cell_cycle_corr_melted,  aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + coord_fixed(ratio = 0.5) +
    geom_text(aes(label = text_values), size = 3.5) +
    xlab("scVI Latent Dimensions") + ylab("Variable Cell Cycle Genes") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(corr_limits),
                         limit = corr_limits, space = "Lab", name=paste(md, "Correlation\n(Absolute Value)"), oob = scales::squish) +
    theme(axis.title.x = element_text(face = "bold", size = 22),
          axis.title.y = element_text(face = "bold", size = 22),
          axis.text.x = element_text(angle = 90, vjust = 0.3, face = "bold", size = 12),
          axis.text.y = element_text(face = "bold.italic", vjust = 0.3, size = 11),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(face = "bold", size = 12))
  print(p)
  # ggsave(filename = paste0("./LatentSpace-CellCycle-", md,"CorrHeatmap.png"), device = "png", width = 10, height = 9, dpi = 500)
}

#we often find that only one latent dimension displays high correlation to cell cycle genes,
#and exclude this dimensions from our downstream analysis