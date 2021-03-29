# R code is for generating a virtual disc model using Distmap (Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 

# will need a geometry file (supplemental file 4)
# reference expression patterns (supplemental file 5) 
# and single-cell RNA seq data (we used imputed values from scVI for mapping the cells to the locations). 
# The cells will be mapped to the reference patterns, and then can be used to generate virtual in situ or predicted gene expression patterns. 


# Import libraries

#DistMap (Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 
library(DistMap)

#load in the geometry coordinates (supplemental_file_4.csv), the coordinates for reference genes (supplemental_file_5.csv), and
#the imputed expression data from scVI.
geometryfile = as.matrix(read.csv("./supplemental_file_4.csv", row.names = 1, stringsAsFactor = FALSE))
discExp = as.matrix(read.csv("./supplemental_file_5.csv", header = TRUE, row.names = 1, stringsAsFactor = FALSE))
imputedData = as.matrix(read.csv(path_to_imputed_data, stringsAsFactor = FALSE))

# Builing DistMap model 

dmDiscImputed = new("DistMap",
         data = imputedData,
         insitu.matrix = discExp,
         geometry= geometryfile)
dmDiscImputed <- binarizeSingleCellData(object = dmDiscImputed, quantiles = seq(0, 0.5, 0.01))
dmDiscImputed <- mapCells(object = dmDiscImputed)

# to calculate virtual in situ:

# gene of interest 
gene = gene 
VISH <- computeVISH( dmDiscImputed, gene, threshold = 0)
