# R code is for generating a virtual disc model using Distmap (Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 

# will need a geometry file (supplemental file 4) : this will have x,y,z coordinate for the location of the reference gene expression patterns 
#                 x   y  z
#        1    -11.5  38  0
#        2    -10.5  38  0
#        3     -9.5  38  0

# reference expression patterns (supplemental file 5) 
#        For the following genes: Ance, ap, brk, caup, ci, ct, Dad, Dll, Doc1, dpp, dve, en, eya, eyg, Fas3, fj, hh, Him, hth, htl, mirr, Nrt, 
#        nub, pdm2, pnr, ptc, rho, rn, salm, Sox15, SPARC, ths, tsh, tup, twi, Ubx, vg, vkg, vvl, wg, Wnt6, zfh2

# and single-cell RNA seq data (we used imputed values from scVI for mapping the cells to the locations). 
# The cells will be mapped to the reference patterns, and then can be used to generate virtual in situ or predicted gene expression patterns. 


# Import libraries

#DistMap (Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 
library(DistMap)

#load in the geometry coordinates file (supplemental_file_4.csv), the coordinates for reference gene expression patterns (supplemental_file_5.csv), and
#the imputed expression data from scVI.
geometryfile = as.matrix(read.csv("./supplemental_file_4.csv", row.names = 1, header = TRUE, stringsAsFactor = FALSE))
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

# this virtual in situ can then be plotted on the same geometry as povided above (see plottingDiscmap.R) 

saveRDS( dmDiscImputed , file = paste0("dmDiscImputed", ".rds")) # Save DistMap model
