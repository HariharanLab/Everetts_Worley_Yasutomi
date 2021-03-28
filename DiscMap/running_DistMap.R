# R code is for generating a virtual disc model using Distmap (Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 

# will need a geometry file (supplemental file 4)
# reference expression patterns (supplemental file 5) 
# and single-cell RNA seq data (we used imputed values from ScVI for mapping the cells to the locations). 
# The cells will be mapped to the reference patterns, and then can be used to generate virtual in situ or predicted gene expression patterns. 


# Import libraries

#Distmap ( Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 
library(DistMap)


geometryfile = "~/supplemtnal_file_4.csv"
discExp = "~/supplemtnal_file_5.csv"
sc_object = "~/object_imputed"

imputedData = as.matrix(sc_object@assays$scVIxNB_Imputed@data)

# Running DistMap 

dmDiscImputed = new("DistMap",
         data = imputedData,
         insitu.matrix = discExp,
         geometry= geometryfile )

dmDiscImputed <- binarizeSingleCellData(object = dmDiscImputed, quantiles = seq(0, 0.5, 0.01))

dmDiscImputed <- mapCells(object = dmDiscImputed )


# to calculate virtual in situ :  

# gene of interest 
gene = gene 

VISH <- computeVISH( dmDiscImputed, gene, threshold = 0)



