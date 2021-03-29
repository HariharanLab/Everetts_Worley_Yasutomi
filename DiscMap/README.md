# Virtual three-layered wing  disc from single-cell data mapped to a reference map


![alt text](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/DiscMap/discmap_image.jpg?raw=true)


* Distmap ( Karaiskos et al., 2017 - https://github.com/rajewsky-lab/distmap) 


We assembled reference gene expression patterns from a number of sources:
* [Held Jr, 2002](https://www.sdbonline.org/sites/fly/lewheld/00idheld.htm)
* [Butler et al., 2003](https://dev.biologists.org/content/130/4/659.long)
* and based our starting geometry on the disc proper from images in [Bageritz et al., 2019](https://www.nature.com/articles/s41592-019-0492-x?proof=t)

* The images were processed in Adobe Photoshop and assembled in R with EBimage ([Pau et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20338898/)) to generate binarized gene expression reference for the AMPs, disc proper, and peripodial epithelium. The geometry of the three-layered model is provided in [Supplementary file 4](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/DiscMap/Supplementary%20file%204.csv) and the binarized reference gene expression patterns are provided in [Supplementary file 5](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/DiscMap/Supplementary%20file%205.csv). 

* We used DistMap ([Karaiskos et al., 2017](https://science.sciencemag.org/content/358/6360/194)) to statistically map single cells back to the reference. With this virtual wing disc model we used DistMap to calculate a ‘virtual _in situ_’ or a prediction of gene expression patterns. This is based on the detected gene expression with the single-cell data and the mapping location to calculate relative expression values for our model. 

* We mapped the AMP and epithelial cells separately, as this improved how well the model predicted genes with known expression patterns. In addition, we used the scVI imputed gene expression values when mapping the cells to the reference [(see section on scVI)](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/tree/master/scVI).

>
> We provide example R code for generating the virtual disc model with DistMap [(running_DistMap.R)](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/DiscMap/running_DistMap.R)
>
> and for generating the plots of the three-layered disc from the virtual _in situ_ data [(plottingDiscmap.R)](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/DiscMap/plottingDiscmap.R)

![alt text](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/DiscMap/Three_layerDisc.jpg)
