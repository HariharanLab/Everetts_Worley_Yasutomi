## Here is code for visulizing the predicted expression patterns from the Virtual Wing disc (from the DiscMap model) 


# import libraries 
library(DistMap)
library(ggplot2)
library(gridExtra)
library(viridis)
viridis = viridis(n = 16, begin = 0, end = 1, option = "D")

# need to have DistMap model loaded and the geometry from supplemental_file_4:  
dmDiscImputed = readRDS( file = paste0("dmDiscImputed", ".rds")) # load DistMap model
geometryfile = as.matrix(read.csv("./supplemental_file_4.csv", row.names = 1, header = TRUE, stringsAsFactor = FALSE))

#Functions for visulization
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
pickcolor <- function(Vish_data, lowlim, uplim){
  range01(pmax(pmin(range01(Vish_data), uplim), lowlim))
}


# gene of interest 
gene = 'pyr'
VISH <- computeVISH( dmDiscImputed  , gene, threshold = 0)

# This changes the min/max of the in situ.  
VISHcolor = pickcolor(VISH, 0.5, 1)

#Create a data fram of geometry and in situ data togeter 
toGraph <- as.data.frame(cbind(geometryfile, VISHcolor))
Nzlayers <- unique(toGraph$z)
plotsHere <- list()
counter = 1 
for(n in Nzlayers){
  tosave <- toGraph[toGraph$z == n,]
  p1 <- ggplot(data = tosave, mapping = aes( x, y, color= VISHcolor )) + xlim( min(toGraph$x), max(toGraph$x)) + ylim( min(toGraph$y), max(toGraph$y)) + theme_void() + geom_point(size= 5, show.legend = FALSE , alpha = 1) + scale_colour_gradientn( colours = viridis, limits=c( min(toGraph$VISHcolor) , max(toGraph$VISHcolor)))
  plotsHere[[counter]] <-  p1
  counter <- counter + 1
}

# This is to plot the three layers (AMPs, Epithelium, PE) next to each other 
grid.arrange( plotsHere[[1]] , plotsHere[[2]], plotsHere[[3]], ncol = 3, top= paste(gene))
