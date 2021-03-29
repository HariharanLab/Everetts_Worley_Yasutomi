
gene = "pyr"

VISH <- computeVISH( dmDiscMyoImputed  , gene, threshold = 0)

VISHcolor = pickcolor(VISH, 0.5, 1)
toGraph <- as.data.frame(cbind(geometry, VISHcolor))

Nzlayers <- unique(toGraph$z)
plotsHere <- list()
counter = 1 

for(n in Nzlayers){
  tosave <- toGraph[toGraph$z == n,]
  p1 <- ggplot(data = tosave, mapping = aes( x, y, color= VISHcolor )) + xlim( min(toGraph$x), max(toGraph$x)) + ylim( min(toGraph$y), max(toGraph$y)) + theme_void() + geom_point(size= 5, show.legend = FALSE , alpha = 1) + scale_colour_gradientn( colours = viridis, limits=c( min(toGraph$VISHcolor) , max(toGraph$VISHcolor)))
  plotsHere[[counter]] <-  p1
  counter <- counter + 1
}

p <- grid.arrange( plotsHere[[1]] , plotsHere[[2]], plotsHere[[3]], ncol = 3, top= paste(gene))
