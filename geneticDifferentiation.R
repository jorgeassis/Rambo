## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
##
##                                      #####
##                                ####  #####
##                                ####       
##          ####                         
##         ##################             
##           ##################           
##       #######################
##   ##################################   
##  #######################################
##  ######################################
##  ###################################### 
##  ####################################
##  ##################################     
##  ####################                   
##  ###################                    
##  ##################                     
##  #################                      
##  ###############                                     
##      
##  theMarineDataScientist
##
##  github.com/jorgeassis
##  medium.com/themarinedatascientist
##  medium.com/@jorgemfa
##
## -------------------------------------------------------------------------------
##
##  rambo 2.0
##  R Pipelines for Genetic Analyses
##
## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

library(adegenet) 
library(PopGenReport)
library(reshape2)
library(hierfstat)
library(raster)
library(mmod)

## --------------------------------------

mainDataFile <- "/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Genetic/All Genetic Data 6L Genepop Clean Final Genetix.gtx"
data <- read.genetix(mainDataFile) 
D <- pairwise_D(data)
D <- as.matrix(D)

Fst <- pairwise.fst(data)
Fst <- pairwise.fst(data)
Fst <- as.matrix(Fst)

clusteringVector <- c( rep(1,72) , rep(2,43) )
clusteringVector <- c(rep(1,7),2,2,1,rep(2,9),rep(3,18),rep(4,35),5,5,rep(6,7),rep(7,23),rep(6,5),7,7,7,6,6,7)

differentiation <- data.frame(matrix(ncol=length(unique(clusteringVector)),nrow=length(unique(clusteringVector))))
  
for( i in 1:ncol(differentiation)) {
  
  for( j in 1:ncol(differentiation)) {
    
    differentiation[i,j] <- mean(D[  which( clusteringVector == i) , which( clusteringVector == j) ])

  }
}
  
colnames(differentiation) <- 1:ncol(differentiation)
rownames(differentiation) <- 1:ncol(differentiation)
differentiation <- differentiation[-5,-5]

dendogram <- hclust(as.dist(differentiation), method = "complete", members = NULL)
dev.off()
plot(dendogram,hang = -1, cex = 0.8)
plot(as.phylo(dendogram), cex = 0.6, label.offset = 0.5)
plot(as.phylo(dendogram), type = "unrooted", cex = 0.6,no.margin = TRUE)

diag(differentiation) <- NA
differentiation[lower.tri(differentiation)] <- NA
differentiation

library(RColorBrewer)
pal <- colorRampPalette(c("#9ad1ee","#b10303"))
breakpoints <- seq(min(differentiation,na.rm=T),max(differentiation,na.rm=T),length.out = 10)
plot(raster(as.matrix(differentiation)),breaks=breakpoints,col=pal(10))

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
