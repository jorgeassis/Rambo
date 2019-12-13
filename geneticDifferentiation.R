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
library(mmod)

## --------------------------------------

mainDataFile <- "/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Genetic/All Genetic Data 6L Genepop Clean Final Genetix.gtx"
data <- read.genetix(mainDataFile) 
D <- pairwise_D(data)
D <- as.matrix(D)

clusteringVector <- c( rep(1,72) , rep(2,43) )

differentiation <- data.frame(matrix(ncol=length(unique(clusteringVector)),nrow=length(unique(clusteringVector))))
  
for( i in 1:ncol(differentiation)) {
  
  for( j in 1:ncol(differentiation)) {
    
    differentiation[i,j] <- mean(D[  which( clusteringVector == i) , which( clusteringVector == j) ])
    
  }
}
  
differentiation

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
