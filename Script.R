## -----------------------------------------------------
## Rambo V1.0
## Assis et al., 2018
## -----------------------------------------------------
## -----------------------------------------------------

## -----------------------------------------------------
## -----------------------------------------------------

Rambo <- function (main.data.file, missing.data, ncode , replace, resample.number.auto, resample.number,
                   discard.pops, number.iteractions, alfa.test, clustering.vector, savefile)
  
{
  ## -----------------------------------------------------
  ## Installing and Reading Dependences
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  ## ----------------
  
  installing <- FALSE
  
  if( !"adegenet" %in% as.character(installed.packages()[,"Package"]) ) {
    
    cat( "\r" , "\r" , " Installing Dependences...", "\r", "\r")
    suppressMessages(suppressWarnings(
      try(install.packages("adegenet" , verbose=FALSE , quiet = TRUE) , silent = TRUE) ))
    
    installing <- TRUE
    
  }
  
  suppressMessages(suppressWarnings(  try( library(adegenet) , silent = TRUE) ))
  
  if( !"PopGenReport" %in% as.character(installed.packages()[,"Package"]) ) {
    
    if( ! installing ) {
      
      cat( "\r" , "\r" , " Installing Dependences...", "\r", "\r")
      
    }
    suppressMessages(suppressWarnings(
      try(install.packages("PopGenReport" , verbose=FALSE , quiet = TRUE) , silent = TRUE) ))
  }
  
  suppressMessages(suppressWarnings(  try( library(PopGenReport) , silent = TRUE) ))
  
  ## -----------------------------------------------------
  ## Reading Data
  
  # A genepop object
  
  if( class(main.data.file) == "genind" ) {
    
    data <- main.data.file
    
  }
  
  # A genepop file
  
  if( class(main.data.file) == "character" ) {
    
      suppressMessages(suppressWarnings( try( data <- read.genetix(main.data.file) , silent = FALSE) ))
   
    
  }
  
  if( ! exists("data") ) {
    
    stop("Error or Uncoded situation [!]\nGenetix Data is unreadable or the path is uncorrect")
    
  }
  
  ## -----------------------------------------------------
  ## Main Objets and Defenitions
  
  if( is.null(resample.number.auto) | ! exists("resample.number.auto") ) { resample.number.auto <- TRUE }
  if( is.null(clustering.vector) | ! exists("clustering.vector") ) { clustering.vector <- NULL }
  if( is.null(resample.number) | ! exists("resample.number") ) { resample.number <- NULL }
  if( is.null(savefile) | ! exists("savefile") ) { savefile <- NULL }
  if( is.null(replace) | ! exists("replace") ) { replace <- FALSE }
  if( is.null(resample.number.auto) | ! exists("resample.number.auto") ) { resample.number.auto <- TRUE }
  if( is.null(discard.pops) | ! exists("discard.pops") ) { discard.pops <- NULL }
  if( is.null(number.iteractions) | ! exists("number.iteractions") ) { number.iteractions <- 999 }
  if( is.null(alfa.test) | ! exists("alfa.test") ) { alfa.test <- 0.05 }
  
  loci.names <- as.character(unique(data$loc.fac))
  loci.names.individuals <- data$loc.fac
  loci.n <- length(loci.names)
  
  cat("\n")
  cat("\n")
  cat("\n")
  cat(" ------------------------------------------------------------------------------------------------")
  cat("\n")
  cat(" Rambo . Random reductions for genetic data")
  cat("\n")
  cat("\n")
  cat(" Data file: ", main.data.file)
  cat("\n")
  cat(" Loci: ", loci.n)
  cat("\n")
  cat(" Ploidy: ",unique(data$ploidy))
  
  ## -----------------------------------------------------
  ## Discard Populations
  
  if( exists("discard.pops") & ! is.null(discard.pops) ) {
    
    data <- data[ ! data@pop %in% unique(data@pop)[discard.pops]  , ]
    cat("\n")
    cat(" Removed populations: ",discard.pops)
    
  }
  
  cat("\n")
  cat(" Populations: ",length(unique(data@pop)))
  
  cat("\n")
  cat("\n")
  cat(" ------------------------------------------------------------------------------------------------")
  cat("\n")
  
  ## -----------------------------------------------------
  ## Other Objets and Defenitions
  
  population.names <- as.character(unique(data@pop))
  population.names.results <- as.character(substr(population.names, 1, 10))
  population.n <- length(unique(data@pop))
  
  individual.n.per.pop <- vector()
  
  for (f in 1:population.n){
    
    individual.n.per.pop <- c(individual.n.per.pop,nrow(data[ data@pop == population.names[f] , ]@tab))
    
  }
  
  ## -----------------------------------------------------
  ## Clustering
  
  prior.variable.clustering <- TRUE
  
  if( is.null(clustering.vector) | ! exists("clustering.vector") ) {
    
    clustering.vector <- 1:population.n
    cat("\n")
    cat(" Clustering: FALSE")
    cat("\n")
    cat("\n")
    
    prior.variable.clustering <- FALSE
    
  }
  
  if ( prior.variable.clustering & ! is.null(clustering.vector) | prior.variable.clustering & exists("clustering.vector") ) {
    
    if( length(clustering.vector) != population.n) { stop("Error or Uncoded situation [!]\nThe number of populations does not match the length of the clustering vector.")   }
    
    clustering.array <- vector()
    
    for (i in 1:population.n ) {
      
      clustering.array <- c(clustering.array,rep(paste("POP.",as.character(clustering.vector[i]),sep=""),individual.n.per.pop[i]))
      
    }
    
    data@pop <- as.factor(clustering.array)

    cat("\n")
    cat("\n")
    
  }
  
  population.names <- unique( data@pop )
  individual.names <- data@pop
  population.n <- length(population.names)
  individual.n <- length(individual.names)
  
  if ( resample.number.auto ) {
    
    samples <- vector()
    
    for (i in 1:population.n) {
      
      samples <- c(samples, sum( !is.na(which(individual.names %in% population.names[i]))))
    }
    
    resamp <- samples[which.min(samples)]
    
  }
  
  if ( ! resample.number.auto | is.null(resample.number.auto) ) {
    
    resamp <- resample.number
    
  }
  
  ## -----------------------------------------------------
  ## Main Random Matrices to produce Results
  
  results.unique <- matrix(NA, ncol = number.iteractions, nrow = population.n)
  results.unique.others <- matrix(NA, ncol = number.iteractions, nrow = population.n)
  
  results.richness <- matrix(NA, ncol = number.iteractions, nrow = population.n)
  results.richness.others <- matrix(NA, ncol = number.iteractions, nrow = population.n)
  
  results.he <- matrix(NA, ncol = number.iteractions, nrow = population.n)
  results.he.others <- matrix(NA, ncol = number.iteractions, nrow = population.n)
  
  for (int in 1:number.iteractions) {
    
    progress <- round((int/number.iteractions) * 100, digits = 0)
    
    cat("     \r")
    cat(paste(" Rambo in progress: ", progress, " % | ", number.iteractions,
              " randomizations for ", population.n, " populations with ",
              resamp, " individuals                                       ",
              sep = ""), "\r")
    flush.console()
    
    resample.array <- vector()
    resample.array.others <- vector()
    
    for( i in 1:population.n) {
      
      if( sum(individual.names == population.names[i]) < resamp ) {  replace.pop <- TRUE  } else { replace.pop <- replace  }
      resample.array <- c(resample.array,sample(which(individual.names == population.names[i]), resamp, replace = replace.pop))
      resample.array.others <- c(resample.array.others,sample(which(individual.names != population.names[i]), resamp, replace = replace.pop))
    }
    
    data.resamp <- data[ resample.array , ]
    data.resamp.otehrs <- data[ resample.array.others , ]
    
    unique <- vector()
    unique.others <- vector()
    alleles <- vector()
    alleles.other <- vector()
    
    # Allele Richness & Unique Alleles
    
    for (j in 1:population.n ) {
      
      all.alleles.perlocus.j <- apply( data.resamp[ data.resamp@pop == population.names[j] , ]@tab , 2 , sum , na.rm=T)
      all.alleles.perlocus.j <- data.frame(all.alleles.perlocus.j)
      all.alleles.perlocus.j <- all.alleles.perlocus.j[substrRight(rownames(all.alleles.perlocus.j), 3) != missing.data,]
      all.alleles.perlocus.j[all.alleles.perlocus.j > 0] <- 1
      
      all.alleles.perlocus.others.j <- apply( data.resamp[ data.resamp@pop != population.names[j] , ]@tab , 2 , sum, na.rm=T)
      all.alleles.perlocus.others.j <- data.frame(all.alleles.perlocus.others.j)
      all.alleles.perlocus.others.j <- all.alleles.perlocus.others.j[substrRight(rownames(all.alleles.perlocus.others.j), 3) != missing.data,]
      all.alleles.perlocus.others.j[all.alleles.perlocus.others.j > 0] <- 1
      
      merged <- rbind(all.alleles.perlocus.j , all.alleles.perlocus.others.j )
      
      unique <- c(unique , length(which(merged[1,] == 1 & merged[2,] == 0)) )
      unique.others <- c(unique.others , length(which(merged[1,] == 0 & merged[2,] == 1)) )
      
      alleles <- c( alleles , (sum(all.alleles.perlocus.j) ) / loci.n )
      alleles.other <- c(alleles.other , (sum(all.alleles.perlocus.others.j > 0) ) / loci.n )
      
    }
    
    he <- Hs(data.resamp)
    he.other <- Hs(data.resamp.otehrs)
    
    # Place results
    
    results.unique[, int] <- unique
    results.unique.others[, int] <- unique.others
    results.richness[, int] <- alleles
    results.richness.others[, int] <- alleles.other
    results.he[, int] <- he
    
    results.he.others[ which(unique(data@pop) %in% names(he)[which(names(he) %in% names(he.other))]) , int] <- he.other
    
  }
  
  results.significance.richness <- numeric(population.n)
  results.significance.unique <- numeric(population.n)
  results.significance.he <- numeric(population.n)
  
  for (j in 1:population.n ) {
    
    results.significance.richness[j] <- signif(wilcox.test(results.richness[j, ] , results.richness.others[j, ], alternative = "greater")$p.value,3)
    results.significance.unique[j] <- signif(wilcox.test(results.unique[j, ] , results.unique.others[j, ], alternative = "greater")$p.value,3)
    results.significance.he[j] <- signif(wilcox.test(results.he[j, ] , results.he.others[j, ], alternative = "greater")$p.value,3)
    
  }
  
  results.significance.richness[results.significance.richness > alfa.test] <- "NS"
  results.significance.unique[results.significance.unique > alfa.test] <- "NS"
  results.significance.he[results.significance.he > alfa.test] <- "NS"
  
  
  if( length(unique(clustering.vector)) == length(individual.n.per.pop) ) {
    
    main.results.data.frame.file <-  data.frame( Population = ( population.names.results ) ,
                                                 n = individual.n.per.pop ,
                                                 A = round(apply(results.richness, 1, mean), 2) ,
                                                 SD.A = round(apply(results.richness, 1, sd), 2) ,
                                                 Signif.A = results.significance.richness ,
                                                 PA = round(apply(results.unique, 1, mean), 2) ,
                                                 SD.PA = round(apply(results.unique, 1, sd), 2) ,
                                                 Signif.PA = results.significance.unique ,
                                                 He = round(apply(results.he, 1, mean), 2),
                                                 SD.He = round(apply(results.he, 1, sd), 2),
                                                 Signif.He = results.significance.he
                                                 
    )
    
    main.results.data.frame <-  data.frame( Population = ( population.names.results ) ,
                                            . = rep("        ", length(individual.n.per.pop) ) ,
                                            n = individual.n.per.pop ,
                                            . = rep("        ", length(individual.n.per.pop) ) ,
                                            A = round(apply(results.richness, 1, mean), 2) ,
                                            SD.A = round(apply(results.richness, 1, sd), 2) ,
                                            Signif.A = results.significance.richness ,
                                            .. = rep("        ", length(individual.n.per.pop) ) ,
                                            PA = round(apply(results.unique, 1, mean), 2) ,
                                            SD.PA = round(apply(results.unique, 1, sd), 2) ,
                                            Signif.PA = results.significance.unique ,
                                            ... = rep("        ", length(individual.n.per.pop) ) ,
                                            He = round(apply(results.he, 1, mean), 2),
                                            SD.He = round(apply(results.he, 1, sd), 2),
                                            Signif.He = results.significance.he
                                            
    ) }
  
  if( length(unique(clustering.vector)) != length(individual.n.per.pop) ) {
    
    main.results.data.frame.file <-  data.frame( Population = 1:length(unique(clustering.vector)) ,
                                                 A = round(apply(results.richness, 1, mean), 2) ,
                                                 SD.A = round(apply(results.richness, 1, sd), 2) ,
                                                 Signif.A = results.significance.richness ,
                                                 PA = round(apply(results.unique, 1, mean), 2) ,
                                                 SD.PA = round(apply(results.unique, 1, sd), 2) ,
                                                 Signif.PA = results.significance.unique ,
                                                 He = round(apply(results.he, 1, mean), 2),
                                                 SD.He = round(apply(results.he, 1, sd), 2),
                                                 Signif.He = results.significance.he
                                                 
    )
    
    main.results.data.frame <-  data.frame( Population = 1:length(unique(clustering.vector)) ,
                                            . = rep("        ", length(unique(clustering.vector)) ) ,
                                            A = round(apply(results.richness, 1, mean), 2) ,
                                            SD.A = round(apply(results.richness, 1, sd), 2) ,
                                            Signif.A = results.significance.richness ,
                                            .. = rep("        ", length(unique(clustering.vector)) ) ,
                                            PA = round(apply(results.unique, 1, mean), 2) ,
                                            SD.PA = round(apply(results.unique, 1, sd), 2) ,
                                            Signif.PA = results.significance.unique ,
                                            ... = rep("        ", length(unique(clustering.vector)) ) ,
                                            He = round(apply(results.he, 1, mean), 2),
                                            SD.He = round(apply(results.he, 1, sd), 2),
                                            Signif.He = results.significance.he
                                            
    )
    
  }
  
  if (savefile) {
    write.table(  main.results.data.frame.file , quote = FALSE , file = paste(save.filename, ".txt", sep = ""), sep = ";", row.names = FALSE, col.names = TRUE, na = "NA", dec = ".")
    
  }
  
  cat("\n")
  cat("\n")
  cat(" Rambo Results for", number.iteractions," randomizations using", population.n, "populations with", resamp, "individuals" )
  cat("\n")
  cat(" ------------------------------------------------------------------------------------------------")
  cat("\n")
  cat("\n")
  print(main.results.data.frame)
  cat("\n")
  cat(" ------------------------------------------------------------------------------------------------")
  cat("\n")
  cat(" NS for non-significant Wilcoxon tests p-values for differences in mean diversity each population and the global pool (i.e. it is a paired difference test).")
  cat("\n")
  cat("\n")
  
}
