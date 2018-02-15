# Rambo

Genetic diversity estimation per site or group of sites as standardised allelic richness, standardised number of private alleles and expected heterozygosity determined to the smallest site or group in terms of individuals under a set of randomizations.


Citation:
Assis, J., Coelho, N. C., Lamy, T., Valero, M., Alberto, F., & Serrão, E. A. (2016). Deep reefs are climatic refugia for genetic diversity of marine forests. Journal of Biogeography, (43), 833–844. https://doi.org/10.1111/jbi.12677


# Usage

main.data.file <- "example.file.gen"
missing.data <- 0
ncode <- 3
replace <- FALSE
resample.number.auto <- FALSE
resample.number <- 20
discard.pops <- NULL
number.iteractions <- 999
alfa.test <- 0.05
clustering.vector <- NULL
savefile <- TRUE
save.filename <- "richness"

Rambo(main.data.file, 
      missing.data, 
      ncode, 
      replace, 
      resample.number.auto, 
      resample.number,
      discard.pops, 
      number.iteractions, 
      alfa.test, 
      clustering.vector, 
      savefile,save.filename)
      
