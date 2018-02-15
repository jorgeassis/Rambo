# Rambo

Genetic diversity estimation per site or group of sites as standardised allelic richness, standardised number of private alleles and expected heterozygosity determined to the smallest site or group in terms of individuals under a set of randomizations.


Citation:
Assis, J., Coelho, N. C., Lamy, T., Valero, M., Alberto, F., & Serrão, E. A. (2016). Deep reefs are climatic refugia for genetic diversity of marine forests. Journal of Biogeography, (43), 833–844. https://doi.org/10.1111/jbi.12677


# Usage

main.data.file <- "example.file.gen" <br />
missing.data <- 0 <br />
ncode <- 3 <br />
replace <- FALSE <br />
resample.number.auto <- FALSE <br />
resample.number <- 20 <br />
discard.pops <- NULL <br />
number.iteractions <- 999 <br />
alfa.test <- 0.05 <br />
clustering.vector <- NULL <br />
savefile <- TRUE <br />
save.filename <- "richness" <br /> <br />

Rambo(main.data.file,  <br />
      missing.data,  <br />
      ncode,  <br />
      replace,  <br />
      resample.number.auto,  <br />
      resample.number, <br />
      discard.pops,  <br />
      number.iteractions,  <br />
      alfa.test,  <br />
      clustering.vector,  <br />
      savefile,save.filename) <br /> <br />
      
