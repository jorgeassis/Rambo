# Standardized population genetic diversity

Genetic diversity estimation per site or group of sites as standardised allelic richness, standardised number of private alleles and expected heterozygosity determined to the smallest site or group in terms of individuals under a set of randomizations.

## Getting Started

Instructions to get the project up and running on your local machine.

### Prerequisites

Install the last verion of R available at [The Comprehensive R Archive Network](https://cran.r-project.org/).

### Running the code

Open R and set the working directory (path to) <br>
Load the main function "Rambo" into memory <br>

```
# Source main function
source("https://raw.githubusercontent.com/jorgeassis/rambo/master/Script.R")
main.data.file <- "../Data/Genetics/laminariaPallidaGenPop.gtx"
missing.data <- 0
replace <- FALSE
ncode <- 3
resample.number.auto <- TRUE
resample.number <- 20
discard.pops <- NULL
number.iteractions <- 999
alfa.test <- 0.05
savefile <- FALSE
save.filename <- "richnessK2"
# Define a vector to cluster populations
clustering.vector <- 1:26 # c(1,1,2,2)
clustering.vector <- c(rep(1,11),2,2,rep(3,13))
clustering.vector <- c(rep(1,11),2,2,rep(3,6),rep(4,7))

# Run function
Rambo(main.data.file, missing.data, ncode, replace, resample.number.auto, resample.number, discard.pops, number.iteractions, alfa.test, clustering.vector, savefile)
```

## Authors

* **Jorge Assis** @ [theMarineDataScientist](https://medium.com/@jorgemfa)

## License

Except where otherwise noted, the content on this repository is licensed under a [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).

## Citation

Assis, J., Coelho, N. C., Lamy, T., Valero, M., Alberto, F., & Serrão, E. A. (2016). Deep reefs are climatic refugia for genetic diversity of marine forests. Journal of Biogeography, (43), 833–844. https://doi.org/10.1111/jbi.12677

