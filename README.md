# Rambo

Genetic diversity estimation per site or group of sites as standardised allelic richness, standardised number of private alleles and expected heterozygosity determined to the smallest site or group in terms of individuals under a set of randomizations.

## How it works?

A high-resolution polygon is converted to an infinite resistance surface. <br>
Minimum distances between sites are computed with a shortest path algorithm considering the infinite resistance of land and null resistance throughout the sea. <br>
The outcomes are a matrix of pairwise distances, a figure to visualize if sites are well represented in the study area and a figure depicting an example of a shortest marine distance. <br>
The main file with the sites should be structured as “Name Lon Lat” or “Name Lat Lon”. Coordinates must be in decimal degrees.

## Getting Started

Instructions to get the project up and running on your local machine.

### Prerequisites

Install the last verion of R available at [The Comprehensive R Archive Network](https://cran.r-project.org/).

### Running the code

Open R and set the working directory (path to) <br>
Load the main function "Rambo" into memory <br>

```
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

Rambo(main.data.file, missing.data, ncode, replace, resample.number.auto, resample.number, discard.pops, number.iteractions, alfa.test, clustering.vector, savefile, save.filename)
```

## Authors

* **Jorge Assis** @ [Bears Studio](https://www.bears.studio)

## License

Except where otherwise noted, the content on this repository is licensed under a [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).

## Citation

Assis, J., Coelho, N. C., Lamy, T., Valero, M., Alberto, F., & Serrão, E. A. (2016). Deep reefs are climatic refugia for genetic diversity of marine forests. Journal of Biogeography, (43), 833–844. https://doi.org/10.1111/jbi.12677

