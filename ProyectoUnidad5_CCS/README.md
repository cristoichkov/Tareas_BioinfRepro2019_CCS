# README

## Interactive visualization using Shiny package and Arteaga et al. (2016) maize dataset

This repository contains the development of a desktop application using [Shiny Dashboard](https://rstudio.github.io/shinydashboard/). The application generates an interactive plot, map, and a table with an example of an exploratory genomic analysis of populations with the [Arteaga et al. (2016)](https://www.sciencedirect.com/science/article/pii/S2213596015300714?via%3Dihub) maize dataset. Data deposited in the Dryad repository: [http://dx.doi.org/10.5061/dryad.4t20n](http://dx.doi.org/10.5061/dryad.4t20n)

#### Prerequisites
##### Software:
- R 3.5.3

##### R packages:
- shiny 1.2.0
- ggplot2 3.1.0
- dplyr 0.8.0.1
- raster 2.8.19
- shinydashboard 0.7.1
- plotly 4.8.0
- leaflet 2.0.2
- SNPRelate 1.16.0

#### Directories:
###### bin
Contains scripts `bin/plink_to_gds.R` and `bin/app.R` used for the analysis.
- First you have to execute the `bin/plink_to_gds.R` script to generate a gds file that is used in the following script
- The `bin/app.R` script generates the shiny dashboard app. The first part of the script is used to generate and filter the databases used to create three interactive objects:
   * PCA plot.- the user can select the eigenvalue for the x-axis and the and y-axis, and the variable with which the points are colored.
   * Maize distribution map.- this map is associated with the PCA and the points of distribution will be colored according to the variable selected in the graph.   
   * Tables.- the user can download the database of eigenvalues, coordinates and 19 WorldClim variables en formato .csv y .tsv

###### data
Contains the genomic data of maize in format plink:
 - `data/maicesArtegaetal2015.bed`
 - `data/maicesArtegaetal2015.bim`
 - `data/maicesArtegaetal2015.fam`

###### meta
- The file `meta/maizteocintle_SNP50k_meta_extended.txt` contains information about the samples maize, we mainly use latitude and longitude coordinates.

#### Notes
- We run once the `bin/plink_to_gds.R` script to generate `data/maicesArtegaetal2015.gds`
- WorldClim variables are downloaded at a resolution of 2.5 and stored in `data/wc2-5`
- Therefore, in this repository these files are not available, they will be available when the user executes the script


#### Credits
##### Cristian Cervantes
