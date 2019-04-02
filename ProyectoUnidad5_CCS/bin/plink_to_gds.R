#### Convert plink data format to gds data format ####
library(SNPRelate)


## load the data in plink format and convert them to gds
snpgdsBED2GDS("../data/maicesArtegaetal2015.bed", 
              "../data/maicesArtegaetal2015.fam", 
              "../data/maicesArtegaetal2015.bim", 
              out.gdsfn="../data/maicesArtegaetal2015.gds", 
              option = snpgdsOption(Z=10)) # 10 cromosomas


## show the summary of the gds file
snpgdsSummary("../data/maicesArtegaetal2015.gds")
