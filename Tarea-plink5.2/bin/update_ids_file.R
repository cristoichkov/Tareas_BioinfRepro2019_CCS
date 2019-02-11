### script to create update file (.txt) of FID (Family ID) and IID (Individual ID) for plink

rm(list = ls())
library(dplyr)
library(tidyr)

##### Get data #####
## Info of maices
maices_meta <- read.delim("./meta/maizteocintle_SNP50k_meta_extended.txt")

## Info of old Ids
maices_old_ids <- read.delim("./meta/maices_old_ids.txt", header = FALSE)

## Separate in two columns old FID and old IID
maices_old_ids <- separate(maices_old_ids, col = "V1", 
                             into = c("OFID", "OIID"), sep = " ") 

## Get the column Categ.Altitud for new FID
NFID <- maices_meta %>%
  select(Categ.Altitud) %>%
  rename(NFID = Categ.Altitud)

## Combine old IDs and new ID
maices_old_ids <-  cbind(maices_old_ids, NFID) 

## Clone column old IID to obteind new IID
maices_old_ids <- mutate(maices_old_ids, NIID = OIID)

## Write the result in txt file, separate for space, without name of columns and rows.
write.table(maices_old_ids, "./meta/maices_update_ids_plink.txt", sep = " ", col.names=FALSE, row.names=FALSE, quote=FALSE)


