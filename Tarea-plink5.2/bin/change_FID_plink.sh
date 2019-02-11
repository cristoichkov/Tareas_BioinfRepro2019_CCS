#!/bin/bash

## change ids
plink --file data/maicesArtegaetal2015 --update-ids meta/maices_update_ids_plink.txt --recode --out data/maicesArtegaetal2015
