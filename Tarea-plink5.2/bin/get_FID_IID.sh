#!/bin/bash

##### get FID (Family ID) and IID (Individual ID) #####

cut -d " " -f 1-2 data/maicesArtegaetal2015.ped > meta/maices_old_ids.txt
