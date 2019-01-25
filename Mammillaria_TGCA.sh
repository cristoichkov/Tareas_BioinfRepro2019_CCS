#!/bin/bash

curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&id=AY545239.1,AY545256.1,AY545268.1,AY545272.1,AY545323.1" > Mamm.fasta
grep -oE ' \Mammillaria \w+ |TCGA' Mamm.fasta
