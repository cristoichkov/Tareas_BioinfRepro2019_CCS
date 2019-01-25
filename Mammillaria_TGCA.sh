#!/bin/bash
## Este script baja 5 secuencias de Mammillaria de NCBI
## Baja 5 secuencias en un archivo llamado Mamm.fasta
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&id=AY545239.1,AY545256.1,AY545268.1,AY545272.1,AY545323.1" > Mamm.fasta
## Busca el genero Mammillaria y la especie, así como la secuencia 'TGCA' que contiene cada una
grep -oE ' \Mammillaria \w+ |TGCA' Mamm.fasta
