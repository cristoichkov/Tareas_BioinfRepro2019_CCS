#### 4.1.3 Running process_radtags

Here is how single-end data received from an Illumina sequencer might look:

```
% ls ./raw
lane3_NoIndex_L003_R1_001.fastq.gz  lane3_NoIndex_L003_R1_006.fastq.gz  lane3_NoIndex_L003_R1_011.fastq.gz
lane3_NoIndex_L003_R1_002.fastq.gz  lane3_NoIndex_L003_R1_007.fastq.gz  lane3_NoIndex_L003_R1_012.fastq.gz
lane3_NoIndex_L003_R1_003.fastq.gz  lane3_NoIndex_L003_R1_008.fastq.gz  lane3_NoIndex_L003_R1_013.fastq.gz
lane3_NoIndex_L003_R1_004.fastq.gz  lane3_NoIndex_L003_R1_009.fastq.gz
lane3_NoIndex_L003_R1_005.fastq.gz  lane3_NoIndex_L003_R1_010.fastq.gz
```
Then you can run <span style="color:green;">**process_radtags**</span> in the following way:

```
% process_radtags -p ./raw/ -o ./samples/ -b ./barcodes/barcodes_lane3 \
                  -e sbfI -r -c -q
```
I specify the directory containing the input files, ```./raw```, the directory I want <span style="color:green">**process_radtags**</span> to enter the output files, ./samples, and a file containing my barcodes, ```./barcodes/barcodes_lane3```, along with the restrction enzyme I used and instructions to clean the data and correct barcodes and restriction enzyme cutsites (```-r```, ```-c```, ```-q```).

Here is a more complex example, using paired-end double-digested data (two restriction enzymes) with combinatorial barcodes, and gzipped input files. Here is what the raw Illumina files may look like:
```
% ls ./raw
GfddRAD1_005_ATCACG_L007_R1_001.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_001.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_002.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_002.fastq.gz
GfddRAD1_005_ATCACG_L```007_R1_003.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_003.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_004.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_004.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_005.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_005.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_006.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_006.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_007.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_007.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_008.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_008.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_009.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_009.fastq.gz
```
Now we specify both restriction enzymes using the ```--renz_1``` and ```--renz_2``` flags along with the type combinatorial barcoding used. Here is the command:
```
% process_radtags -P -p ./raw  -b ./barcodes/barcodes -o ./samples/ \
                  -c -q -r --inline_index --renz_1 nlaIII --renz_2 mluCI
```

#### The output of process_radtags
The output of the <span style="color:green">**process_radtags**</span> differs depending if you are processing single-end or paired-end data. In the case of single-end reads, the program will output one file per barcode into the output directory you specify. If the data do not have barcodes, then the file will retain its original name.

If you are processing paired-end reads, then you will get four files per barcode, two for the single-end read and two for the paired-end read. For example, given barcode ACTCG, you would see the following four files:
```
sample_ACTCG.1.fq
sample_ACTCG.rem.1.fq
sample_ACTCG.2.fq
sample_ACTCG.rem.2.fq
```

The process_radtags program wants to keep the reads in phase, so that the first read in the ```sample_XXX.1.fq``` file is the mate of the first read in the ```sample_XXX.2.fq``` file. Likewise for the second pair of reads being the second record in each of the two files and so on. When one read in a pair is discarded due to low quality or a missing restriction enzyme cut site, the remaining read can't simply be output to the ```sample_XXX.1.fq``` or ```sample_XXX.2.fq``` files as it would cause the remaining reads to fall out of phase. Instead, this read is considered a remainder read and is output into the ```sample_XXX.rem.1.fq``` file if the paired-end was discarded, or the ```sample_XXX.rem.2.fq``` file if the single-end was discarded.

#### Modifying how process_radtags executes
The <span style="color:green">**process_radtags**</span> program can be modified in several ways. If your data do not have barcodes, omit the barcodes file and the program will not try to demultiplex the data. You can also disable the checking of the restriction enzyme cut site, or modify what types of quality are checked for. So, the program can be modified to only demultiplex and not clean, clean but not demultiplex, or some combination.

There is additional information available in <span style="color:green">**process_radtags**</span> <a href="http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php" title="Hobbit lifestyles">manual page.</a>
