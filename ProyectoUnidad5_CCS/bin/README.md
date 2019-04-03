# README
### /bin

#### Admixture

- Installation on Linux

Download software for Linux x86 [64]: [admixture_linux-1.3.0.tar.gz](http://software.genetics.ucla.edu/admixture/download.html)
```
$tar -xzvf admixture_linux-1.3.0.tar.gz
$cd admixture_linux-1.3.0
$sudo mv admixture /usr/local/bin/
```

- Run `bin/1-runadmixture.sh` script admixture for K 1-5 and the results are generated in the `data/admixture` folder.

```
$chmod chmod u+x 1-runadmixture.sh
$bash 1-runadmixture.sh
```

#### Plink to gds

- `bin/plink_to_gds.R` script convert plink files to gds files


#### Shiny app

The `bin/app.R` script generates the shiny dashboard app. Before run Shiny app you have to create the admixture and gds files.
