### Clone Repository
```bash

$ git clone https://github.com/gurinovich/GWASpipeline

```

### Conda Envrionment
```bash

conda create -n gwas python=2.7
source activate gwas
conda list

# Installs popular R packages
conda install r-essentials r-base

# R packages not included in essentials
conda config --add channels bioconda
conda install -c bioconda/label/gcc7 bioconductor-seqarray

```

### Download Nextflow Executable
```

$ cd GWASpipeline/pipeline
$ curl -s https://get.nextflow.io | bash

```

### Locally Run Example Data
```bash

$ ./nextflow gwas.nf -c gwas.config

```

### Expected Output
```bash

N E X T F L O W  ~  version 19.04.1
Launching `gwas.nf` [festering_hamilton] - revision: 782fb9722c
-

G W A S  ~  P I P E L I N E

================================
indir     : /Users/anthonyfederico/GitHub/GWASpipeline/pipeline/data
outdir    : /Users/anthonyfederico/GitHub/GWASpipeline/pipeline/results
vcf       : /Users/anthonyfederico/GitHub/GWASpipeline/pipeline/data/toy_vcf.csv

-
[warm up] executor > local
executor >  local (23)
[b1/9519b9] process > vcf_to_gds [100%] 22 of 22 ✔
[8b/663edf] process > merge_gds  [100%] 1 of 1 ✔
Completed at: 02-Jun-2019 18:59:47
Duration    : 31.1s
CPU hours   : (a few seconds)
Succeeded   : 23

```