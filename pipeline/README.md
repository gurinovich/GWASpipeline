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

conda install -c bioconda vcftools
conda install -c bioconda bcftools
conda install -c bioconda plink
conda install -c bioconda annovar

```

### Clone Repository
```bash

$ git clone https://github.com/gurinovich/GWASpipeline

```
### Initalize Paths to Test Data
```bash

$ cd GWASpipeline/pipeline 
$ python utils/paths.py  

```

### Download Nextflow Executable
```

$ curl -s https://get.nextflow.io | bash

```

### Locally Run Example Data
```bash
$ ./module load R/3.6.0
$ ./module load vcftools
$ ./module load bcftools
$ ./module load plink/2.00a1LM
$ ./module load annovar/2018apr
$ ./mkdir results

$ ./nextflow gwas.nf -c gwas.config

```

### Expected Output
```bash

N E X T F L O W  ~  version 19.04.1
Launching `gwas.nf` [irreverent_spence] - revision: 8a796af47a
-

G W A S  ~  P I P E L I N E

================================
indir     : /restricted/projectnb/necs/Zeyuan_Analysis/GWASpipeline/GWASpipeline/pipeline/data/
outdir    : /restricted/projectnb/necs/Zeyuan_Analysis/GWASpipeline/GWASpipeline/pipeline/results

vcf       : /restricted/projectnb/necs/Zeyuan_Analysis/GWASpipeline/GWASpipeline/pipeline/data//toy_vcf.csv
pheno     : /restricted/projectnb/necs/Zeyuan_Analysis/GWASpipeline/GWASpipeline/pipeline/data//pheno_file_logistic.csv
snpset    : /restricted/projectnb/necs/Zeyuan_Analysis/GWASpipeline/GWASpipeline/pipeline/data//snpset.txt

phenotypes: outcome
covars    : 'age,sex,PC1,PC2,PC3,PC4'
model     : logistic
test      : Score

-
[warm up] executor > local
executor >  local (92)
[48/a51761] process > qc_miss    [100%] 22 of 22 ✔
[bd/b4c149] process > qc_mono    [100%] 22 of 22 ✔
[87/979c47] process > vcf_to_gds [100%] 22 of 22 ✔
[29/592b55] process > merge_gds  [100%] 1 of 1 ✔
[0e/e0ef74] process > pcair      [100%] 1 of 1 ✔
[ce/334838] process > pcrelate   [100%] 1 of 1 ✔
[23/442a05] process > nullmod    [100%] 1 of 1 ✔
[b3/da7597] process > gwas       [100%] 22 of 22 ✔
Completed at: 11-May-2020 15:29:54
Duration    : 1m 6s
CPU hours   : 0.1
Succeeded   : 92

```
