# GWASpipeline
The goal of this project is to create a shareable GWAS pipeline. The initial focus of this effort is to define the feature list, output format and the tools/languages used in the development of this pipeline. The proposed steps in this pipeline are shown below.

1. VCF -> GDS
2. GDS - > PCA
3. PCA <-> GRM
4. GRM -> Null model + GWAS <- Phenotype
5. GWAS -> Output
6. Output -> Bayesian analysis -> updated output
