# HyPrColoc <img src="man/figures/logo.png" align="right" height="139"/>
Hypothesis Prioritisation in multi-trait Colocalization (HyPrColoc).  

Genome-wide association studies (GWAS) have identified thousands of genomic regions affecting complex diseases. The next challenge is to elucidate the causal genes and mechanisms involved. One approach is to use statistical colocalization to assess shared genetic aetiology across multiple related traits (e.g. molecular traits, metabolic pathways and complex diseases) to identify causal pathways, prioritize causal variants and evaluate pleiotropy.  

HyPrColoc is an efficient deterministic Bayesian divisive clustering algorithm using GWAS summary statistics that can detect colocalization across vast numbers of traits simultaneously.  

## Functions
* `hyprcoloc`: identifies clusters of colocalized traits and candidate causal SNPs using the HyPrColoc Bayesian divisive clustering algorithm.  
* `sensitivity_plot`: plots a heatmap to visualise how stable the clusters of colocalized traits are to variations in the algorithms input parameters.  
* `cred.sets`: computes credible sets of SNPs for each cluster of colocalized traits.  

## Installation
```
install.packages("remotes")
remotes::install_version("RcppEigen", version = "0.3.3.9.3")
remotes::install_github("jrs95/hyprcoloc", build_vignettes = TRUE)
```

## Example
```
# Library
library(hyprcoloc)

# Regression coefficients from ten studies (Traits 1-5, 6-8 & 9-10 colocalize)
betas <- hyprcoloc::test.betas
head(betas)
ses <- hyprcoloc::test.ses
head(ses)

# Trait names and SNP IDs
traits <- paste0("T", 1:10)
rsid <- rownames(betas)

# Colocalization analysis
hyprcoloc(betas, ses, trait.names = traits, snp.id = rsid)
```

## Citation
Foley CN, Staley JR, *et al.* A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. [Nat Commun](https://pubmed.ncbi.nlm.nih.gov/33536417/) 2021;12(1):764.  
