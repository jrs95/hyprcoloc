# hyprmtc
Hypothesis Prioritisation Multi-Trait Colocalization.

# Functions
* hyprmtc - performs multi-trait colocalization across numerous traits.  

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/hyprmtc", build_vignettes=T)
4. library(hyprmtc)

# Example
\# Regression coefficients and standard errors from ten GWAS studies (Traits 1-5, 6-8 & 9-10 colocalize)
betas <- hyprmtc::test.betas
head(betas)
ses <- hyprmtc::test.ses
head(ses)
  
\# Trait names and SNP IDs
traits <- paste0("T", 1:10)
rsid <- rownames(betas)

\# Colocalization analysis
hyprmtc(betas, ses, trait.names=traits, snp.id=rsid)  
