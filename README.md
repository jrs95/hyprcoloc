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
betas <- hyprmtc::test.betas  
head(betas)  
ses <- hyprmtc::test.ses  
head(ses)  
traits <- paste0("T", 1:10)  
rsid <- rownames(betas)    
hyprmtc(betas, ses, trait.names=traits, snp.id=rsid)  
