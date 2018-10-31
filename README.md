# hyprmtc
Hypothesis Prioritisation Multi-Trait Colocalization.

# Functions
* hyprmtc - performs multi-trait colocalization across numerous traits.  

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/hyprmtc")
4. library(hyprmtc)

# Example
betas <- hyprmtc::betas  
head(betas)  
ses <- hyprmtc::ses  
head(ses)  
traits <- paste0("T", 1:10)  
rsid <- rownames(betas)    
ld <- hyprmtc::ld  
corr <- diag(10)  
hyprmtc(betas, ses, trait.names=traits, snp.id=rsid, ld.matrix=ld, trait.cor=corr, n.cvs=1, bb.alg=TRUE)  
