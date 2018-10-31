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
print(traits)   
rsid <- rownames(betas)
print(rsid)  
ld <- hyprmtc::ld  
print(ld[1:5,1:5])   
corr <- diag(10)
print(corr)  
hyprmtc(betas, ses, trait.names=traits, snp.ind=rsid, ld.matrix=ld, trait.corr=corr, n.cvs=1, bb.alg=TRUE)  
