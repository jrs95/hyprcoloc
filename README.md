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
traits <- hyprmtc::traits  
print(traits)   
rsid <- hyprmtc::rsid
print(rsid)  
ld <- hyprmtc::ld
print(ld[1:5,1:5])   
corr <- hyprmtc::corr
print(corr)  
hyprmtc(betas, ses, trait.names=traits, snp.ind=rsid, ld.matrix=ld, trait.corr=corr, n.cvs=1, bb.alg=TRUE)  
