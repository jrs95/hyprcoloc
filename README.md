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
ses <- hyprmtc::ses  
corr <- hyprmtc::corr   
traits <- hyprmtc::traits  
rsid <- hyprmtc::rsid
hyprmtc(betas, ses, traits, rsid, corr)  
