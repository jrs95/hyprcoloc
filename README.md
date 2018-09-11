# hyprmtc
Hypothesis prioritisation multi-trait colocalisation.

# Functions
* hyprmtc - performs multi-trait colocalisation across numerous traits.  

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/hyprmtc")
4. library(hyprmtc)

# Example
betas <- hyprmtc::betas   
ses <- hyprmtc::ses  
corr <- hyprmtc::corr   
hyprmtc(betas, ses, corr)  
