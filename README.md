# SpongeModule
## Introduction
Predicting competing endogenous RNA (ceRNA) or microRNA (miRNA) sponge modules is a challenge and meaningful task in revealing ceRNA regulation mechanism at the module level. Modules upon this backdrop represent groups of miRNA sponges with mutual competition, which are prone to be functional units for achieving biological processes. The existing in silico methods of identifying miRNA sponge modules are categorized into three types: (i) network-based clustering, (ii) matrix factorization, and (iii) step-wise evaluation. 

## Five representative module discovery methods for comparison

-Network-based clustering: SC+MCL, SPONGE+MCL

-Matrix factorization: jointNMF, LAceModule

-Step-wise evaluation: LMSM

## The usage of five representative module discovery methods
Paste all files including scripts and datasets into a single folder (set the folder as the directory of Matlab and R environment), the scripts of five representative module discovery methods are implemented in five folders (jointNMF, LAceModule, SC+MCL, SPONGE+MCL, LMSM). It is noted that some scripts are running on Matlab, and some scripts are running on R. For example, users can simply run the R scripts to identify miRNA sponge modules by using LMSM method as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("LMSM.R")
```
