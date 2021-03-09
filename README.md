# SpongeModule
## Introduction
Predicting competing endogenous RNA (ceRNA) or microRNA (miRNA) sponge modules is a challenge and meaningful task in revealing ceRNA regulation mechanism at the module level. Modules upon this backdrop represent groups of miRNA sponges with mutual competition, which are prone to be functional units for achieving biological processes. The existing in silico methods of identifying miRNA sponge modules are categorized into three types: (i) network-based clustering, (ii) matrix factorization, and (iii) step-wise evaluation. 

## Five representative module discovery methods for comparison

-Network-based clustering: SC+MCL, SPONGE+MCL

-Matrix factorization: CeModule, LAceModule

-Step-wise evaluation: LMSM

## The usage of five representative module discovery methods
Paste all files including scripts and datasets into a single folder (set the folder as the directory of Matlab and R environment), the scripts of five representative module discovery methods are implemented in the 'Scripts' folder. It is noted that some scripts are running on Matlab, and some scripts are running on R. For example, users can simply run the R scripts to identify miRNA sponge modules using LMSM method as follows.

```{r echo=FALSE, results='hide', message=FALSE}
# Load functions and packages
source('Scripts/Preprocess.R')
library(miRSM)
library(SummarizedExperiment)
library(GSEABase)
library(WGCNA)

# Construct SummarizedExperiment object
Pancancer_miRNA_Exp_DEG_SummarizedExperiment <- Pancancer_miRNA_Exp_DEG[[2]]
rownames(Pancancer_miRNA_Exp_DEG_SummarizedExperiment) <- rownames(Pseudo_sample)
Pancancer_miRNA_Exp_DEG_SummarizedExperiment <- SummarizedExperiment(assays=list(Pancancer_miRNA_Exp_DEG_SummarizedExperiment = Pancancer_miRNA_Exp_DEG_SummarizedExperiment))

Pancancer_lncRNA_Exp_DEG_SummarizedExperiment <- Pancancer_lncRNA_Exp_DEG[[2]]
rownames(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment) <- rownames(Pseudo_sample)
Pancancer_lncRNA_Exp_DEG_SummarizedExperiment <- SummarizedExperiment(assays=list(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment = Pancancer_lncRNA_Exp_DEG_SummarizedExperiment))

Pancancer_pseudogene_Exp_DEG_SummarizedExperiment <- Pancancer_pseudogene_Exp_DEG[[2]]
rownames(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment) <- rownames(Pseudo_sample)
Pancancer_pseudogene_Exp_DEG_SummarizedExperiment <- SummarizedExperiment(assays=list(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment = Pancancer_pseudogene_Exp_DEG_SummarizedExperiment))

Pancancer_mRNA_Exp_DEG_SummarizedExperiment <- Pancancer_mRNA_Exp_DEG[[2]]
rownames(Pancancer_mRNA_Exp_DEG_SummarizedExperiment) <- rownames(Pseudo_sample)
Pancancer_mRNA_Exp_DEG_SummarizedExperiment <- SummarizedExperiment(assays=list(Pancancer_mRNA_Exp_DEG_SummarizedExperiment = Pancancer_mRNA_Exp_DEG_SummarizedExperiment))

miRTarget_lncR_vs_mR_SummarizedExperiment <- as.matrix(miRTarget_lncR_vs_mR)
miRTarget_lncR_vs_mR_SummarizedExperiment <- SummarizedExperiment(assays=list(miRTarget_lncR_vs_mR_SummarizedExperiment = miRTarget_lncR_vs_mR_SummarizedExperiment))

miRTarget_pseudo_vs_mR_SummarizedExperiment <- as.matrix(miRTarget_pseudo_vs_mR)
miRTarget_pseudo_vs_mR_SummarizedExperiment <- SummarizedExperiment(assays=list(miRTarget_pseudo_vs_mR_SummarizedExperiment = miRTarget_pseudo_vs_mR_SummarizedExperiment))

miRTarget_lncR_vs_pseudo_SummarizedExperiment <- as.matrix(miRTarget_lncR_vs_pseudo)
miRTarget_lncR_vs_pseudo_SummarizedExperiment <- SummarizedExperiment(assays=list(miRTarget_lncR_vs_pseudo_SummarizedExperiment = miRTarget_lncR_vs_pseudo_SummarizedExperiment))

# LMSM method for miRNA sponge modules
WGCNA_lncR_vs_mR_Coexpression <- module_WGCNA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, Pancancer_mRNA_Exp_DEG_SummarizedExperiment, RsquaredCut = 0.8)
miRSM_ceRNA_module_lncR_vs_mR_WGCNA <- miRSM(Pancancer_miRNA_Exp_DEG_SummarizedExperiment, 
                                             Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
				             Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
                                             miRTarget_lncR_vs_mR_SummarizedExperiment, 
				             WGCNA_lncR_vs_mR_Coexpression, 
				             method = "SCC")
```
