# LMSM method for miRNA sponge modules
source(¡®Scripts/Preprocess.R')
library(miRSM)
library(SummarizedExperiment)
library(GSEABase)
library(WGCNA)

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

WGCNA_lncR_vs_mR_Coexpression <- module_WGCNA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, Pancancer_mRNA_Exp_DEG_SummarizedExperiment, RsquaredCut = 0.8)
miRSM_ceRNA_module_lncR_vs_mR_WGCNA <- miRSM(Pancancer_miRNA_Exp_DEG_SummarizedExperiment, 
                                             Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
				             Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
                                             miRTarget_lncR_vs_mR_SummarizedExperiment, 
				             WGCNA_lncR_vs_mR_Coexpression, 
				             method = "SCC", MC.cutoff = 0.6)

WGCNA_pseudo_vs_mR_Coexpression <- module_WGCNA(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, Pancancer_mRNA_Exp_DEG_SummarizedExperiment, RsquaredCut = 0.8)
miRSM_ceRNA_module_pseudo_vs_mR_WGCNA <- miRSM(Pancancer_miRNA_Exp_DEG_SummarizedExperiment, 
                                               Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, 
				               Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
                                               miRTarget_pseudo_vs_mR_SummarizedExperiment, 
				               WGCNA_pseudo_vs_mR_Coexpression, 
				               method = "SCC", MC.cutoff = 0.6)


