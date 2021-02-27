## Identifying miRNA sponge interactions using Sparse Partial correlation ON Gene Expression (SPONGE) method
# miRExp: miRNA expression data, rows are samples and columns are miRNAs.
# ceRExp: ceRNA expression data, rows are samples and columns are ceRNAs (lncRNAs, pseudogenes, etc.).
# mRExp: miRNA expression data, rows are samples and columns are mRNAs.
# miRTarget: Putative miRNA-target interactions.
# coefficient.threshold: Threshold to cross for a regression coefficient to be called significant.
# min.cor: Consider only ceRNA-mRNA pairs with a minimum correlation specified here.
# minSharedmiR: Minimum number of sharing miRNAs for rach ceRNA-mRNA pair.
# pvaluecutoff: Cutoff value for p-values.
# number_of_datasets: The number of datesets defining the precision of the p-value when building null model for p-value computation.
SPONGE <- function(miRExp, ceRExp, mRExp, miRTarget, coefficient.threshold = -0.05, 
                   min.cor = 0.1, minSharedmiR = 1, pvaluecutoff = 0.05, 
		   null_model, num.cores = 6){   
        
	miRTarget <- as.matrix(miRTarget)	
        miRTarget <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(cbind(ceRExp, mRExp)))), ]
        
	miRlist <- unique(miRTarget[, 1])
	Tarlist <- unique(miRTarget[, 2])
        miRTarget_graph <- make_graph(c(t(miRTarget)), directed = FALSE)
	miRTarget_adjacency <- as_adjacency_matrix(miRTarget_graph, sparse = FALSE)
        miRTarget_adjacency_matrix <- miRTarget_adjacency[which(rownames(miRTarget_adjacency) %in% Tarlist), which(colnames(miRTarget_adjacency) %in% miRlist)]

	# get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)
        
	miRTarget_candidates <- sponge_gene_miRNA_interaction_filter(gene_expr = cbind(ceRExp, mRExp),
                                                                     mir_expr = miRExp,
                                                                     mir_predicted_targets = miRTarget_adjacency_matrix,
								     coefficient.threshold = coefficient.threshold)
	
	ceRlist <- names(miRTarget_candidates)[which(names(miRTarget_candidates) %in% colnames(ceRExp))]
	mRlist <- names(miRTarget_candidates)[which(names(miRTarget_candidates) %in% colnames(mRExp))]
	gene.combinations <- expand.grid(ceRlist, mRlist)

        ceRNA_interactions <- sponge(gene_expr = cbind(ceRExp, mRExp),
                                     mir_expr = miRExp,
                                     mir_interactions = miRTarget_candidates,
				     gene.combinations = gene.combinations,
				     min.cor = min.cor)

	# shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

	ceRNA_interactions_filter <- ceRNA_interactions[which(ceRNA_interactions$df >= minSharedmiR), ]
	
	ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions_filter, 
                                                   null_model = null_model)

	Res <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < pvaluecutoff),]

return(Res)

}

## SPONGE method for miRNA sponge interaction networks
source(¡®Scripts/Preprocess.R')
library(igraph)
library(SPONGE)
library(doParallel)
miRlncR <- read.csv("LncBase_v2.0+NPInter_4.0.csv", header = TRUE, sep = ",")
miRmR <- read.csv("miRTarBase_v8.0+TarBase_v8.0.csv", header = TRUE, sep = ",")
miRpseudo <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
colnames(miRlncR) <- colnames(miRmR) <- colnames(miRpseudo) <- c("miRNA", "target")
miRTarget_lncR_vs_mR <- unique(rbind(miRlncR, miRmR))
miRTarget_pseudo_vs_mR <- unique(rbind(miRpseudo, miRmR))
miRTarget_lncR_vs_pseudo <- unique(rbind(miRlncR, miRpseudo))
pre_null_model <- sponge_build_null_model(number_of_datasets = 1e4, number_of_samples = nrow(Pancancer_miRNA_Exp_DEG[[2]]))

SPONGE_ceRNA_network_lncR_vs_mR <- SPONGE(Pancancer_miRNA_Exp_DEG[[2]], Pancancer_lncRNA_Exp_DEG[[2]], 
                                          Pancancer_mRNA_Exp_DEG[[2]], miRTarget_lncR_vs_mR, 
				          pvaluecutoff = 0.05, null_model = pre_null_model)

SPONGE_ceRNA_network_pseudo_vs_mR <- SPONGE(Pancancer_miRNA_Exp_DEG[[2]], Pancancer_pseudogene_Exp_DEG[[2]], 
                                            Pancancer_mRNA_Exp_DEG[[2]], miRTarget_pseudo_vs_mR, 
				            pvaluecutoff = 0.05, null_model = pre_null_model)

## MCL method for miRNA sponge modules
library(miRspongeR)
library(GSEABase)
SPONGE_ceRNA_module_lncR_vs_mR_MCL <- netModule(SPONGE_ceRNA_network_lncR_vs_mR[, 1:2], method = "MCL", modulesize = 3)
SPONGE_ceRNA_module_pseudo_vs_mR_MCL <- netModule(SPONGE_ceRNA_network_pseudo_vs_mR[, 1:2], method = "MCL", modulesize = 3)

SPONGE_ceRNA_finalmodule_lncR_vs_mR_MCL <- geneIds(CandModgenes(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], SPONGE_ceRNA_module_lncR_vs_mR_MCL))
SPONGE_ceRNA_finalmodule_pseudo_vs_mR_MCL <- geneIds(CandModgenes(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], SPONGE_ceRNA_module_pseudo_vs_mR_MCL))

