## Identifying miRNA sponge interactions using Sensitivity Correlation (SC) method
# miRExp: miRNA expression data, rows are samples and columns are miRNAs.
# ceRExp: ceRNA expression data, rows are samples and columns are ceRNAs (lncRNAs, pseudogenes, etc.).
# mRExp: miRNA expression data, rows are samples and columns are mRNAs.
# miRTarget: Putative miRNA-target interactions.
# minSharedmiR: Minimum number of sharing miRNAs for rach ceRNA-mRNA pair.
# pvaluecutoff: Cutoff value for p-values.
# poscorcutoff: Cutoff value of positive correlation.
# senscorcutoff: Cutoff value of sensitivity correlation.
SC <- function(miRExp, ceRExp, mRExp, miRTarget, minSharedmiR = 1, 
                 pvaluecutoff = 0.05, poscorcutoff = 0, senscorcutoff = 0.3,
		 num.cores = 6){   
         
        miRTarget <- as.matrix(miRTarget)	
        miRceR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(ceRExp))), ]
        miRmR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(mRExp))), ]

	ceRSym <- unique(miRceR[, 2])
        mRSym <- unique(miRmR[, 2])
	miRSym <- unique(c(miRceR[, 1], miRmR[, 1]))

	m <- length(ceRSym)
	n <- length(mRSym)

        ceRExp_query <- ceRExp[, which(colnames(ceRExp) %in% ceRSym)]
        mRExp_query <- mRExp[, which(colnames(mRExp) %in% mRSym)]
	
        Cor.Pvalue <- corAndPvalue(ceRExp_query, mRExp_query)
        
	index <- which(Cor.Pvalue$cor > poscorcutoff & Cor.Pvalue$p < pvaluecutoff, arr.ind = TRUE)        
	
        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        Res <- foreach(i = seq_len(nrow(index)), .packages = "corpcor") %dopar% {
	    
	        Interin1 <- miRceR[which(miRceR[, 2] %in% ceRSym[index[i, 1]]), 1]
                Interin2 <- miRmR[which(miRmR[, 2] %in% mRSym[index[i, 2]]), 1]

		M1 <- length(Interin1)
                M2 <- length(Interin2)
                M3 <- length(intersect(Interin1, Interin2))
                M4 <- length(miRSym)
                M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

                if (M3 >= minSharedmiR & M5 < pvaluecutoff) {
                    
                    C1 <- ceRSym[index[i, 1]]
                    C2 <- mRSym[index[i, 2]]

                    ceRExpIdx <- which(colnames(ceRExp) %in% ceRSym[index[i, 1]])
                    mRExpIdx <- which(colnames(mRExp) %in% mRSym[index[i, 2]])
                    miRExpIdx <- which(colnames(miRExp) %in% intersect(Interin1, Interin2))

		    M6 <- Cor.Pvalue$cor[index[i, 1], index[i, 2]]
		    M7 <- Cor.Pvalue$p[index[i, 1], index[i, 2]]
                    M8 <- M6 - corpcor::pcor.shrink(cbind(ceRExp[, ceRExpIdx], mRExp[, mRExpIdx],
                                                    miRExp[, miRExpIdx]), verbose = FALSE)[1, 2]
                } else {

                C1 <- NA; C2 <- NA; M6 <- NA; M7 <- NA; M8 <- NA 
	       
                }
	       
	        tmp <- c(C1, C2, M3, M5, M6, M7, M8)    
                return(tmp)
	}
        
	# shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
        
	Res <- do.call(rbind, Res)
        Res <- Res[which(as.numeric(Res[, 7]) > senscorcutoff), ]

        colnames(Res) <- c("sponge_1", "sponge_2", "#Shared miRNAs", 
                           "Sig. p.value of sharing miRNAs", 
			   "Correlation", 
                           "Sig. p.value of correlation", 			   
			   "Sensitivity correlation")
			   
       rownames(Res) <- seq_len(nrow(Res))

return(Res)

}

## SC method for miRNA sponge interaction networks
source(¡®Scripts/Preprocess.R')
library(igraph)
library(WGCNA)
library(corpcor)
library(doParallel)
miRlncR <- read.csv("LncBase_v2.0+NPInter_4.0.csv", header = TRUE, sep = ",")
miRmR <- read.csv("miRTarBase_v8.0+TarBase_v8.0.csv", header = TRUE, sep = ",")
miRpseudo <- read.csv("ENCORI.csv", header = TRUE, sep = ",")
colnames(miRlncR) <- colnames(miRmR) <- colnames(miRpseudo) <- c("miRNA", "target")
miRTarget_lncR_vs_mR <- unique(rbind(miRlncR, miRmR))
miRTarget_pseudo_vs_mR <- unique(rbind(miRpseudo, miRmR))
miRTarget_lncR_vs_pseudo <- unique(rbind(miRlncR, miRpseudo))

SC_ceRNA_network_lncR_vs_mR <- SC(Pancancer_miRNA_Exp_DEG[[2]], Pancancer_lncRNA_Exp_DEG[[2]], 
                                      Pancancer_mRNA_Exp_DEG[[2]], miRTarget_lncR_vs_mR, 
		                      pvaluecutoff = 0.05, senscorcutoff = 0.3)

SC_ceRNA_network_pseudo_vs_mR <- SC(Pancancer_miRNA_Exp_DEG[[2]], Pancancer_pseudogene_Exp_DEG[[2]], 
                                      Pancancer_mRNA_Exp_DEG[[2]], miRTarget_pseudo_vs_mR, 
		                      pvaluecutoff = 0.05, senscorcutoff = 0.3)

## MCL method for miRNA sponge modules
library(miRspongeR)
library(GSEABase)
SC_ceRNA_module_lncR_vs_mR_MCL <- netModule(SC_ceRNA_network_lncR_vs_mR[, 1:2], method = "MCL", modulesize = 3)
SC_ceRNA_module_pseudo_vs_mR_MCL <- netModule(SC_ceRNA_network_pseudo_vs_mR[, 1:2], method = "MCL", modulesize = 3)

SC_ceRNA_finalmodule_lncR_vs_mR_MCL <- geneIds(CandModgenes(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], SC_ceRNA_module_lncR_vs_mR_MCL))
SC_ceRNA_finalmodule_pseudo_vs_mR_MCL <- geneIds(CandModgenes(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], SC_ceRNA_module_pseudo_vs_mR_MCL))
