## Import generic EMT signatures from the reference: Tan TZ, Miow QH, Miki Y, et al. 
## Epithelial-mesenchymal transition spectrum quantification and its efficacy in deciphering survival and drug responses of cancer patients. 
## EMBO Mol Med. 2014;6(10):1279-93. 
## Generic_EMT_signature: a list of generic EMT signatures, first column is generic EMT signatures and second column is EMT types: Epi and Mes.
## Exp_RNA: gene expression data, the rows are genes and the columns are samples.
## Output: samples with EMT types: Epi and Mes.
Pseudo_EMT <- function(Generic_EMT_signature, Exp_RNA){
    
    Epi_ind <- which(Generic_EMT_signature[, 2] %in% "Epi")
    Mes_ind <- which(Generic_EMT_signature[, 2] %in% "Mes")
    Epi <- as.character(Generic_EMT_signature[Epi_ind, 1])
    Mes <- as.character(Generic_EMT_signature[Mes_ind, 1])
    Epi_Update <- Epi[Epi %in% rownames(Exp_RNA)]
    Mes_Update <- Mes[Mes %in% rownames(Exp_RNA)]
    
    # Estimate GSVA enrichment scores (ES) for two lists of EMT signatures: Epi and Mes. 
    # The method is ssGSVA (tau = 0.25 in Barbie et al, 2009, and tau = 0.75 in Verhaak et al, 2013 and Tan et al, 2014).
    gsva_es_list <- gsva(Exp_RNA, list(Epi_Update, Mes_Update), method = "ssgsea", tau = .75)
    # The difference of GSVA enrichment scores between Mes and Epi 
    EMT_score <- gsva_es_list[2, ] - gsva_es_list[1, ]
    
    ## Using two-sample KS test, compute the p-values of the difference of GSVA enrichment scores between Epi and Mes
    Epi_signature <- lapply(seq(Epi_Update), function(i) Epi_Update[i])
    Mes_signature <- lapply(seq(Mes_Update), function(i) Mes_Update[i])
    gsva_es_Epi_signature <- gsva(Exp_RNA, Epi_signature, method = "ssgsea", tau = .75)
    gsva_es_Mes_signature <- gsva(Exp_RNA, Mes_signature, method = "ssgsea", tau = .75)
    EMT_p.value <- unlist(lapply(seq(dim(Exp_RNA)[2]), function(i) ks.test(gsva_es_Epi_signature[, i], gsva_es_Mes_signature[, i])$p.value))

    # We divide the samples to four types: 
    # Mes (EMT_p.value[i] < 0.05 & EMT_score[i] > 0), 
    # Intermediate Mes (EMT_p.value[i] >= 0.05 & EMT_score[i] > 0), 
    # Epi (EMT_p.value[i] < 0.05 & EMT_score[i] < 0),
    # Intermediate Epi (EMT_p.value[i] >= 0.05 & EMT_score[i] < 0)
    Sample_Phenotype <- c()
    for(i in seq(dim(Exp_RNA)[2])){
        if(EMT_p.value[i] < 0.05 & EMT_score[i] > 0){
            Sample_Phenotype[i] = "Mes"}
        else if(EMT_p.value[i] >= 0.05 & EMT_score[i] > 0){
            Sample_Phenotype[i] = "Intermediate Mes"}
        else if(EMT_p.value[i] < 0.05 & EMT_score[i] < 0){
            Sample_Phenotype[i] = "Epi"}
        else if(EMT_p.value[i] >= 0.05 & EMT_score[i] < 0){
            Sample_Phenotype[i] = "Intermediate Epi"}
    }

return(cbind(EMT_score, EMT_p.value, Sample_Phenotype))
}

## Gene differential expression analysis between Epithelial (E) and Mesenchymal (M) using limma
# Exp: Input gene expression data, rows are samples and columns are genes.
# EMT_class: EMT classes with EMT types (Epi, Intermediate Epi, Intermediate Mes, and Mes) by using generic EMT signatures in GSVA R package.
# trend: Logical value. TRUE for log2(FPKM+1) and log2(CPM+1) values, FALSE for read counts.
# p.value: Cutoff value for adjusted p-values. Only genes with lower p-values are listed.
# lfc: Minimum absolute log2-fold-change required.
# adjust.method: Method used to adjust the p-values for multiple testing, including "none", "BH", "BY" and "holm". The default adjustment method is "BH".
limma_DEG <- function(Exp, EMT_class, trend = TRUE, p.value = 0.01, lfc = 1, adjust.method = "BH") {
    
    Exp_E <- Exp[which(EMT_class[, 3] == "Epi"), ]
    Exp_M <- Exp[which(EMT_class[, 3] == "Mes"), ]
    design <- model.matrix(~ -1+factor(c(rep(1, nrow(Exp_E)), rep(2, nrow(Exp_M)))))
    colnames(design) <- c("E", "M")
    contrast.matrix <- makeContrasts(E-M, levels = design)
    Exp_EM <- t(rbind(Exp_E, Exp_M))
    fit1 <- lmFit(Exp_EM, design)    
    fit1 <- eBayes(fit1, trend = TRUE)
    fit2 <- contrasts.fit(fit1, contrast.matrix)
    fit2 <- eBayes(fit2, trend = TRUE)
    res <- topTable(fit2, number = nrow(Exp_EM), p.value = p.value, lfc = lfc, sort.by = "p", adjust.method = adjust.method)
    Exp_DEG <- Exp[, which(colnames(Exp) %in% rownames(res))]

    return(list(res, Exp_DEG))
}

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

# Identify EMT phenotypes for each TCGA sample by using generic EMT signatures
library(GSVA)
Generic_EMT_signature <- read.csv("Generic_EMT_signature.csv", header=FALSE, sep=",")
Exp_RNA <- t(cbind(Pancancer_miRNA_Exp, Pancancer_lncRNA_Exp, Pancancer_mRNA_Exp, Pancancer_pseudogene_Exp))
Pseudo_sample <- Pseudo_EMT(Generic_EMT_signature, Exp_RNA)

# Gene differential expression analysis
library(limma)
Pancancer_miRNA_Exp_DEG <- limma_DEG(Pancancer_miRNA_Exp, Pseudo_sample, p.value = 0.01, lfc = log2(1.5), adjust.method = "BH")
Pancancer_lncRNA_Exp_DEG <- limma_DEG(Pancancer_lncRNA_Exp, Pseudo_sample, p.value = 0.01, lfc = log2(1.5), adjust.method = "BH")
Pancancer_pseudogene_Exp_DEG <- limma_DEG(Pancancer_pseudogene_Exp, Pseudo_sample, p.value = 0.01, lfc = log2(1.5), adjust.method = "BH")
Pancancer_mRNA_Exp_DEG <- limma_DEG(Pancancer_mRNA_Exp, Pseudo_sample, p.value = 0.01, lfc = log2(1.5), adjust.method = "BH")

## SPONGE method for miRNA sponge interaction networks
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

