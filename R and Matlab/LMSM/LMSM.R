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

# LMSM method for miRNA sponge modules
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


