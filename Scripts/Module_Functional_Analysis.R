## Load required packages
library(miRSM)
library(SummarizedExperiment)
library(GSEABase)
library(survival)
library(utiml)
library(e1071)

## EMT enrichment analysis
dbEMT <- read.csv("dbEMT_v2.0.csv", header = TRUE, sep = ",")
dbEMT_SummarizedExperiment <- as.matrix(dbEMT[, 1])
dbEMT_SummarizedExperiment <- SummarizedExperiment(assays=list(dbEMT_SummarizedExperiment = dbEMT_SummarizedExperiment))

SPONGE_ceRNA_module_lncR_vs_mR_MCL_CEA <- module_CEA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
                                                     Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						     dbEMT_SummarizedExperiment, 
						     SPONGE_ceRNA_finalmodule_lncR_vs_mR_MCL)
SC_ceRNA_module_lncR_vs_mR_MCL_CEA <- module_CEA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
                                                 Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						 dbEMT_SummarizedExperiment, 
						 SC_ceRNA_finalmodule_lncR_vs_mR_MCL)
miRSM_ceRNA_module_lncR_vs_mR_WGCNA_CEA <- module_CEA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
                                                   Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						   dbEMT_SummarizedExperiment, 
						   miRSM_ceRNA_module_lncR_vs_mR_WGCNA[[2]])
CeModule_finalmodule_lncR_vs_mR_CEA <- module_CEA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
                                                   Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						   dbEMT_SummarizedExperiment,
						   CeModule_finalmodule_lncR_vs_mR)
LAceModule_finalmodule_lncR_vs_mR_CEA <- module_CEA(Pancancer_lncRNA_Exp_DEG_SummarizedExperiment, 
                                                   Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						   dbEMT_SummarizedExperiment,
						   LAceModule_finalmodule_lncR_vs_mR)

SPONGE_ceRNA_module_pseudo_vs_mR_MCL_CEA <- module_CEA(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, 
                                                     Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						     dbEMT_SummarizedExperiment, 
						     SPONGE_ceRNA_finalmodule_pseudo_vs_mR_MCL)
SC_ceRNA_module_pseudo_vs_mR_MCL_CEA <- module_CEA(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, 
                                                 Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						 dbEMT_SummarizedExperiment, 
						 SC_ceRNA_finalmodule_pseudo_vs_mR_MCL)
miRSM_ceRNA_module_pseudo_vs_mR_WGCNA_CEA <- module_CEA(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, 
                                                   Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						   dbEMT_SummarizedExperiment, 
						   miRSM_ceRNA_module_pseudo_vs_mR_WGCNA[[2]])
CeModule_finalmodule_pseudo_vs_mR_CEA <- module_CEA(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, 
                                                   Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						   dbEMT_SummarizedExperiment,
						   CeModule_finalmodule_pseudo_vs_mR)
LAceModule_finalmodule_pseudo_vs_mR_CEA <- module_CEA(Pancancer_pseudogene_Exp_DEG_SummarizedExperiment, 
                                                   Pancancer_mRNA_Exp_DEG_SummarizedExperiment, 
						   dbEMT_SummarizedExperiment,
						   LAceModule_finalmodule_pseudo_vs_mR)

## Disease enrichment analysis
SPONGE_ceRNA_module_lncR_vs_mR_MCL_DEA <- module_FA(SPONGE_ceRNA_finalmodule_lncR_vs_mR_MCL, Analysis.type = 'DEA')
SC_ceRNA_module_lncR_vs_mR_MCL_DEA <- module_FA(SC_ceRNA_finalmodule_lncR_vs_mR_MCL, Analysis.type = 'DEA')
miRSM_ceRNA_module_lncR_vs_mR_WGCNA_DEA <- module_FA(miRSM_ceRNA_module_lncR_vs_mR_WGCNA[[2]], Analysis.type = 'DEA')
CeModule_finalmodule_lncR_vs_mR_DEA <- module_FA(CeModule_finalmodule_lncR_vs_mR, Analysis.type = 'DEA')
LAceModule_finalmodule_lncR_vs_mR_DEA <- module_FA(LAceModule_finalmodule_lncR_vs_mR, Analysis.type = 'DEA')

SPONGE_ceRNA_module_pseudo_vs_mR_MCL_DEA <- module_FA(SPONGE_ceRNA_finalmodule_pseudo_vs_mR_MCL, Analysis.type = 'DEA')
SC_ceRNA_module_pseudo_vs_mR_MCL_DEA <- module_FA(SC_ceRNA_finalmodule_pseudo_vs_mR_MCL, Analysis.type = 'DEA')
miRSM_ceRNA_module_pseudo_vs_mR_WGCNA_DEA <- module_FA(miRSM_ceRNA_module_pseudo_vs_mR_WGCNA[[2]], Analysis.type = 'DEA')
CeModule_finalmodule_pseudo_vs_mR_DEA <- module_FA(CeModule_finalmodule_pseudo_vs_mR, Analysis.type = 'DEA')
LAceModule_finalmodule_pseudo_vs_mR_DEA <- module_FA(LAceModule_finalmodule_pseudo_vs_mR, Analysis.type = 'DEA')

## Survival analysis
moduleSurvival <- function(Modulelist, ExpData, SurvData, devidePercentage = 0.5, plot = FALSE) {

    ExpDataNames <- colnames(ExpData)    
    myfit <- list()
    LogRank <- list()

    for (i in seq_along(Modulelist)) {
        Interin_Data <- cbind(SurvData[, seq(2, 3)], ExpData[, which(ExpDataNames %in% Modulelist[[i]])])
        Interin_Data <- na.omit(Interin_Data)

        try_mm <- try(coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data)),
            silent = TRUE)
        if ("try-error" %in% class(try_mm))
            next

        mm <- coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data))

        Risk_score <- predict(mm, newdata = data.frame(Interin_Data), type = "risk")

        group <- rep("NA", dim(Interin_Data)[1])
        group[Risk_score > quantile(Risk_score, probs = devidePercentage)] <- "High"
        group[Risk_score <= quantile(Risk_score, probs = devidePercentage)] <- "Low"

        Data <- cbind(Interin_Data[, seq_len(2)], group)
        myfit[[i]] <- survfit(survival::Surv(time, status) ~ group, data = Data)

        sdf <- survdiff(survival::Surv(time, status) ~ group, data = Data)
        sdf.p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        HR <- (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])
        HRlow95 <- exp(log(HR) - qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))
        HRup95 <- exp(log(HR) + qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))

        LogRank[[i]] <- c(sdf$chisq, sdf.p.val, HR, HRlow95, HRup95)
    }

    if (plot) {
        for (i in seq_along(myfit)) {
            if (!is.null(LogRank[[i]])) {
                dev.new()
                plot(myfit[[i]], lty = 1, col = c("red", "green"), main = paste("Module", i), xlab = "Time (Months)",
                    ylab = "Probability of survival")

                legend("topright", legend = c("High risk group", "Low risk group"), lty = seq_len(2),
                    col = c("red", "green"))
            }
        }
    }

    LogRank_res <- do.call(rbind, LogRank)

    if (length(myfit) >= 1) {
        colnames(LogRank_res) <- c("Chi-square", "p-value", "HR", "HRlow95", "HRup95")
        names(LogRank) <- seq_along(Modulelist)
        LogRank[sapply(LogRank, is.null)] <- NULL
        rownames(LogRank_res) <- paste("Module", names(LogRank))
    }

    return(LogRank_res)
}

ExpData1 <- cbind(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]])
SPONGE_ceRNA_module_lncR_vs_mR_MCL_sur <- moduleSurvival(SPONGE_ceRNA_finalmodule_lncR_vs_mR_MCL, ExpData1, Pancancer_survival)
SC_ceRNA_module_lncR_vs_mR_MCL_sur <- moduleSurvival(SC_ceRNA_finalmodule_lncR_vs_mR_MCL, ExpData1, Pancancer_survival)
miRSM_ceRNA_module_lncR_vs_mR_WGCNA_sur <- moduleSurvival(miRSM_ceRNA_module_lncR_vs_mR_WGCNA[[2]], ExpData1, Pancancer_survival)
CeModule_finalmodule_lncR_vs_mR_sur <- moduleSurvival(CeModule_finalmodule_lncR_vs_mR, ExpData1, Pancancer_survival)
LAceModule_finalmodule_lncR_vs_mR_sur <- moduleSurvival(LAceModule_finalmodule_lncR_vs_mR, ExpData1, Pancancer_survival)

ExpData2 <- cbind(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]])
SPONGE_ceRNA_module_pseudo_vs_mR_MCL_sur <- moduleSurvival(SPONGE_ceRNA_finalmodule_pseudo_vs_mR_MCL, ExpData2, Pancancer_survival)
SC_ceRNA_module_pseudo_vs_mR_MCL_sur <- moduleSurvival(SC_ceRNA_finalmodule_pseudo_vs_mR_MCL, ExpData2, Pancancer_survival)
miRSM_ceRNA_module_pseudo_vs_mR_WGCNA_sur <- moduleSurvival(miRSM_ceRNA_module_pseudo_vs_mR_WGCNA[[2]], ExpData2, Pancancer_survival)
CeModule_finalmodule_pseudo_vs_mR_sur <- moduleSurvival(CeModule_finalmodule_pseudo_vs_mR, ExpData2, Pancancer_survival)
LAceModule_finalmodule_pseudo_vs_mR_sur <- moduleSurvival(LAceModule_finalmodule_pseudo_vs_mR, ExpData2, Pancancer_survival)

## Classification analysis
# Evaluate the performance of each module for classifying EMT types
module.classify <- function(ceRExp, mRExp, EMT_type, Modulelist, method = "br", base.algorith = "SVM", cv.folds = 10, 
	                    cv.sampling = "stratified", cv.seed = 123) {

    module_ceRExp <- lapply(seq_along(Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]])])
    module_mRExp <- lapply(seq_along(Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist[[i]])])
    Epi <- as.numeric(EMT_type[, 3] == "Epi")
    Mes <- as.numeric(EMT_type[, 3] == "Mes")
    InterEpi <- as.numeric(EMT_type[, 3] == "Intermediate Epi")
    InterMes <- as.numeric(EMT_type[, 3] == "Intermediate Mes")        
    module_classify <- list()

    for (i in seq_along(Modulelist)){        
	
        temp <- as.data.frame(cbind(module_ceRExp[[i]], module_mRExp[[i]], Epi, Mes, InterEpi, InterMes))
	Indices <- ncol(temp)
	temp_mldr <- mldr_from_dataframe(temp, labelIndices = c(Indices-3, Indices-2, Indices-1, Indices), name = "TEMPMLDR")
        temp_res <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds, 
	               cv.sampling = cv.sampling, cv.seed = cv.seed)
        module_classify[[i]] <- temp_res

    }

    return(module_classify)
}

SPONGE_ceRNA_module_lncR_vs_mR_MCL_classify <- module.classify(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, SPONGE_ceRNA_finalmodule_lncR_vs_mR_MCL)
SC_ceRNA_module_lncR_vs_mR_MCL_classify <- module.classify(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, SC_ceRNA_finalmodule_lncR_vs_mR_MCL)
miRSM_ceRNA_module_lncR_vs_mR_WGCNA_classify <- module.classify(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, miRSM_ceRNA_module_lncR_vs_mR_WGCNA[[2]])
CeModule_finalmodule_lncR_vs_mR_classify <- module.classify(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, CeModule_finalmodule_lncR_vs_mR)
LAceModule_finalmodule_lncR_vs_mR_classify <- module.classify(Pancancer_lncRNA_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, LAceModule_finalmodule_lncR_vs_mR)  

SPONGE_ceRNA_module_pseudo_vs_mR_MCL_classify <- module.classify(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, SPONGE_ceRNA_finalmodule_pseudo_vs_mR_MCL)
SC_ceRNA_module_pseudo_vs_mR_MCL_classify <- module.classify(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, SC_ceRNA_finalmodule_pseudo_vs_mR_MCL)
miRSM_ceRNA_module_pseudo_vs_mR_WGCNA_classify <- module.classify(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, miRSM_ceRNA_module_pseudo_vs_mR_WGCNA[[2]])
CeModule_finalmodule_pseudo_vs_mR_classify <- module.classify(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, CeModule_finalmodule_pseudo_vs_mR)
LAceModule_finalmodule_pseudo_vs_mR_classify <- module.classify(Pancancer_pseudogene_Exp_DEG[[2]], Pancancer_mRNA_Exp_DEG[[2]], Pseudo_sample, LAceModule_finalmodule_pseudo_vs_mR)  
