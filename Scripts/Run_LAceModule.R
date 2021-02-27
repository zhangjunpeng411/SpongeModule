%% The scripts are adapted from Xiao et al. (Wen X, Gao L, Hu Y. LAceModule: Identification of Competing Endogenous RNA Modules by Integrating Dynamic Correlation. Front Genet. 2020;11:235.)
%% LAceModule is available at https://github.com/GaoLabXDU/LAceModule

############################################################################## Note: Following codes are running on R ###########################################################################################################
# Loding Data and Codes, LAceModule_data.RData can be downloaded from https://drive.google.com/drive/folders/1VDWq7y_9PTV1cTECHTJ_bY-ZoeLDv_Bs?usp=sharing
source('Scripts/src/function.R')
load('LAceModule_data.RData')
# Set Prameters
core <- 2
Klist <- seq(2,10)
geneset_lncR_vs_mR <- rownames(target.exp_lncR_vs_mR)
geneset_pseudo_vs_mR <- rownames(target.exp_pseudo_vs_mR)
useGPU <- F
resultpath_lncR_vs_mR <- paste(getwd(),'/tempresult_lncR_vs_mR/',sep="")
resultpath_pseudo_vs_mR <- paste(getwd(),'/tempresult_pseudo_vs_mR/',sep="")

# Calculation MS, PCC and LA
mirwalk <- edge_lncR_vs_mR
MS_lncR_vs_mR <- gene_pair(geneset = geneset_lncR_vs_mR, mirwalk = mirwalk, core = core)
PCC_lncR_vs_mR <- gene_correlation(geneset = geneset_lncR_vs_mR, rna.exp = target.exp_lncR_vs_mR)
PCC_lncR_vs_mR <- PCC_lncR_vs_mR$correlation
LA_lncR_vs_mR <- liquid_association(geneset = geneset_lncR_vs_mR, rna.exp = target.exp_lncR_vs_mR, miRNA.exp = micro.exp_lncR_vs_mR, mirwalk = mirwalk, core = core)

mirwalk <- edge_pseudo_vs_mR
MS_pseudo_vs_mR <- gene_pair(geneset = geneset_pseudo_vs_mR, mirwalk = mirwalk, core = core)
PCC_pseudo_vs_mR <- gene_correlation(geneset = geneset_pseudo_vs_mR, rna.exp = target.exp_pseudo_vs_mR)
PCC_pseudo_vs_mR <- PCC_pseudo_vs_mR$correlation
LA_pseudo_vs_mR <- liquid_association(geneset = geneset_pseudo_vs_mR, rna.exp = target.exp_pseudo_vs_mR, miRNA.exp = micro.exp_pseudo_vs_mR, mirwalk = mirwalk, core = 1)

# Remove Invalid Paris
PCC_lncR_vs_mR[MS_lncR_vs_mR>=0.05] <- 0
LA_lncR_vs_mR[MS_lncR_vs_mR>=0.05] <- 0
diag(PCC_lncR_vs_mR) <- 0
diag(PCC_lncR_vs_mR) <- 0

PCC_lncR_vs_mR[PCC_lncR_vs_mR<0] <- 0
LA_lncR_vs_mR[LA_lncR_vs_mR<0] <- 0

PCC_lncR_vs_mR[PCC_lncR_vs_mR==0|LA_lncR_vs_mR==0] <- 0
LA_lncR_vs_mR[PCC_lncR_vs_mR==0|LA_lncR_vs_mR==0] <- 0

PCC_pseudo_vs_mR[MS_pseudo_vs_mR>=0.05] <- 0
LA_pseudo_vs_mR[MS_pseudo_vs_mR>=0.05] <- 0
diag(PCC_pseudo_vs_mR) <- 0
diag(PCC_pseudo_vs_mR) <- 0

PCC_pseudo_vs_mR[PCC_pseudo_vs_mR<0] <- 0
LA_pseudo_vs_mR[LA_pseudo_vs_mR<0] <- 0

PCC_pseudo_vs_mR[PCC_pseudo_vs_mR==0|LA_pseudo_vs_mR==0] <- 0
LA_pseudo_vs_mR[PCC_pseudo_vs_mR==0|LA_pseudo_vs_mR==0] <- 0

# Decide Final K
finalK <- 100
library(R.matlab)
writeMat(LA_lncR_vs_mR=LA_lncR_vs_mR,PCC_lncR_vs_mR=PCC_lncR_vs_mR,Klist=Klist, useGPU=useGPU, resultpath_lncR_vs_mR=resultpath_lncR_vs_mR, finalK=finalK, con = 'LAceModule_dataset_lncR_vs_mR.mat')
writeMat(LA_pseudo_vs_mR=LA_pseudo_vs_mR,PCC_pseudo_vs_mR=PCC_pseudo_vs_mR,Klist=Klist, useGPU=useGPU, resultpath_pseudo_vs_mR=resultpath_pseudo_vs_mR, finalK=finalK, con = 'LAceModule_dataset_pseudo_vs_mR.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Note: Following codes are running on Matlab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set matlab code path
addpath('Scripts/src/mulNMF/');
% Detect LAceModule with candidate Klist for one time
load('LAceModule_dataset_lncR_vs_mR.mat');
load('LAceModule_dataset_pseudo_vs_mR.mat');
runMultiNMF(PCC_lncR_vs_mR,LA_lncR_vs_mR,Klist,useGPU,resultpath_lncR_vs_mR);
runMultiNMF(PCC_pseudo_vs_mR,LA_pseudo_vs_mR,Klist,useGPU,resultpath_pseudo_vs_mR);
% Detect LAceModule with finalK for 10 times and get the consensed module
finalC_lncR_vs_mR=repeatMultiNMF(PCC_lncR_vs_mR,LA_lncR_vs_mR,finalK,useGPU,resultpath_lncR_vs_mR);
xlswrite('finalC_lncR_vs_mR.xlsx',finalC_lncR_vs_mR);

finalC_pseudo_vs_mR=repeatMultiNMF(PCC_pseudo_vs_mR,LA_pseudo_vs_mR,finalK,useGPU,resultpath_pseudo_vs_mR);
xlswrite('finalC_pseudo_vs_mR.xlsx',finalC_pseudo_vs_mR);

############################################################################## Note: Following codes are running on R ###########################################################################################################
# Extracting candidate module genes
CandModgenes <- function(ceRExp, mRExp, Modulegenes, num.ModuleceRs = 2, 
    num.ModulemRs = 2){
  
    ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(ceRExp))))
    mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(mRExp))))

    index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
    CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])
    CandidateModulegenes <- lapply(seq_along(index), function(i) GeneSet(CandidateModulegenes[[i]], 
        setName = paste("Module", i, sep=" ")))
    CandidateModulegenes <- GeneSetCollection(CandidateModulegenes)
    
    return(CandidateModulegenes)
}

library(GSEABase)
finalC_lncR_vs_mR <- read.csv("finalC_lncR_vs_mR.csv", header = FALSE, sep = ",") # Final module
colnames(finalC_lncR_vs_mR) <- rownames(target.exp_lncR_vs_mR)

finalC_pseudo_vs_mR <- read.csv("finalC_pseudo_vs_mR.csv", header = FALSE, sep = ",") # Final module
colnames(finalC_pseudo_vs_mR) <- rownames(target.exp_pseudo_vs_mR)

LAceModule_module_lncR_vs_mR <- lapply(seq(finalK), function(i) colnames(finalC_lncR_vs_mR)[which(finalC_lncR_vs_mR == i)])
LAceModule_module_pseudo_vs_mR <- lapply(seq(finalK), function(i) colnames(finalC_pseudo_vs_mR)[which(finalC_pseudo_vs_mR == i)])

LAceModule_finalmodule_lncR_vs_mR <- CandModgenes(lncRNA.exp_lncR_vs_mR, mRNA.exp_lncR_vs_mR, LAceModule_module_lncR_vs_mR)
LAceModule_finalmodule_pseudo_vs_mR <- CandModgenes(pseudo.exp_pseudo_vs_mR, mRNA.exp_pseudo_vs_mR, LAceModule_module_pseudo_vs_mR)
LAceModule_finalmodule_lncR_vs_mR <- geneIds(LAceModule_finalmodule_lncR_vs_mR)
LAceModule_finalmodule_pseudo_vs_mR <- geneIds(LAceModule_finalmodule_pseudo_vs_mR)
