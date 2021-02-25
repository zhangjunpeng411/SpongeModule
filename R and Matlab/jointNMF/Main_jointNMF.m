%% The jointNMF Matlab scripts are adapted from Xiao et al. (Xiao Q, Luo J, Liang C, Cai J, Li G, Cao B. CeModule: an integrative framework for discovering regulatory patterns from genomic data in cancer. BMC Bioinformatics. 2019;20(1):67.)
%% One of the jointNMF method CeModule is available at https://github.com/xiaoqiu2018/CeModule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Note: Following codes are running on Matlab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%% dataset_lncR_vs_mR.mat and dataset_pseudo_vs_mR.mat can be downloaded from https://drive.google.com/drive/folders/1VDWq7y_9PTV1cTECHTJ_bY-ZoeLDv_Bs?usp=sharing
load('dataset_lncR_vs_mR.mat');
load('dataset_pseudo_vs_mR.mat');

%% Generate W,H1,H2,H3
[W_lncR_vs_mR,H1_lncR_vs_mR,H2_lncR_vs_mR,H3_lncR_vs_mR] = jointNMF(lncRNA_exp_lncR_vs_mR, micro_exp_lncR_vs_mR, mRNA_exp_lncR_vs_mR, miR2lncR_lncR_vs_mR, miR2mR_lncR_vs_mR, mR2mR_lncR_vs_mR, 0.0001, 0.01, 0.01, 0.01, 10, 10, 100);
[W_pseudo_vs_mR,H1_pseudo_vs_mR,H2_pseudo_vs_mR,H3_pseudo_vs_mR] = jointNMF(pseudo_exp_pseudo_vs_mR, micro_exp_pseudo_vs_mR, mRNA_exp_pseudo_vs_mR, miR2pseudo_pseudo_vs_mR, miR2mR_pseudo_vs_mR, mR2mR_pseudo_vs_mR, 0.0001, 0.01, 0.01, 0.01, 10, 10, 100);

%% Generate modules
mu = 0;
sigma = 1;
p = 0.05;
T = -icdf('Normal', p, mu, sigma);
Co_module_lncR_vs_mR = jointNMF_modules(W_lncR_vs_mR,H1_lncR_vs_mR,H2_lncR_vs_mR,H3_lncR_vs_mR, T, lncR_name_lncR_vs_mR, miR_name_lncR_vs_mR, mR_name_lncR_vs_mR);
Co_module_pseudo_vs_mR = jointNMF_modules(W_pseudo_vs_mR,H1_pseudo_vs_mR,H2_pseudo_vs_mR,H3_pseudo_vs_mR, T, pseudo_name_pseudo_vs_mR, miR_name_pseudo_vs_mR, mR_name_pseudo_vs_mR);

for i =1:100
    Module_lncR_vs_mR{i} = [Co_module_lncR_vs_mR{i,1}; Co_module_lncR_vs_mR{i,3}];
end

for i =1:100
    basename = ['Module_', num2str(i), '_lncR_vs_mR.xlsx'];
    xlswrite(basename, Module_lncR_vs_mR{1,i})
end

for i =1:100
    Module_pseudo_vs_mR{i} = [Co_module_pseudo_vs_mR{i,1}; Co_module_pseudo_vs_mR{i,3}];
end

for i =1:100
    basename = ['Module_', num2str(i), '_pseudo_vs_mR.xlsx'];
    xlswrite(basename, Module_pseudo_vs_mR{1,i})
end

############################################################################## Note: Following codes are running on R ###########################################################################################################
## Output lncRNA and pseudogene related miRNA sponge modules identified by jointNMF
library(xlsx)
jointNMF_module_lncR_vs_mR <- lapply(seq(100), function(i) c(read.xlsx(paste("Module_", i, "_lncR_vs_mR.xlsx", sep = ""), sheetName = "Sheet1", encoding = 'UTF-8', header = FALSE)))
jointNMF_module_pseudo_vs_mR <- lapply(seq(100), function(i) c(read.xlsx(paste("Module_", i, "_pseudo_vs_mR.xlsx", sep = ""), sheetName = "Sheet1", encoding = 'UTF-8', header = FALSE)))
 
jointNMF_finalmodule_lncR_vs_mR <- lapply(seq(jointNMF_module_lncR_vs_mR), function(i) unique(jointNMF_module_lncR_vs_mR[[i]][[1]]))
jointNMF_finalmodule_pseudo_vs_mR <- lapply(seq(jointNMF_module_pseudo_vs_mR), function(i) unique(jointNMF_module_pseudo_vs_mR[[i]][[1]]))
