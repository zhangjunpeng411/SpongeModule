%% The jointNMF Matlab scripts are adapted from Xiao et al. (Xiao Q, Luo J, Liang C, Cai J, Li G, Cao B. CeModule: an integrative framework for discovering regulatory patterns from genomic data in cancer. BMC Bioinformatics. 2019;20(1):67.)
%% One of the jointNMF method CeModule is available at https://github.com/xiaoqiu2018/CeModule
%% jointNMF (running on MATLAB 2010 or newer version)

Function Descriptions
1. jointNMF.m     
input:     
X1 (S,N1): S (sample) x N1 (lncRNA or pseudogene) non negative input matrix;     
X2 (S,N2): S (sample) x N2 (miRNA) non negative input matrix;     
X3 (S,N3): S (sample) x N3 (mRNA) non negative input matrix;     
A,B,C: miRNA-lncRNA/pseudogene (rows are miRNAs, cols are lncRNAs or pseudogenes),miRNA-mRNA (rows are miRNAs, cols are mRNAs)and mRNA-mRNA adjacent matrices;     
a        : control the trade-off of Hi;     
L1, L2, L3      : weight the must link constraints in A,B,C;     
r1       : limit the growth of W;     
r2       : limit the growth of H1, H2 and H3;     
K        : Number of components;     
output:     
W        : S x K matrix;     
H1       : N1 x K matrix;     
H2       : N2 x K matrix;     
H3       : N3 x K matrix;     

2. jointNMF_modules.m     
input:     
W          : common basis matrix;     
H1,H2,H3   : coefficient matrices;     
T          : a given threshold for z-score;     
lncRNAs (pseudogenes), miRNAs, mRNAs        :list of all lncRNAs (pseudogenes), miRNAs, mRNAs;          
output:     
Co_module : the index list of lncRNAs (pseudogenes), miRNAs and mRNAs in a module.    

