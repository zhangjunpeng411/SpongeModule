function [Co_module] = jointNMF_modules(W, H1, H2, H3, T, lncRNAs, miRNAs, mRNAs)
%
% INPUT
% W          : common basis matrix
% H1,H2,H3   : coefficient matrices
% T          : a given threshold for z-score.
% lncRNAs, miRNAs, mRNAs        :list of all lncRNAs, miRNAs, mRNAs  
%
% OUTPUT
% Co_module : the index list of lncRNAs, miRNAs and mRNAs in a module.
%

n1 = size(H1,1);
n2 = size(H2,1);
n3 = size(H3,1);
K = size(W,2);
Co_module = cell(K,3);

MH1 = mean(H1,1);  
VH1 = std(H1,0,1);
MH2 = mean(H2,1);  
VH2 = std(H2,0,1);
MH3 = mean(H3,1);
VH3 = std(H3,0,1);

for i = 1:K
    c1 = find(H1(:,i) > MH1(i) + T*VH1(i));
    c1(c1 > n1/2) = c1(c1 > n1/2) - n1/2;% tranform the double lncRNA index into origin index
    c1 = round(unique(c1));        
      
    c2 = find(H2(:,i) > MH2(i) + T*VH2(i));
    c2(c2 > n2/2) = c2(c2 > n2/2) - n2/2; % tranform the double microRNA index into origin index
    c2 = round(unique(c2));
        
    c3 = find(H3(:,i) > MH3(i) + T*VH3(i));
    c3(c3 > n3/2) = c3(c3 > n3/2) - n3/2;% tranform the double mRNA index into origin index
    c3 = round(unique(c3));
        
    Co_module{i,1}=lncRNAs(c1); Co_module{i,2}=miRNAs(c2); Co_module{i,3}=mRNAs(c3);
end
