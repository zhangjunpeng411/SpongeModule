function [W,H1,H2,H3] = jointNMF(X1, X2, X3, A, B, C, a, L1, L2, L3, r1, r2, K)
%
% INPUT:
% X1 (S,N1): S (sample) x N1 (lncRNA) non negative input matrix
% X2 (S,N2): S (sample) x N2 (miRNA) non negative input matrix
% X3 (S,N3): S (sample) x N3 (mRNA) non negative input matrix
% A,B,C: miRNA-lncRNA,miRNA-mRNA and mRNA-mRNA adjacent matrices
% a        : control the trade-off of Hi
% L1, L2, L3      : weight the must link constraints in A,B,B
% r1       : limit the growth of W
% r2       : limit the growth of H1, H2 and H3
% K        : Number of components

% OUTPUT:
% W        : S x K matrix
% H1       : N1 x K matrix
% H2       : N2 x K matrix
% H3       : N3 x K matrix

% avoid this kind of colomn or row: sum == 0
index = find(sum(X1,1) == 0); 
X1(:,index) = X1(:,index) + eps;

index = find(sum(X2,1) == 0);
X2(:,index) = X2(:,index) + eps;

index = find(sum(X3,1) == 0);
X3(:,index) = X3(:,index) + eps;

index = find(sum(A,1) == 0);
A(:,index) = A(:,index) + eps;

index = find(sum(B,1) == 0);
B(:,index) = B(:,index) + eps;

index = find(sum(C,1) == 0);
C(:,index) = C(:,index) + eps;

nloop = 5; 
verbose=1;
n1 = size(X1,2);
n2 = size(X2,2);
[s,n3] = size(X3);

bestW=zeros(s,K);
bestH1=zeros(n1,K);
bestH2=zeros(n2,K);
bestH3=zeros(n3,K);

bestobj1=10000000;
bestobj2=10000000;
bestobj3=10000000;
fid = fopen(['record_K' int2str(K) '_a=' num2str(a) '_L1=' num2str(L1) '_L2=' num2str(L2) '_L3=' num2str(L3) '_r1=' num2str(r1) '_r2=' num2str(r2) '.txt'],'wt+');

for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 
    disp(['loop: ', num2str(iloop)]);    
    maxiter=500; %User adjustable parameters      
    speak=1; 
    print_iter = iloop

	if (min(min(X1)) < 0) || (min(min(X2)) < 0) || (min(min(X3)) < 0) 
	    error('Input matrix elements can not be negative');
	    return
	end

	[s1,n1] = size(X1);
	[s2,n2] = size(X2);
	[s3,n3] = size(X3);
	if s1 ~= s2 || s2 ~= s3 || s1 ~= s3
	    error('Input matrices should have the same rows');
	    return
	end
	s = s1;


        % initialize random W, H1, H2 and H3
	W = rand(s,K);
	H1 = rand(n1,K);
	H2 = rand(n2,K);
	H3 = rand(n3,K);

        % use W*H to test for convergence
	Xr_old1 = W*H1';
	Xr_old2 = W*H2';
	Xr_old3 = W*H3';

	for iter = 1 : maxiter
	    W = W.*([X1 X2 X3]*[H1; H2; H3])./(W*([H1; H2; H3]'*[H1; H2; H3]+r1*eye(K))+eps);    %Update rule-1;   
	    H1_new = H1.*(X1'*W + a*H1 + L1/2*A'*H2)./(H1*(W'*W + a*H1'*H1) + r2/2*ones(n1,K)+eps);
	    H2_new = H2.*(X2'*W + a*H2 + L1/2*A*H1 + L2/2*B*H3)./(H2*(W'*W + a*H2'*H2) + r2/2*ones(n2,K)+eps);
	    H3 = H3.*(X3'*W + a*H3 + L2/2*B'*H2 + L3*C*H3)./(H3*(W'*W+a*H3'*H3) + r2/2*ones(n3,K)+eps); 
	    H1 = H1_new;
	    H2 = H2_new;
	   
	    if (rem(iter,print_iter) == 0) && speak
		Xr1 = W*H1'+eps;            
		Xr2 = W*H2'+eps;
		Xr3 = W*H3'+eps;
		diff_iter = sum(sum(abs(Xr_old1-Xr1)))+sum(sum(abs(Xr_old2-Xr2)))+sum(sum(abs(Xr_old3-Xr3)));
		       
		eucl_dist1 = nmf_euclidean_dist(X1,W*H1');
		eucl_dist2 = nmf_euclidean_dist(X2,W*H2');
		eucl_dist3 = nmf_euclidean_dist(X3,W*H3');        
		diff1 = eucl_dist1 + eucl_dist2 + eucl_dist3;  
		
		diff2 = a/2*nmf_euclidean_dist(H1'*H1, eye(K))+ a/2*nmf_euclidean_dist(H2'*H2, eye(K))+ a/2*nmf_euclidean_dist(H3'*H3, eye(K));
		diff3 = -L1*trace(H2'*A*H1);
		diff4 = -L2*trace(H2'*B*H3);
		diff5 = -L3*trace(H3'*C*H3);
		diff6 = r1*sum(sum(W.*W));
		diff7 = r2*(sum(sum(H1))+sum(sum(H2))+sum(sum(H3)));
		total_obj = diff1 + diff2 + diff3 + diff4 + diff5 + diff6  + diff7;
		
		errorx1 = mean(mean(abs(X1-W*H1')))/mean(mean(X1));
		errorx2 = mean(mean(abs(X2-W*H2')))/mean(mean(X2));
		errorx3 = mean(mean(abs(X3-W*H3')))/mean(mean(X3));
		errorx = errorx1 + errorx2 + errorx3;
		
		disp(['iter = ', num2str(iter),'  relative error = ', num2str(errorx)]);
		disp(['diff_iter = ', num2str(diff_iter), '  total_obj = ', num2str(total_obj)]);

		Xr_old1 = Xr1;
		Xr_old2 = Xr2;
		Xr_old3 = Xr3;
		fprintf(fid,'%s\n',[sprintf('Iter = \t'),int2str(iter),...
		    sprintf('\t relative error = \t'),num2str(errorx),...
		    sprintf('\t diff_iter = \t'),num2str(diff_iter),...
		    sprintf('\t total_obj = \t'), num2str(total_obj),...
		    sprintf('\t diff1 = \t'), num2str(diff1),...
		    sprintf('\t diff2 = \t'), num2str(diff2),...
		    sprintf('\t diff3 = \t'), num2str(diff3),...
			    sprintf('\t diff4 = \t'), num2str(diff4),...
		    sprintf('\t diff5 = \t'), num2str(diff5),...
		    sprintf('\t diff6 = \t'), num2str(diff6),...
			    sprintf('\t diff7 = \t'), num2str(diff7)]);
		if errorx < 10^(-5), break, end
	    end
	 end

    newobj1 = sum(sum((X1-W*H1').^2));
    newobj2 = sum(sum((X2-W*H2').^2));
    newobj3 = sum(sum((X3-W*H3').^2));
       
    if (newobj1<bestobj1)||(newobj2<bestobj2)||(newobj3<bestobj3)          
        bestobj1 = newobj1;
        bestobj2 = newobj2;  
        bestobj3 = newobj3;
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
        bestH3 = H3;
    end
end
fclose(fid);
W = bestW; H1 = bestH1; H2 = bestH2; H3 = bestH3;

function err = nmf_euclidean_dist(X,Y)
err = sum(sum((X-Y).^2));

