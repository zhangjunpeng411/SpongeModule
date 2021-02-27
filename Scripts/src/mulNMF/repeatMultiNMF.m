function finalC=repeatMultiNMF(COR,LA,K,useGPU,resultPath)
    %% parameter setting
    options = [];
    options.maxIter = 200;
    options.error = 1e-6;
    options.nRepeat = 30;
    options.minIter = 30;
    options.meanFitRatio = 0.1;
    options.rounds = 30;
    options.alpha = [1 1];
    %%
    COR=(COR+COR')/2;
	LA=(LA+LA')/2;
	label=zeros(size(COR,1),1);
	data{1} = single(COR);
	data{2} = single(LA);

	%% normalize data matrix
	for i = 1:length(data)
		data{i} = data{i} / sum(sum(data{i}));
	end
	if(useGPU)
		data{1}=gpuArray(data{1});
		data{2}=gpuArray(data{2});
	end


%%
    for i = 1:10
        [A, B, C, log] = MultiNMF(data, K, label, options); 
        if(useGPU)
			VC=gather(C);
		else
			VC=C;
        end
        save([resultPath,'/repeat-',num2str(i),'.mat'],'VC')
        clear('A','B','C') 
    end
    cluster=zeros(10,size(VC,1));
    for i=1:10
        load([resultPath,'/repeat-',num2str(i),'.mat']);
        [M,I]=max(VC,[],2);
        cluster(i,:)=I;
    end
    
    [x,y]=size(cluster);
    consenseM=zeros(y,y);
    for i=1:x
        c=unique(cluster(i,:));
        for k=c
            index=find(cluster(i,:)==k);
            consenseM(index,index)=consenseM(index,index)+1;
        end
    end
    [group, eigengap] = SpectralClustering(consenseM, K);
    finalC=group';
end