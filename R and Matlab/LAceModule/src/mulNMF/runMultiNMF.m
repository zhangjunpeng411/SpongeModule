function runMultiNMF(COR,LA,Klist,useGPU,resultPath)
  disp(1)
    if exist(resultPath,'dir')==0
       mkdir(resultPath);
    end
	%% parameter setting
	options = [];
	options.maxIter = 200;
	options.error = 1e-6;
	options.nRepeat = 30;
	options.minIter = 30;
	options.meanFitRatio = 0.1;
	options.rounds = 30;
    options.alpha = [1 1];
    
	COR=(COR+COR')/2;
	LA=(LA+LA')/2;
	label=zeros(size(COR,1),1);
	data{1} = single(COR);
	data{2} = single(LA);

	%% normalize data matrix
    disp(2)
	for i = 1:length(data)
		data{i} = data{i} / sum(sum(data{i}));
	end
	if(useGPU)
		data{1}=gpuArray(data{1});
		data{2}=gpuArray(data{2});
	end

	%U_final = cell(1,length(Klist));
	%V_final = cell(1,length(Klist));
	%V_centroid = cell(1,length(Klist));
	for i =1:length(Klist)
		[A, B, C, log] = MultiNMF(data, Klist(i), label, options);
		if(useGPU)
			VC=gather(C);
		else
			VC=C;
		end
	  save([resultPath,'/VC-',num2str(Klist(i)),'.mat'],'VC')
	  clear('A','B','C')
	end
end
