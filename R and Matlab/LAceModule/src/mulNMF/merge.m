basedir='D:/¿ÎÌâ/liquidassociation/LIHC/repeat3';
num=10;
all=cell(1,num);
cluster=zeros(num,9404);
for i=11:20
    load([basedir,'/result',num2str(i),'-VC.mat'])
    all{i}=VC;
    [m,index]=max(VC,[],2);
    cluster(i-10,:)=index;
end
VC=all;
save([basedir,'/result1-50-VC.mat'],'VC')
save([basedir,'/cluster1-50.mat'],'cluster')

load([basedir,'/result1-50-VC.mat'])
load([basedir,'/cluster1-50.mat'])
error=zeros(num,num);
for i=1:10
    for j=1:10
        dis=VC{i}-VC{j};
        error(i,j)=sum(sum(dis.*dis));
    end
end
[m,i]=min(sum(error))
cluster=cluster(i,:);
save([basedir,'-1/result1-50-VC.mat'],'VC');
save([basedir,'-1/cluster1-50.mat'],'cluster');

load([basedir,'/result1-50-VC.mat'])
load([basedir,'/cluster1-50.mat'])
[x,y]=size(cluster);
consenseM=zeros(y,y);
for i=1:x
    c=unique(cluster(i,:));
    for k=c
        index=find(cluster(i,:)==k);
        consenseM(index,index)=consenseM(index,index)+1;
    end
end
[group, eigengap] = SpectralClustering(consenseM, 360)
cluster=group';
save([basedir,'/cluster.final.mat'],'cluster');

load([basedir,'/result1-50-VC.mat'])
result=zeros(size(VC{1}));
for i=1:num
    result=result+VC{i};
end
VC=result/num;
[m,index]=max(VC,[],2);
cluster=index';
save([basedir,'-1-2/result1-50-VC.mat'],'VC')
save([basedir,'-1-2/cluster1-50.mat'],'cluster')