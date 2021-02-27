% load cluster1-30.mat
% load Code_multiNMF\data.mat
% distance=max(max(COR))-COR;
% distance(logical(eye(size(distance))))=0;
% distance=triu(distance);
% binary=zeros(size(distance));
function [Modularity,Density]=net_eva(cluster,network)
    Density=density(cluster,network);
    Modularity=modu(cluster,network);
end

function den=density(clusters, network)
    den=[];
    for i=1:size(clusters,1)
        denc=[];
        disp(i)
        vector=clusters(i,:);
        cluster=unique(vector);
        for c=cluster
            index=find(vector==c);
            n=length(index);
            dc=sum(sum(network(index,index)))/n^2;
            denc=[denc,dc];
        end
        den=[den,mean(denc)];
    end
end
function md=modu(cluster,network)
    md=[];
    binary=zeros(size(network));
    for i=1:size(cluster,1)
        disp(i)
        M=modularity(cluster(i,:),network,binary)
        md=[md,M];
    end
end