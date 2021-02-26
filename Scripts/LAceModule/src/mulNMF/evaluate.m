% load cluster1-30.mat
% load Code_multiNMF\data.mat
% distance=max(max(COR))-COR;
% distance(logical(eye(size(distance))))=0;
% distance=triu(distance);
% binary=zeros(size(distance));
function [c_index,McClain_Rao,Point_biserial,Modularity]=evaluate(cluster,similar)
    distance=max(max(similar))-similar;
    distance(logical(eye(size(distance))))=0;
    distance=triu(distance);
    binary=zeros(size(distance));

    c_index=zeros(1,size(cluster,1));
    McClain_Rao=c_index;
    Point_biserial=c_index;
    Modularity=c_index;
    D=c_index;
    for i=1:size(cluster,1)
        disp(i);
        [NT,NB,NW,SW,SB,Smin,Smax]=basevalue(cluster(i,:),distance,binary);
        c_index(i)=(SW-Smin)/(Smax-Smin);
        McClain_Rao(i)=(SW/NW)/(SB/NB);
        Point_biserial(i)=(SW/NW-SB/NB)*sqrt(NW*NB)/NT;
        Modularity(i)=modularity(cluster(i,:),similar,binary);
    end
end

function [NT,NB,NW,SW,SB,Smin,Smax]=basevalue(vector,distance,binary)
    NT=length(vector)*(length(vector)-1)/2;
    count=tabulate(vector);
    count=count(:,1:2);
    NW=sum(count(:,2).*(count(:,2)-1)/2);
    NB=NT-NW;
    for c=count(:,1)'
        index=find(vector==c);
        binary(index,index)=1;
    end
    binary(logical(eye(size(binary))))=0;
    binary=triu(binary);
    SW=sum(sum(distance.*binary));
    SB=sum(sum(triu(~binary).*distance));
    t=sort(distance(:));
    Smin=sum(t(1:NW));
    t=sort(distance(:),'descend');
    Smax=sum(t(1:NW));
end

function M=modularity(vector,similar,binary)
    for c=unique(vector)
        index=find(vector==c);
        binary(index,index)=1;
    end
    binary(logical(eye(size(binary))))=0;

    rowsum=sum(similar,1);
    total=sum(sum(similar));
    M=sum(sum((similar-rowsum'*rowsum/total).*binary))/total;
end

function D=Dunn_index(vector,distance)
    cluster=unique(vector);
    distance=distance+diag(repmat(NaN,size(distance,1),1));
    dmin=[];
    dmax=[];
    for c=cluster
        index1=find(vector==c);
        index2=find(vector~=c);
        dmin=[dmin,min(min(distance(index1,index2)))];
        dmax=[dmax,max(max(distance(index1,index1)))];
    end
    D=min(dmin)/max(dmax);
end