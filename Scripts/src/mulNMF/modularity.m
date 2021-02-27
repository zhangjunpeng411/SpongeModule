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

