clc
clear all
% RMA normalization
DataMatrix = affyrma('*', 'HG-U133_Plus_2.cdf','CELPath', '/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/Microarray/GEO/celFile','CDFPath', '/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/Microarray/GEO/A-AFFY-44.cdf','Output', 'log2');
Data=double(DataMatrix);
Information=get(DataMatrix);
GeneName=Information.RowNames;
SampleName=Information.ColNames;

clear DataMatrix Information

load GPL570
% match the gene names with the normalized features
[c idx idx1]=intersect(GPL570(1:41639,1),GeneName);
for i=1:41639
    DataAfterSort(idx(i,1),:)=Data(idx1(i,1),:);
end
GeneID=GPL570(1:41639,2);
clear c idx idx1 Data i GPL560 GeneName 

k=1;
n=1;
A=ones(20827,1);
for i=1:41638
    if strcmp(GeneID{i,1},GeneID{i+1,1})==1;
        k=k+1;
    else
        A(n,1)=k;
        k=1;
        GeneName{n,1}=GeneID{i,1};
        n=n+1;
        
    end
end

% generate the gene expression values
Index=A;
Data=zeros(size(Index,1),size(SampleName,2));
n=1;
for i=1:size(Index)
    if Index(i,1)==1
        Data(i,:)=DataAfterSort(n,:);
    else 
        Data(i,:)=mean(DataAfterSort(n:n+Index(i,1)-1,:));
    end
    n=n+Index(i,1);
end

clear i n DataAfterSort Index A GPL570 k GeneID

save microarray_GEO
