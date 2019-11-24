clc
clear
load AMatrix;
load Counts_Normalized Count Probe CellLine;

G = [];
T = [];
P = [];
D = [];

for i = 1:length(Data)
    probe = Data(i).probe;
    P = cat(1,P,probe);
    
    tmp1 = Data(i).transcript';
    tmp2 = Data(i).ix;
    iix = strfind(tmp1,'|');
    TT = [];
    IX = [];
    for j = 1:length(iix)
        if ~isempty(iix{j,1})
            for k = 1:length(iix{j,1})
                if k == 1
                    TT = cat(1,TT,{tmp1{j,1}(1:iix{j,1}(k)-1)});
                    IX = cat(2,IX,tmp2(:,j));
                else
                    TT = cat(1,TT,{tmp1{j,1}(iix{j,1}(k-1)+1:iix{j,1}(k)-1)});
                    IX = cat(2,IX,tmp2(:,j));
                end
                
            end
            TT = cat(1,TT,{tmp1{j,1}(iix{j,1}(end)+1:end)});
            IX = cat(2,IX,tmp2(:,j));
        else
            TT = cat(1,TT,tmp1(j,1));
            IX = cat(2,IX,tmp2(:,j));
        end
    end
    

    T = cat(1,T,TT);
    
    G = cat(1,G,repmat({Data(i).gene},length(TT),1));
    
    tmp = [];
    X = [];
    for j = 1:length(probe)
        idx = find(strcmp(Probe,probe{j,1}));
        tmp = cat(1,tmp,Count(idx,:));
    end
    n = length(TT);
    for j = 1:size(tmp,2)
        y = tmp(:,j);
        
        cvx_begin quiet
            cvx_solver sedumi
            cvx_precision('best')
            variable x(n)
            minimize(norm(y-IX*x,'fro'));
            subject to
            x>=0
        cvx_end
        X = cat(2,X,x);
    end
    D = cat(1,D,X);
    
end

save IsoformExpression D G T CellLine P;

