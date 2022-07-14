function [W]=SSMF(Wall,K,t,ALPHA)
s=0;
C = length(Wall);
[m,n]=size(Wall{1});
I=eye(m,n);

%% 选择性归一化
for i = 1 : C
    newW1{i} = Wall{i}./repmat(sum(Wall{i},1),n,1); 
    TF = isnan(newW1{i});
    a = [];
    for j=1:m
        if TF(1,j) == 1
            s = s+1;
            a(s) = j;
        end
    end
    s;
    newW1{i}(find(TF(:,:)==1)) = 0;       
end
%%
sumW1 = zeros(m,n);
%求和 P_{m,l}
for i = 1 : C
    sumW1= sumW1 + newW1{i};
end



%找邻居
for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));
end
%% 选择性归一化处理
for i = 1 : C
    Wall{i} = Wall{i}./repmat(sum(Wall{i},1),n,1);

    TF = isnan(Wall{i});
    Wall{i}(find(TF(:,:)==1)) = 0;  
   
end
%%
Wsum = zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};
end
for ITER=1:t
    for i = 1 : C           
         Wall{i}=ALPHA*newW{i}*(Wsum - Wall{i} + I)*newW{i}' + (1-ALPHA)*(sumW1 - newW1{i});
    end   
    Wsum = zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
    end     
end
%W = Wsum ./ (repmat(sum(Wsum,2),1,n));
%W = (W + W' + eye(n,m) / 2;
W = Wsum/C;
end

%% 内核邻域归一化
function newW = FindDominateSet(W,K)
[m,n]=size(W);
%按行降序排序
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
TF = isnan(newW);
newW(find(TF(:,:)==1)) = 0;
clear IW1;
clear IW2;
clear temp;
end

