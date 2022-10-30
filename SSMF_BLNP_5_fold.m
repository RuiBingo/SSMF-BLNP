clc;clear;
tic
% format long
%% 原始数据
%y = importdata('5-cv(9x9).mat');
%y1=importdata('156-158.mat');
%y2=importdata('159-161.mat');
%y3=importdata('162-164.mat');
%y=[y1;y2;y3];
%show(y);
lncSim = load('DATA\B.mat');   %lncRNA 表达相似性
disSim = importdata('DATA\DiseaseSimilarityModel.xlsx');  % disease 语义相似性
interaction = load('DATA\DiseaseAndRNABinary.csv');
interaction = interaction';
[nl,nd] = size(interaction);
interaction_ori = interaction;
%% 数据处理方法
lncSim = lncSim.B;
 n = 1;
 m=1;
 %for alpha=0:0.1:0.9
        alpha =0.1;
        beta = 0.1;

    K = 180;
    K3 = 194;

    t=1;
    km = Label_Propagation(interaction_ori,0,K,'regulation2');
    kd = Label_Propagation(interaction_ori',0,K3,'regulation2'); 

 K1 = [];
 K1(:,:,1)=km;
 K1(:,:,2)=lncSim;
 
 K2 = [];
 K2(:,:,1)=kd;
 K2(:,:,2)=disSim;
 
 KL=SSMF({K1(:,:,1),K1(:,:,2)},K,t,alpha);
 KD=SSMF({K2(:,:,1),K2(:,:,2)},K3,t,alpha);
   
%% 主方法
   F_ori = BLNP(interaction_ori,KD,KL,nl,nd,beta);
 
 %%
 F_ori_ori= F_ori;
 index=find(interaction_ori==1);
%%  5-fold 交叉验证
auc = zeros(1,100);
for k = 1:100
    K;
    K3;
    alpha
    k
    indices = crossvalind('Kfold', length(index), 5 );
    interaction = interaction_ori;
    F_ori=F_ori_ori;
for cv = 1:5 
       cv;
       index_2 = find(cv == indices);
       %%%移除已知关系
       interaction(index(index_2)) = 0;
       %%%计算得分矩阵
%%  数据处理
 km = Label_Propagation(interaction_ori,0,K,'regulation2');
 kd = Label_Propagation(interaction_ori',0,K3,'regulation2'); 

K1 = [];
 K1(:,:,1)=km;
 K1(:,:,2)=lncSim;
 
 K2 = [];
 K2(:,:,1)=kd;
 K2(:,:,2)=disSim;
 
 KL=SSMF({K1(:,:,1),K1(:,:,2)},K,t,alpha);
 KD=SSMF({K2(:,:,1),K2(:,:,2)},K3,t,alpha);

%% 主方法
   F = BLNP(interaction,KD,KL,nl,nd,beta); 
 
%%
      F_ori(index(index_2)) = F(index(index_2));
      interaction = interaction_ori;
end
%% 画auc曲线 
    pre_label_score = F_ori(:);
    label_y = interaction_ori(:);
    auc(k) = roc_1(pre_label_score,label_y,'red');
    a=1;
end
%% 
 auc_ave = mean(auc);
 auc_std = std(auc);
 

   y(m,n)=auc_ave;
   n=n+1;
 %end 
 m=m+1;
 n=1;
 %end 
 o=toc;
 o
