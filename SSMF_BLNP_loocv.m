clc;clear;
%% 原始处理
%{
lncSim = load('DATA\lncRNAsimilarity.txt');   %lncRNA 表达相似性
disSim = load('DATA\diseasesimilarity.txt');  % disease 语义相似性
interaction = load('DATA\known_lncRNA_disease_interaction.txt');
%}

lncSim = importdata('DATA\B.mat');   %lncRNA 表达相似性
disSim = importdata('DATA\DiseaseSimilarityModel.xlsx');  % disease 语义相似性
interaction = load('DATA\DiseaseAndRNABinary.csv');
interaction = interaction';
[nl,nd] = size(interaction);
interaction_ori = interaction; 

n=1;
m=1;
[~,l] = size(lncSim);
[~,d] = size(disSim);
%% CSMF算法部分
%for alpha=0.1:0.1:0.9
 alpha = 0.1;      
 beta = 0.1;
%for beta=0.1:0.1:0.9
%for K=20:20:100
%for K3=26:24:218
%for t=1:3
K = 180;
t=1;
K3 = 194;
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
%% BLNP主方法
F_ori = BLNP(interaction_ori,KD,KL,nl,nd,beta);
F_ori0 = F_ori;
%%  留一交叉验证（LOOCV）
index=find(interaction_ori==1); 
for u=1:length(index)
     K;
     K3;
     t;
     alpha;
     u
     
    interaction(index(u))=0;  
% computing liner neighborhood interaction profile kernel of lncRNAs and disease
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
F_ori(index(u))=F(index(u));
interaction = interaction_ori;    
end
pre_label_score = F_ori(:);
label_y = interaction_ori(:);
auc=roc_1(pre_label_score,label_y,'red');

[ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.95' );
aupr=pr_cure(pre_label_score,label_y,'red');
x(m,n)=aupr;
y(m,n)=auc;
n=n+1;
%end
n=1;
m=m+1;
%end  
%end  
%grid on;



