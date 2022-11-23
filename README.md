# SSMF-BLNP
Using computer technology to study bioinformatics. Code and Datasets of Paper--Predicting lncRNA-disease associations based on combining selective similarity matrix fusion and bidirectional linear neighborhood label propagation
Input:
B.mat (lncRNA expression similarity metrix)
DiseaseAndRNABinary.csv (lncRNA-disease associations)
DiseaseSimilarityModel.xlsx (disease semantic similarity metrix)

Output:
Predictive scoring matrix

Label_Propagation.m -- Calculating linear neighborhood similarity
SSMF.m -- Selective similarity matrix fusion
BLNP.m -- Bidirectional linear neighborhood label propagation
SSMF_BLNP_LOOCV.m -- LOOCV validation process
SSMF_BLNP_5_fold.m -- 5-fold validation process
