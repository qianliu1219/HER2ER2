---
This is the source code for our *Development and validation of a prognostic 15-gene signature for stratifying HER2+/ER+ breast cancer* paper.

![Overall workflow:](https://github.com/qianliu1219/HER2ER2/blob/main/fig/1.png "Title")
---

#### **Overall workflow:** 15,850 genes are in common among TCGA-BRCA, METABRIC, and GSE149283 HER2+/ER+ patients.  12,236 of them with at least one count in one sample are kept and input into a Cox regression-based feature selection step which results in 549 significant genes based on the criteria of p-value < 0.01. Consensus clustering are then performed to stratify TCGA-BRCA HER2+/ER+ patients based on gene expression profile of these 549 significant genes. Gene differential analysis are done among the identified subtypes to identify most differentially expressed genes. Genes that are significant in both Cox regression analysis and gene expression differential analysis are selected to form the proposed gene signature. Validations of this gene signature are performed on METABRIC and GSE149283 HER2+/ER+ cohorts. An XGBoost classifier is trained using the proposed gene signature on TCGA-BRCA data, then applied to assign METABRIC and GSE149283 BCs into two subgroups. For METABRIC, survival difference of the predicted subgroups is tested. For GSE149283, the drug response difference between the predicted subgroups is tested. 

