rm(list = ls())

library(CancerSubtypes)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
library(limma)
library(survival)
library(survminer)
library(lolR)
library(ggplot2)
library(MASS)
library(ggrepel)
library(edgeR)
library(ggfortify)
library(sva)
library(caret)
library(RTCGA)
library(broom)
library(rms)
library(dplyr)
library(ComplexHeatmap)

####################################   Read data  ###################################################

## set data working directory
setwd("~/Desktop/projects/HER2+ER2+Subyping/data")

## read gene expression data
pp_TC_gxp<-read.csv("TC_gxp_pp.csv",row.names=1) #55368 genes, 123 patients 
pp_MB_gxp<-read.csv("MB_gxp_pp.csv",row.names = 1) #24368 genes, 104 patients

## match TCGA and METABRIC genes
pp_MB_gxp<-pp_MB_gxp[intersect(row.names(pp_MB_gxp),row.names(pp_TC_gxp)),] #15850 matched genes, 123 patients
pp_TC_gxp<-pp_TC_gxp[intersect(row.names(pp_MB_gxp),row.names(pp_TC_gxp)),] #15850 matched genes, 104 patients

## read clinical data
pp_TC_clinic<-read.csv("TC_clinic.csv",row.names = 1)
pp_MB_clinic<-read.csv("MB_clinic.csv",row.names = 1)

## read survial data
tc.surv.data<-read.csv("TC_survival.csv",row.names = 1)
mb.surv.data<-read.csv("MB_survival.csv",row.names = 1)


###############################################################################################################

#####################################      Consensus clustering on TCGA    ###########################################

## CancerSybtype feature selection on TCGA gene expression data
dge <- DGEList(counts=pp_TC_gxp) # 15850   123
keep <- filterByExpr(dge, design = NULL)
### remove genes that consistently have 0 counts
dge <- dge[keep,,keep.lib.sizes=F] # 12236   123
dge <- calcNormFactors(dge) # normalization factor
logCPM_cc <- cpm(dge, log=TRUE, prior.count=3) #12236   123
data.checkDistribution(logCPM_cc)
### remove genes that are not significant (p>0.01 in Cox regression analysis)
logCPM_cox<-FSbyCox(as.matrix(logCPM_cc),tc.surv.data$T,tc.surv.data$vital_status,cutoff=0.01) #549 123
data.checkDistribution(logCPM_cox)
#write.csv(logCPM_cox,"Supplementary Table 1.csv") # save the preprocessed TCGA gene expression data

## run CNMF subtyping on preprocessed TCGA gene expression data
set.seed(100)
result<-ExecuteCNMF(logCPM_cox,clusterNum=2,nrun=30)
sil<-silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
plot(sil)

## survival analysis on the clustered groups of TCGA
TC_group<-result$group
distanceMatrix<-result$distanceMatrix
p_value<-survAnalysis(mainTitle="TCGA",tc.surv.data$T,tc.surv.data$vital_status,TC_group,
                     distanceMatrix,similarity=TRUE)
names(TC_group)<-row.names(tc.surv.data)
print(tc.surv.data[1:5,])
### another way to visualize the survival results
KM <- survfit(Surv(T, vital_status) ~ TC_group, data = as.data.frame(tc.surv.data))
see<-ggsurvplot(KM, risk.table = TRUE, pval = TRUE, main="TCGA 549 genes",
                xlab="Survival time (days)", ylab= "Survival probability",
                legend.labs = c('Subgroup 1','Subgroup 2'))
print(see)




## check clinical relevance of the TCGA subgroups
### age
summary(pp_TC_clinic$age_at_diagnosis[TC_group==1])
sd(pp_TC_clinic$age_at_diagnosis[TC_group==1])
summary(pp_TC_clinic$age_at_diagnosis[TC_group==2])
sd(pp_TC_clinic$age_at_diagnosis[TC_group==2])
t.test(pp_TC_clinic$age_at_diagnosis[TC_group==1],pp_TC_clinic$age_at_diagnosis[TC_group==2])
### T
table(pp_TC_clinic$ajcc_tumor_pathologic_pt[TC_group==1])
table(pp_TC_clinic$ajcc_tumor_pathologic_pt[TC_group==2])
chisq.test(c(14,41,7,1),c(10,39,7,4)) # T1(T1,T1b,T1c); T2; T3; T4
### N
table(pp_TC_clinic$ajcc_nodes_pathologic_pn[TC_group==1])
table(pp_TC_clinic$ajcc_nodes_pathologic_pn[TC_group==2])
chisq.test(c(31,23,6,3,0),c(23,19,11,6,1)) # N0(N0,N0(i-),N0(i+)); N1(N1,N1a,N1b); N2(N2,N2a); N3(N3,N3a); NX
### M
table(pp_TC_clinic$ajcc_metastasis_pathologic_pm[TC_group==1])
table(pp_TC_clinic$ajcc_metastasis_pathologic_pm[TC_group==2])
chisq.test(c(50,1,12),c(50,1,9)) # M0(cM0(i+),M0); M1, MX
### stage
table(pp_TC_clinic$ajcc_pathologic_tumor_stage[TC_group==1])
table(pp_TC_clinic$ajcc_pathologic_tumor_stage[TC_group==2])
chisq.test(c(7,43,12,1,0),c(9,28,21,1,1)) # I(I,IA); II(II,IIA,IIB); III(IIIA,IIIB,IIIC); IV; X
### surgery
table(pp_TC_clinic$surgical_procedure_first[TC_group==1])
table(pp_TC_clinic$surgical_procedure_first[TC_group==2])
chisq.test(c(8,17,12,19,7),c(4,20,4,28,4)) #Lumpectomy; Modified Radical Mastectomy; Simple Mastectomy; Other; Not Available

### Pam50 subtype
TC_Pam50<-read.csv("TC_Pam50.csv",row.names=1)
row.names(TC_Pam50)<-gsub('-','.',row.names(TC_Pam50))
TC_Pam50<-TC_Pam50[names(TC_group),]
row.names(TC_Pam50)<-names(TC_group)
TC_Pam50$TC_group<-TC_group
table(TC_Pam50$PAM50[TC_group==1])
table(TC_Pam50$PAM50[TC_group==2])



## plot the heatmap of TCGA subgroups
TC_ht<-rbind(TC_group,logCPM_cox)
TC_ht_sort<-TC_ht[,order(TC_ht[1,])]
col<-ifelse(TC_ht_sort[1,]==1,"#F8766D","#00BFC4")
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
hm <- heatmap.3(TC_ht_sort[-1,],scale = "row",dendrogram = "none",
                key=T,keysize=2,hclustfun=myclust, distfun=mydist,
                ColSideColors=as.matrix(col),
                ColSideColorsSize=1,Colv = F,Rowv = T,
                col = colorRampPalette(c("black","darkblue","blue3","blue2","blue1",
                                         "yellow","red","red1","red2","red3","darkred")))

###############################################################################################################

#####################################      Differential analysis on TCGA     ###########################################

## read limma data
limma.data <- read.csv('TC_limma_data.csv',row.names = 1) #54247   236
subgroup <- read.csv("TC_limma_design.csv",row.names = 1) #236
subgroup <- factor(subgroup$TC_group)
design <- model.matrix(~0+subgroup)

## remove genes that consistently have 0 counts
dge <- DGEList(counts=limma.data)
keep <- filterByExpr(dge, design = NULL)
dge <- dge[keep,,keep.lib.sizes=FALSE] #18284   236
limma.data <- dge$counts

## remove batch effect
batch <- c(rep("1",sum(subgroup!=0)),rep("2",sum(subgroup==0)))
pca_res <- prcomp(t(limma.data), center = F,scale = F)
autoplot(pca_res, data =t(limma.data), colour = batch)
combat.data <- ComBat_seq(counts=limma.data, batch=batch)
limma.data.br<-as.data.frame(combat.data)
pca_res <- prcomp(t(limma.data.br), scale = F,center=F)
autoplot(pca_res, data = t(limma.data.br), colour = batch)

## normalization on batch effect removed data
dge <- DGEList(counts=limma.data.br) #18284   236
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3) #18284   236

## contrast fit limma model
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
cont.matrix <- makeContrasts(C1vsNormal=subgroup1-subgroup0,
                             C2vsNormal=subgroup2-subgroup0,
                             C1vsC2=subgroup1-subgroup2,
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## extract limma results
### 1vs0
res1<-topTable(fit2,coef=1, adjust="BH",number=nrow(limma.data))
res<-res1
#write.csv(res,"TC_limma_1vs0.csv")
res<-as.data.frame(res)
res$diffexpressed <- "NO"
res$diffexpressed[res$logFC > 0.5 & res$adj.P.Val < 0.05] <- "UP"
res$diffexpressed[res$logFC < -0.5 & res$adj.P.Val < 0.05] <- "DOWN"
res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- row.names(res)[res$diffexpressed != "NO"]
see=ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.5, alpha=0.5) +
  theme_minimal() +
  geom_text_repel(size = 2) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
print(see)
### 2vs0
res2<-topTable(fit2,coef=2, adjust="BH",number=nrow(limma.data))
res<-res2
#write.csv(res,"TC_limma_2vs0.csv")
res<-as.data.frame(res)
res$diffexpressed <- "NO"
res$diffexpressed[res$logFC > 0.5 & res$adj.P.Val < 0.05] <- "UP"
res$diffexpressed[res$logFC < -0.5 & res$adj.P.Val < 0.05] <- "DOWN"
res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- row.names(res)[res$diffexpressed != "NO"]
see=ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.5, alpha=0.5) +
  theme_minimal() +
  geom_text_repel(size = 2) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
print(see)
### 1vs2
res3<-topTable(fit2,coef=3, adjust="BH",number=nrow(limma.data))
res<-res3
#write.csv(res,"TC_limma_1vs2.csv")
res<-as.data.frame(res)
res$diffexpressed <- "NO"
res$diffexpressed[res$logFC > 1 & res$adj.P.Val < 0.05] <- "UP"
res$diffexpressed[res$logFC < -1 & res$adj.P.Val < 0.05] <- "DOWN"
res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- row.names(res)[res$diffexpressed != "NO"]
see=ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.5, alpha=0.5) +
  theme_minimal() +
  geom_text_repel(size = 2) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
print(see)

###############################################################################################################

#####################################      construct the gene signature     ###########################################

## Cox and limma overlapped genes
gs<-intersect(row.names(res)[abs(res$logFC) > 1.1 & res$adj.P.Val < 0.05],row.names(logCPM_cox)) #15 genes
### heatmap (TCGA gene signature)
TC<-rbind(t(subgroup),logCPM[gs,])[,1:123]
row.names(TC)[1]<-"group"
TC_sort<-TC[,order(TC[1,])]
col<-TC_sort[1,]
col[col==2]<-"#F8766D"
col[col==3]<-"#00BFC4"
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
hm <- heatmap.3(TC_sort[-1,],scale = "row",dendrogram = "none",
                key=T,keysize=2,hclustfun=myclust, distfun=mydist,
                ColSideColors=as.matrix(col),
                ColSideColorsSize=1,Colv = F,Rowv = T,
                col = colorRampPalette(c("black","darkblue","blue3","blue2","blue1",
                                         "orange","red","red1","red2","red3","darkred")))

###############################################################################################################

#####################################      validation of the gene signature     ###########################################

## Merge the gene expression data of TCGA, METABRIC, and GSE149283, then do batch effect remove.
### read GSE149283 gene signature data and drug response data
GSE<-read.csv("GSE149283_gxp_pp.csv",row.names = 1)
GSE_dg<-read.csv("GSE149283_drug.csv",row.names = 1)
### prepare METABRIC gene signature data
MB<-pp_MB_gxp[gs,]
### merge three gene expression datasets
gs.gxp<-as.data.frame(t(cbind(TC[-1,],MB,GSE)))
### imputation
gs.gxp<-data.imputation(gs.gxp,fun="median")
### batch effect remove
pca_res <- prcomp(gs.gxp, scale. = TRUE)
ntc<-ncol(TC[-1,])
nmb<-ncol(MB)
ngse<-ncol(GSE)
gs.batch<-c(rep("1",ntc),rep("2",nmb),rep("3",ngse))
autoplot(pca_res, data = gs.gxp, colour = gs.batch)
combat.data = ComBat(dat=t(gs.gxp), batch=gs.batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca_res <- prcomp(t(combat.data), scale. = TRUE)
autoplot(pca_res, data = t(combat.data), colour = gs.batch)
### batch effect removed data
mydata<-scale(combat.data, center = FALSE, scale = apply(combat.data, 2, sd, na.rm = TRUE)) #15 241

## prepare separate data from the batch effect removed data
tc<-as.data.frame(t(mydata[,1:ntc]))
mb<-as.data.frame(t(mydata[,(ntc+1):(ntc+nmb)]))
gse<-as.data.frame(t(mydata[,(ntc+nmb+1):(ntc+nmb+ngse)]))

## prepare TCGA group label data
tc.label<-as.factor(as.character(ifelse(TC[1,]==2,"Subgroup1","Subgroup2"))) #Subgroup1: "0"; Subgroup2: "1"

## normalization on TCGA and METABRIC
tc_mean<-colMeans(tc)
tc_sd<-apply(tc, 2, sd)
mb_mean<-colMeans(mb)
mb_sd<-apply(mb, 2, sd)
tc_nor<-as.data.frame(t(apply(tc,1,function(X){(X - tc_mean)/tc_sd})))
mb_nor<-as.data.frame(t(apply(mb,1,function(X){(X - tc_mean)/tc_sd})))
gse_nor<-as.data.frame(t(apply(gse,1,function(X){(X - tc_mean)/tc_sd})))

## supervised classification
tc_nor$label<-tc.label
trc <- trainControl("cv",10,savePred=TRUE)
### xgboost trainning on TCGA data
set.seed(28)
xgb_grid = expand.grid(nrounds = 100,max_depth = c(1, 2, 4, 6, 8,10, 12),
  eta=c(1, 0.5, 0.1, 0.07, 0.05, 0.01, 0.005),gamma = 0.5,colsample_bytree=0.5,
  min_child_weight=1, subsample=0.5)
xgb_trcontrol = trainControl(method = "cv",number = 5,verboseIter = TRUE,
  returnData = FALSE,returnResamp = "all",classProbs = TRUE,                                                  
  summaryFunction = twoClassSummary,allowParallel = TRUE)
xgb_train = train(x = as.matrix(tc_nor[,1:15]),y = tc.label,
  trControl = xgb_trcontrol,tuneGrid = xgb_grid, method = "xgbTree")

library(xgboost)
xgb.train = xgb.DMatrix(data = as.matrix(tc_nor[,1:15]), label = tc.label)
model_xgboost = xgboost(data = xgb.train, max.depth = 2, nrounds = 100, verbose = 0)
importance_matrix = xgb.importance(colnames(tc_nor[,1:15]), model = model_xgboost)
xgb.plot.importance(importance_matrix)



### xgboost testing on MATABRIC data
mb.pred<-predict(xgb_train, newdata = mb_nor)
mb.surv<-cbind(mb.pred,mb.surv.data[row.names(mb_nor),])
ggsurvplot(fit, risk.table = TRUE, pval = TRUE, main="METABRIC",pval.coord = c(0, 0.55),
           xlab="Survival time (days)", ylab= "Survival probability",ylim = c(0, 1),
           legend.labs = c('Subgroup 1','Subgroup 2'))
MB_ht<-cbind(mb.pred,mb)
MB_ht<-MB_ht[order(MB_ht$mb.pred),]
ha = HeatmapAnnotation(bar = sort(MB_ht$mb.pred))
Heatmap(t(scale(MB_ht[,-1])), name = "mb",  top_annotation = ha, 
        cluster_rows = TRUE,cluster_columns = F,
        column_title = NULL)


### gene expression differential analysis on METABRIC data
limma.data<-read.csv("MB_gxp_pp.csv",row.names = 1) #24368 genes, 104 patients
limma.data<-data.imputation(limma.data,fun="median")
subgroup <- ifelse(mb.pred=='Subgroup1','1','2')
design <- model.matrix(~0+subgroup)

limma_nor<-min(logCPM) + ( ((limma.data - min(limma.data)) * (max(logCPM) - min(logCPM))) / (max(limma.data)-min(limma.data)))


## contrast fit limma model
fit <- lmFit(limma_nor, design)
fit <- eBayes(fit, trend=TRUE)
cont.matrix <- makeContrasts(C1vsC2=subgroup1-subgroup2,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## extract limma results
res1<-topTable(fit2,coef=1, adjust="BH",number=nrow(limma.data))
res<-res1
#write.csv(res,"MB_limma_1vs2.csv")
res<-as.data.frame(res)
res$diffexpressed <- "NO"
res$diffexpressed[res$logFC > 1 & res$P.Value < 0.05] <- "UP"
res$diffexpressed[res$logFC < -1 & res$P.Value < 0.05] <- "DOWN"
res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- row.names(res)[res$diffexpressed != "NO"]
see=ggplot(data=res, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.5, alpha=0.5) +
  theme_minimal() +
  geom_text_repel(size = 2) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
print(see)







### check clinical relevance of the predicted METABRIC subgroups
pp_MB_clinic<-pp_MB_clinic[row.names(mb_nor),]
### age
summary(pp_MB_clinic$age_at_diagnosis[mb.pred=='Subgroup1'])
sd(pp_MB_clinic$age_at_diagnosis[mb.pred=='Subgroup1'])
summary(pp_MB_clinic$age_at_diagnosis[mb.pred=='Subgroup2'])
sd(pp_MB_clinic$age_at_diagnosis[mb.pred=='Subgroup2'])
t.test(pp_MB_clinic$age_at_diagnosis[mb.pred=='Subgroup1'],pp_MB_clinic$age_at_diagnosis[mb.pred=='Subgroup2'])
### size
summary(as.numeric(pp_MB_clinic$size[mb.pred=='Subgroup1']))
sd(as.numeric(pp_MB_clinic$size[mb.pred=='Subgroup1']),na.rm = T)
summary(as.numeric(pp_MB_clinic$size[mb.pred=='Subgroup2']))
sd(as.numeric(pp_MB_clinic$size[mb.pred=='Subgroup2']),na.rm = T)
t.test(as.numeric(pp_MB_clinic$size[mb.pred=='Subgroup1']),as.numeric(pp_MB_clinic$size[mb.pred=='Subgroup2']))
### lymph_nodes_positive
summary(pp_MB_clinic$lymph_nodes_positive[mb.pred=='Subgroup1'])
sd(pp_MB_clinic$lymph_nodes_positive[mb.pred=='Subgroup1'])
summary(pp_MB_clinic$lymph_nodes_positive[mb.pred=='Subgroup2'])
sd(pp_MB_clinic$lymph_nodes_positive[mb.pred=='Subgroup2'])
t.test(pp_MB_clinic$lymph_nodes_positive[mb.pred=='Subgroup1'],pp_MB_clinic$lymph_nodes_positive[mb.pred=='Subgroup2'])
### grade
table(pp_MB_clinic$grade[mb.pred=='Subgroup1'])
table(pp_MB_clinic$grade[mb.pred=='Subgroup2'])
chisq.test(c(2,18,35,2),c(1,14,32,0))
### stage
table(pp_MB_clinic$stage[mb.pred=='Subgroup1'])
table(pp_MB_clinic$stage[mb.pred=='Subgroup2'])
chisq.test(c(15,14,13,2,1,12),c(9,8,19,2,0,9))
### Treatment
table(pp_MB_clinic$Treatment[mb.pred=='Subgroup1'])
table(pp_MB_clinic$Treatment[mb.pred=='Subgroup2'])
chisq.test(c(2,0,5,3,12,23,8,4),c(1,3,6,2,9,21,3,2))
### PAM50 subtypes
table(pp_MB_clinic$Pam50Subtype[mb.pred=='Subgroup1'])
table(pp_MB_clinic$Pam50Subtype[mb.pred=='Subgroup2'])
chisq.test(c(1,16,15,21,4),c(1,16,6,21,3))


### xgboost testing on GSE data
gse.pred<-predict(xgb_train, newdata = gse_nor)
GSE_ht<-cbind(gse.pred,gse)
GSE_ht<-GSE_ht[order(GSE_ht$gse.pred),]
ha = HeatmapAnnotation(bar = sort(GSE_ht$gse.pred))
Heatmap(t(scale(GSE_ht[,-1])), name = "mb",  top_annotation = ha, 
        cluster_rows = TRUE,cluster_columns = F,
        column_title = NULL)
ggplot(as.data.frame(cbind(t(GSE_dg),gse.pred)) %>% count(gse.pred, drugresponse), aes(gse.pred, n, fill=drugresponse)) +
  geom_bar(stat="identity")

table(t(GSE_dg),gse.pred)
#####################################  compare subgroups in TCGA and METABRIC        ######################################################


# mutation

##check mutation of the TCGA subgroups
TC1_mut<-read.csv("TC1_Mutated_Genes.txt",sep = '\t',row.names = 1) #3109
TC2_mut<-read.csv("TC2_Mutated_Genes.txt",sep = '\t',row.names = 1) #4488
TC2_mut_ug<-row.names(TC2_mut)[!row.names(TC2_mut) %in% row.names(TC1_mut)] # 3293 unique genes that mutated in subgroup2 not in subgroup1
TC1_mut_ug<-row.names(TC1_mut)[!row.names(TC1_mut) %in% row.names(TC2_mut)] # 1914 unique genes that mutated in subgroup1 not in subgroup2
TC2_mut_ug<-as.data.frame(cbind(TC2_mut_ug,rep('TC',length(TC2_mut_ug)),rep('Subgroup2',length(TC2_mut_ug))))
colnames(TC2_mut_ug)<-c('gene','type','subgroup') 
TC1_mut_ug<-as.data.frame(cbind(TC1_mut_ug,rep('TC',length(TC1_mut_ug)),rep('Subgroup1',length(TC1_mut_ug))))
colnames(TC1_mut_ug)<-c('gene','type','subgroup')
TC_mut_ug<-rbind(TC2_mut_ug,TC1_mut_ug)# uniquely mutated genes in two subgroups of TCGA

MB1_mut<-read.csv("MB1_Mutated_Genes.txt",sep = '\t',row.names = 1) #108
MB2_mut<-read.csv("MB2_Mutated_Genes.txt",sep = '\t',row.names = 1) #87
MB2_mut_ug<-row.names(MB2_mut)[!row.names(MB2_mut) %in% row.names(MB1_mut)]# 18 unique genes that mutated in subgroup2 not in subgroup1
MB1_mut_ug<-row.names(MB1_mut)[!row.names(MB1_mut) %in% row.names(MB2_mut)]# 39 unique genes that mutated in subgroup1 not in subgroup2
MB2_mut_ug<-as.data.frame(cbind(MB2_mut_ug,rep('MB',length(MB2_mut_ug)),rep('Subgroup2',length(MB2_mut_ug))))
colnames(MB2_mut_ug)<-c('gene','type','subgroup')
MB1_mut_ug<-as.data.frame(cbind(MB1_mut_ug,rep('MB',length(MB1_mut_ug)),rep('Subgroup1',length(MB1_mut_ug))))
colnames(MB1_mut_ug)<-c('gene','type','subgroup')
MB_mut_ug<-rbind(MB2_mut_ug,MB1_mut_ug)# uniquely mutated genes in two subgroups of METABRIC

mut_ug<-rbind(TC_mut_ug, MB_mut_ug)# uniquely mutated genes in two subgroups for both TCGA and METABRIC
# plot the number of mutated genes in different subgroups for different cohorts
barplot(table(mut_ug$type, mut_ug$subgroup) ,beside=T,legend=c('METABRIC','TCGA-BRCA')) 

# genes mutated in both TCGA subgroup 2 and METABRIC subgroup2 (for OncoPrint)
intersect((mut_ug %>% filter(subgroup=="Subgroup2") %>% filter(type=="TC"))[,1], 
          (mut_ug %>% filter(subgroup=="Subgroup2") %>% filter(type=="MB"))[,1]) 

intersect((mut_ug %>% filter(subgroup=="Subgroup1") %>% filter(type=="TC"))[,1], 
          (mut_ug %>% filter(subgroup=="Subgroup1") %>% filter(type=="MB"))[,1]) 



### check public gene signature relevance of the TCGA subgroups and predicted METABRIC subgroups
data(sig.endoPredict)
data(sig.oncotypedx)
data(sig.tamr13)
data(sig.gene70)
data(sig.pik3cags)
data(pam50)
data(sig.ggi)
data(sig.genius)
data(scmod1.robust)
data(pam50.robust)

library(genefu)

anno<-read.csv("anno.csv",header = F)
annot<-cbind(unlist(strsplit(anno$V1,"|",fixed = TRUE))[seq(1,36325,by=2)],
             unlist(strsplit(anno$V1,"|",fixed = TRUE))[seq(2,36326,by=2)])
colnames(annot)<-c("prob","EntrezGene.ID")
annot<-as.data.frame(annot)
annot[annot$prob=='SLC35E2',]
#         prob EntrezGene.ID
#14519 SLC35E2        728661
#14520 SLC35E2          9906
annot[14520,'prob']<-"SLC35E2.2"

TC.gxp<-read.csv("TC_gxp_pp.csv",row.names=1) #55368 genes, 123 patients 
dge <- DGEList(counts=TC.gxp)
dge <- calcNormFactors(dge)
TC.gxp <- cpm(dge, log=TRUE, prior.count=3) #55368   123
MB.gxp<-read.csv("MB_gxp_pp.csv",row.names = 1) #24368 genes, 104 patients

TC.gxp<-TC.gxp[intersect(row.names(TC.gxp),annot$prob),]
TC.annot<-annot[annot$prob  %in% intersect(row.names(TC.gxp),annot$prob),]
row.names(TC.annot)<-TC.annot$prob

MB.gxp<-MB.gxp[intersect(row.names(MB.gxp),annot$prob),]
MB.annot<-annot[annot$prob  %in% intersect(row.names(MB.gxp),annot$prob),]
row.names(MB.annot)<-MB.annot$prob

#TCGA
pik3cags.tc <- pik3cags(data=t(TC.gxp), annot=TC.annot, do.mapping=TRUE)
rorS.tc<-rorS(data=t(TC.gxp), annot= TC.annot, do.mapping = F)$score
ggi.tc<-ggi(data=t(TC.gxp), annot= TC.annot, do.mapping = TRUE)$score
gene70.tc<-gene70(data=t(TC.gxp), annot= TC.annot, do.mapping = TRUE,std ="scale")$score
genius.tc<-genius(data=t(TC.gxp), annot= TC.annot, do.mapping = TRUE, do.scale = F)$score


tc.gs<-as.data.frame(cbind(TC_group,pik3cags.tc,rorS.tc,ggi.tc,gene70.tc,genius.tc))


pam50.tc<-molecular.subtyping(sbt.model = "pam50",data=t(TC.gxp),
                    annot=TC.annot,do.mapping=F)$subtype

pam50.tc<-factor(pam50.tc, levels=c('Basal', 'Her2', 'LumA', 'LumB', 'Normal'))

table(pam50.tc[TC_group==1])
table(pam50.tc[TC_group==2])

chisq.test( c(5, 15, 13, 13,17),c( 5, 14 , 21,18, 2) )

barplot(table(pam50.tc, TC_group) ,beside=T) 


ggplot(tc.gs, aes(x=pik3cags.tc, fill=as.factor(TC_group))) +
  geom_density(alpha=0.4)
t.test(tc.gs$pik3cags.tc[tc.gs$TC_group==1],tc.gs$pik3cags.tc[tc.gs$TC_group==2])

ggplot(tc.gs, aes(x=rorS.tc, fill=as.factor(TC_group))) +
  geom_density(alpha=0.4)
t.test(tc.gs$rorS.tc[tc.gs$TC_group==1],tc.gs$rorS.tc[tc.gs$TC_group==2])

ggplot(tc.gs, aes(x=ggi.tc, fill=as.factor(TC_group))) +
  geom_density(alpha=0.4)
t.test(tc.gs$ggi.tc[tc.gs$TC_group==1],tc.gs$ggi.tc[tc.gs$TC_group==2])

ggplot(tc.gs, aes(x=gene70.tc, fill=as.factor(TC_group))) +
  geom_density(alpha=0.4)
t.test(tc.gs$gene70.tc[tc.gs$TC_group==1],tc.gs$gene70.tc[tc.gs$TC_group==2])

ggplot(tc.gs, aes(x=genius.tc, fill=as.factor(TC_group))) +
  geom_density(alpha=0.4)
t.test(tc.gs$genius.tc[tc.gs$TC_group==1],tc.gs$genius.tc[tc.gs$TC_group==2])


#METABRIC
pik3cags.mb <- pik3cags(data=t(MB.gxp), annot=MB.annot, do.mapping=TRUE)
rorS.mb<-rorS(data=t(MB.gxp), annot= MB.annot, do.mapping = F)$score
ggi.mb<-ggi(data=t(MB.gxp), annot= MB.annot, do.mapping = TRUE)$score
gene70.mb<-gene70(data=t(MB.gxp), annot= MB.annot, do.mapping = TRUE,std ="scale")$score
genius.mb<-genius(data=t(MB.gxp), annot= MB.annot, do.mapping = TRUE, do.scale = F)$score
mb.gs<-as.data.frame(cbind(mb.pred,pik3cags.mb,rorS.mb,ggi.mb,gene70.mb,genius.mb))


pam50.mb<-molecular.subtyping(sbt.model = "pam50",data=t(MB.gxp),
                              annot=MB.annot,do.mapping=F)$subtype

table(pam50.mb[mb.pred=="Subgroup1"])
table(pam50.mb[mb.pred=="Subgroup2"])

barplot(table(pam50.mb, mb.pred) ,beside=T,legend=c('Basal','Her2','LumB','LumA','Normal')) 
barplot(table(pp_MB_clinic$Pam50Subtype, mb.pred) ,beside=T) 


ggplot(mb.gs, aes(x=pik3cags.mb, fill=as.factor(mb.pred))) +
  geom_density(alpha=0.4)
t.test(mb.gs$pik3cags.mb[mb.gs$mb.pred==1],mb.gs$pik3cags.mb[mb.gs$mb.pred==2])

ggplot(mb.gs, aes(x=rorS.mb, fill=as.factor(mb.pred))) +
  geom_density(alpha=0.4)
t.test(mb.gs$rorS.mb[mb.gs$mb.pred==1],mb.gs$rorS.mb[mb.gs$mb.pred==2])

ggplot(mb.gs, aes(x=ggi.mb, fill=as.factor(mb.pred))) +
  geom_density(alpha=0.4)
t.test(mb.gs$ggi.mb[mb.gs$mb.pred==1],mb.gs$ggi.mb[mb.gs$mb.pred==2])

ggplot(mb.gs, aes(x=gene70.mb, fill=as.factor(mb.pred))) +
  geom_density(alpha=0.4)
t.test(mb.gs$gene70.mb[mb.gs$mb.pred==1],mb.gs$gene70.mb[mb.gs$mb.pred==2])

ggplot(mb.gs, aes(x=genius.mb, fill=as.factor(ifelse(mb.pred==1,2,1)))) +
  geom_density(alpha=0.4)
t.test(mb.gs$genius.mb[mb.gs$mb.pred==1],mb.gs$genius.mb[mb.gs$mb.pred==2])


### check TIL relevance of the TCGA subgroups and the predicted METABRIC subgroups
pp_TC_timer<-read.csv("pp_TC_timer.csv",row.names = 1)
pp_TC_timer$label<-ifelse(TC_group==1,"Subgroup1","Subgroup2")

ggplot(pp_TC_timer, aes(x=B.cell, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup1"],
       pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup2"])

ggplot(pp_TC_timer, aes(x=T.cell.CD4., fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_TC_timer$T.cell.CD4.[pp_TC_timer$label=="Subgroup1"],
       pp_TC_timer$T.cell.CD4.[pp_TC_timer$label=="Subgroup2"])

ggplot(pp_TC_timer, aes(x=T.cell.CD8., fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_TC_timer$T.cell.CD8.[pp_TC_timer$label=="Subgroup1"],
       pp_TC_timer$T.cell.CD8.[pp_TC_timer$label=="Subgroup2"])

ggplot(pp_TC_timer, aes(x=Neutrophil, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup1"],
       pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup2"])

ggplot(pp_TC_timer, aes(x=Macrophage, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup1"],
       pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup2"])

ggplot(pp_TC_timer, aes(x=Dendritic.cell, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup1"],
       pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup2"])


### check TIL relevance of the predicted METABRIC subgroups
pp_MB_timer<-read.csv("pp_MB_timer.csv",row.names = 1)
pp_MB_timer$label<-mb.pred


ggplot(pp_MB_timer, aes(x=B.cell, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup2"],
       pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup1"])

ggplot(pp_MB_timer, aes(x=T.cell.CD4., fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_MB_timer$T.cell.CD4.[pp_MB_timer$label=="Subgroup2"],
       pp_MB_timer$T.cell.CD4.[pp_MB_timer$label=="Subgroup1"])

ggplot(pp_MB_timer, aes(x=T.cell.CD8., fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_MB_timer$T.cell.CD8.[pp_MB_timer$label=="Subgroup2"],
       pp_MB_timer$T.cell.CD8.[pp_MB_timer$label=="Subgroup1"])

ggplot(pp_MB_timer, aes(x=Neutrophil, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup2"],
       pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup1"])

ggplot(pp_MB_timer, aes(x=Macrophage, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup2"],
       pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup1"])

ggplot(pp_MB_timer, aes(x=Dendritic.cell, fill=label)) +
  geom_density(alpha=0.4)
t.test(pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup1"],
       pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup2"])


par(mai=c(1,1.5,0.5,0.5))
boxplot(pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup1"],
        pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup2"],
        
        pp_TC_timer$T.cell.CD4.[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$T.cell.CD4.[pp_MB_timer$label=="Subgroup1"],
        pp_TC_timer$T.cell.CD4.[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$T.cell.CD4.[pp_MB_timer$label=="Subgroup2"],
        
        pp_TC_timer$T.cell.CD8.[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$T.cell.CD8.[pp_MB_timer$label=="Subgroup1"],
        pp_TC_timer$T.cell.CD8.[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$T.cell.CD8.[pp_MB_timer$label=="Subgroup2"],
        
        pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup1"],
        pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup2"],
        
        pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup1"],
        pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup2"],
        
        pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup1"],
        pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup2"],
        
        at = c(1:4,6:9,11:14,16:19,21:24,26:29),
        main = "Subgroup TIL",
        names = c("B.cell","B.cell","B.cell","B.cell", "T.cell.CD4","T.cell.CD4","T.cell.CD4","T.cell.CD4",
                  "T.cell.CD8","T.cell.CD8","T.cell.CD8","T.cell.CD8","Neutrophil","Neutrophil",
                  "Neutrophil","Neutrophil","Macrophage","Macrophage","Macrophage","Macrophage",
                  "Dendritic.cell","Dendritic.cell","Dendritic.cell","Dendritic.cell"), las = 1,
        border = c(1,2,1,2, 1,2,1,2, 1,2,1,2, 1,2,1,2, 1,2,1,2, 1,2,1,2),
        
        col=c("orange","orange","grey","grey", "orange","orange","grey","grey", "orange","orange","grey","grey",
              "orange","orange","grey","grey", "orange","orange","grey","grey", "orange","orange","grey","grey"),
        horizontal = TRUE)





par(mai=c(1,1.5,0.5,0.5))
boxplot(c(pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup1"]),
        c(pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup2"]),
        
        c(pp_TC_timer$T.cell.CD4.[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$T.cell.CD4.[pp_MB_timer$label=="Subgroup1"]),
        c(pp_TC_timer$T.cell.CD4.[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$T.cell.CD4.[pp_MB_timer$label=="Subgroup2"]),
        
        c(pp_TC_timer$T.cell.CD8.[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$T.cell.CD8.[pp_MB_timer$label=="Subgroup1"]),
        c(pp_TC_timer$T.cell.CD8.[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$T.cell.CD8.[pp_MB_timer$label=="Subgroup2"]),
        
        c(pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup1"]),
        c(pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup2"]),
        
        c(pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup1"]),
        c(pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup2"]),
        
        c(pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup1"],pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup1"]),
        c(pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup2"],pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup2"]),
        
        at = c(1:2,4:5,7:8,10:11,13:14,16:17),
        main = "Subgroup TIL",
        names = c("B.cell","B.cell","T.cell.CD4","T.cell.CD4",
                  "T.cell.CD8","T.cell.CD8","Neutrophil","Neutrophil",
                 "Macrophage","Macrophage","Dendritic.cell","Dendritic.cell"), las = 1,
        border = c(1,2,1,2, 1,2,1,2, 1,2,1,2),
        
        col=c("orange","grey", "orange","grey", "orange","grey",
              "orange","grey", "orange","grey", "orange","grey"),
        horizontal = TRUE)


t.test(c(pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup1"]),
       c(pp_TC_timer$B.cell[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$B.cell[pp_MB_timer$label=="Subgroup2"]))

t.test(c(pp_TC_timer$T.cell.CD4[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$T.cell.CD4[pp_MB_timer$label=="Subgroup1"]),
       c(pp_TC_timer$T.cell.CD4[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$T.cell.CD4[pp_MB_timer$label=="Subgroup2"]))

t.test(c(pp_TC_timer$T.cell.CD8[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$T.cell.CD8[pp_MB_timer$label=="Subgroup1"]),
       c(pp_TC_timer$T.cell.CD8[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$T.cell.CD8[pp_MB_timer$label=="Subgroup2"]))

t.test(c(pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup1"]),
       c(pp_TC_timer$Neutrophil[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$Neutrophil[pp_MB_timer$label=="Subgroup2"]))

t.test(c(pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup1"]),
       c(pp_TC_timer$Macrophage[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$Macrophage[pp_MB_timer$label=="Subgroup2"]))

t.test(c(pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup1"], pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup1"]),
       c(pp_TC_timer$Dendritic.cell[pp_TC_timer$label=="Subgroup2"], pp_MB_timer$Dendritic.cell[pp_MB_timer$label=="Subgroup2"]))




##########################  DeepDep  ############################################

deepdep_61genes<-read.csv("61genes_deepdep.csv",header = F)


TC_deepdep<-read.csv("TC_predicted_data_model_exp_paper_demo.txt",sep = "\t",row.names = 1)


TC_deepdep_ht<-cbind(TC_group,as.data.frame(t(TC_deepdep)))
TC_deepdep_ht$TC_group<-as.factor(TC_deepdep_ht$TC_group)
TC_deepdep_ht<-TC_deepdep_ht[order(TC_deepdep_ht$TC_group),]


# TC_sig<-rep(NA,1298)
# for (i in 2:1299){
#   TC_sig[i-1]<-t.test(TC_deepdep_ht[,i][1:63],TC_deepdep_ht[,i][64:123])$p.value
# }
# TC_deepdep_ht_filter<-cbind(TC_group, TC_deepdep_ht[,-1][,TC_sig<0.0199])
# 

ha = HeatmapAnnotation(bar = sort(TC_deepdep_ht$TC_group))
Heatmap(t(scale((TC_deepdep_ht[,-1]))), name = "tc",  top_annotation = ha, 
        cluster_rows = T,cluster_columns = F,
        column_title = NULL)

TC_deepdep_ht_filter<-cbind(TC_group, TC_deepdep_ht[,-1][,deepdep_61genes$V1])
TC_deepdep_ht_filter$TC_group<-as.factor(TC_deepdep_ht_filter$TC_group)
ha = HeatmapAnnotation(bar = sort(TC_deepdep_ht_filter$TC_group))
Heatmap(t(scale((TC_deepdep_ht_filter[,-1]))), name = "tc",  top_annotation = ha, 
        cluster_rows = F,cluster_columns = F,
        column_title = NULL)


Heatmap(t(scale((TC_deepdep_ht_filter[1:63,-1]))), name = "tc", 
        cluster_rows = F,cluster_columns = T,
        column_title = NULL)
Heatmap(t(scale((TC_deepdep_ht_filter[64:123,-1]))), name = "tc", 
        cluster_rows = F,cluster_columns = T,
        column_title = NULL)


MB_deepdep<-read.csv("MB_predicted_data_model_exp_paper_demo.txt",sep = "\t",row.names = 1)
MB_deepdep_ht<-cbind(mb.pred,as.data.frame(t(MB_deepdep)))
MB_deepdep_ht<-MB_deepdep_ht[order(MB_deepdep_ht$mb.pred),]

# MB_sig<-rep(NA,1298)
# for (i in 2:1299){
#   MB_sig[i-1]<-t.test(MB_deepdep_ht[,i][1:57],MB_deepdep_ht[,i][58:104])$p.value
# }
# MB_deepdep_ht_filter<-cbind(mb.pred, MB_deepdep_ht[,-1][,MB_sig<0.467])


ha = HeatmapAnnotation(bar = sort(MB_deepdep_ht$mb.pred))
Heatmap(t(scale((MB_deepdep_ht[,-1]))), name = "mb",  top_annotation = ha, 
        cluster_rows = T,cluster_columns = F,
        column_title = NULL)

MB_deepdep_ht_filter<-cbind(mb.pred, MB_deepdep_ht[,-1][,deepdep_61genes$V1])

ha = HeatmapAnnotation(bar = sort(MB_deepdep_ht_filter$mb.pred))
Heatmap(t(scale((MB_deepdep_ht_filter[,-1]))), name = "mb",  top_annotation = ha, 
        cluster_rows = F,cluster_columns = F,
        column_title = NULL)

Heatmap(t(scale((MB_deepdep_ht_filter[1:57,-1]))), name = "mb", 
        cluster_rows = F,cluster_columns = T,
        column_title = NULL)
Heatmap(t(scale((MB_deepdep_ht_filter[58:104,-1]))), name = "mb", 
        cluster_rows = F,cluster_columns = T,
        column_title = NULL)


