RNA_COAD <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/RNA_COAD.txt", row.names=1)
miRNA_COAD <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/miRNA_COAD.txt", row.names=1)
cna_COAD <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/cna_COAD.txt", row.names=1)
methyl_COAD <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/methyl_COAD.txt", row.names=1)

##We have separate file that converts column names into required format of cancerin. But we removed columns with "Normal" data as Lasso models only consider "Tumor" samples
sampleID_Edited <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/sampleID_Edited.txt")
sampleID_Edited<-sampleID_Edited[-grep("Normal",sampleID_Edited$Required_format),]
dim(sampleID_Edited)
#1091 4 ; earlier it was 1204 with 113 normal samples
library(textshape)
# we should tally columnnames with sampleID
t_RNA_COAD<-as.data.frame(t(RNA_COAD))
View(t(RNA_COAD))
t_RNA_COAD<-cbind(Row_names=rownames(t_RNA_COAD),t_RNA_COAD)
dim(RNA_COAD)
dim(t_RNA_COAD)

library(dplyr)
###left-join
New_RNA_COAD<-left_join(sampleID_Edited,t_RNA_COAD,by=c("mi_RNA_format"="Row_names"))
New_RNA_COAD<-t(New_RNA_COAD)
New_RNA_COAD2<-New_RNA_COAD[-c(1:3),]
colnames(New_RNA_COAD2)<-New_RNA_COAD2[1,]
New_RNA_COAD2<-New_RNA_COAD2[-1,]


t_miRNA_COAD<-as.data.frame(t(miRNA_COAD))
t_miRNA_COAD<-cbind(Row_names=rownames(t_miRNA_COAD),t_miRNA_COAD)

New_miRNA_COAD<-left_join(sampleID_Edited,t_miRNA_COAD,by=c("mi_RNA_format"="Row_names"))
New_miRNA_COAD<-t(New_miRNA_COAD)
New_miRNA_COAD2<-New_miRNA_COAD[-c(1:3),]
colnames(New_miRNA_COAD2)<-New_miRNA_COAD2[1,]
New_miRNA_COAD2<-New_miRNA_COAD2[-1,]

t_cna_COAD<-as.data.frame(t(cna_COAD))
t_cna_COAD<-cbind(Row_names=rownames(t_cna_COAD),t_cna_COAD)
t_cna_COAD<-as.matrix(t_cna_COAD)
New_cna_COAD<-left_join(sampleID_Edited,as.data.frame(t_cna_COAD),by=c("methyl_cna_format"="Row_names"))
New_cna_COAD<-t(New_cna_COAD)
New_cna_COAD2<-New_cna_COAD[-c(1:3),]
colnames(New_cna_COAD2)<-New_cna_COAD2[1,]
New_cna_COAD2<-New_cna_COAD2[-1,]

t_methyl_COAD<-as.data.frame(t(methyl_COAD))
t_methyl_COAD<-cbind(Row_names=rownames(t_methyl_COAD),t_methyl_COAD)
t_methyl_COAD<-as.matrix(t_methyl_COAD)
New_methyl_COAD<-left_join(sampleID_Edited,as.data.frame(t_methyl_COAD),by=c("methyl_cna_format"="Row_names"))
New_methyl_COAD<-t(New_methyl_COAD)
New_methyl_COAD2<-New_methyl_COAD[-c(1:3),]
colnames(New_methyl_COAD2)<-New_methyl_COAD2[1,]
New_methyl_COAD2<-New_methyl_COAD2[-1,]


Common_colnames<-Reduce(intersect,list(colnames(New_RNA_COAD2),colnames(New_miRNA_COAD2),colnames(New_cna_COAD2),colnames(New_methyl_COAD2)))
length(Common_colnames)
#537
View(Common_colnames)

###We have not convert geneID into symbol yet for RNA dataset
library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(New_RNA_COAD2)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
New_RNA_COAD2<-subset(New_RNA_COAD2,rownames(New_RNA_COAD2) %in% geneID_RNA$GENEID)
rownames(New_RNA_COAD2)<-geneID_RNA$SYMBOL
New_RNA_COAD2<-New_RNA_COAD2[duplicated(rownames(New_RNA_COAD2))==FALSE,]
dim(New_RNA_COAD2) ## this results large number of genes. therefore we filtered ONLY genes significant from GDCRNATools analysis in COAD data

##We filter ONLY differentially expressed genes (dePC) and lncRNAs (deLNC) from our previous study :)
dePC_COAD <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/dePC.txt", row.names=1)
deLNC_COAD <- read.delim("Z:/GENETICS/DNA_Methylation/COAD/deLNC.txt", row.names=1)

New_RNA_COAD2<-subset(New_RNA_COAD2,rownames(New_RNA_COAD2) %in% c(dePC_COAD$symbol,deLNC_COAD$symbol))
dim(New_RNA_COAD2)

###We observed that there is a large number of data lines in RNA, and target interaction files. Therefore we removed some genes that do not report for all 5 file types (except miRNA expression)

#We concerned on gene-lncRNA pairs that are highly correlated from MC_HC_ceOutput_COAD.Then we converted their ensemble ID to symbol

#Then we identified common genes in both lncRNA-based and cancerin analysis's input
New_cna_COAD2<-as.data.frame(New_cna_COAD2)
New_cna_COAD2<-subset(New_cna_COAD2,rownames(New_cna_COAD2) %in% rownames(New_RNA_COAD2))
New_methyl_COAD2<-subset(New_methyl_COAD2,rownames(New_methyl_COAD2) %in% rownames(New_RNA_COAD2))
miRNA.target.interactions<-subset(miRNA_target,miRNA_target$target %in% c(dePC_COAD$symbol,deLNC_COAD$symbol))
tf.target.interactions<-subset(tf_target,tf_target$target %in% c(dePC_COAD$symbol,deLNC_COAD$symbol))

write.table(miRNA.target.interactions,file="Z:/GENETICS/DNA_Methylation/COAD/miRNA.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(tf.target.interactions,file="Z:/GENETICS/DNA_Methylation/COAD/tf.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(New_methyl_COAD2,file="Z:/GENETICS/DNA_Methylation/COAD/methyl.txt",sep = "\t",quote = FALSE)
write.table(New_cna_COAD2,file="Z:/GENETICS/DNA_Methylation/COAD/cna.txt",sep = "\t",quote = FALSE)
write.table(New_miRNA_COAD2,file="Z:/GENETICS/DNA_Methylation/COAD/miRNA.txt",sep = "\t",quote = FALSE)
write.table(New_RNA_COAD2,file="Z:/GENETICS/DNA_Methylation/COAD/RNA.txt",sep = "\t",quote = FALSE)

#we had an issue with duplicated rownames in RNA expression data

New_RNA_COAD2<-New_RNA_COAD2[duplicated(rownames(New_RNA_COAD2))==FALSE,]

## We had an issue when we install "WGCNA" R package in HPC. Therefore we calculate correlation and p-value in desktop for "pair.dt" dataset got from run_Cancerin script

