RNA_READ <- read.delim("Z:/GENETICS/DNA_Methylation/READ/RNA_READ.txt", row.names=1)
miRNA_READ <- read.delim("Z:/GENETICS/DNA_Methylation/READ/miRNA_READ.txt", row.names=1)
cna_READ <- read.delim("Z:/GENETICS/DNA_Methylation/READ/cna_READ.txt", row.names=1)
methyl_READ <- read.delim("Z:/GENETICS/DNA_Methylation/READ/methyl_READ.txt", row.names=1)

##We have separate file that converts column names into required format of cancerin. But we removed columns with "Normal" data as Lasso models only consider "Tumor" samples
sampleID_Edited <- read.delim("Z:/GENETICS/DNA_Methylation/READ/sampleID_Edited.txt")
sampleID_Edited<-sampleID_Edited[-grep("Normal",sampleID_Edited$Required_format),]
dim(sampleID_Edited)
#1091 4 ; earlier it was 1204 with 113 normal samples
library(textshape)
# we should tally columnnames with sampleID
t_RNA_READ<-as.data.frame(t(RNA_READ))
View(t(RNA_READ))
t_RNA_READ<-cbind(Row_names=rownames(t_RNA_READ),t_RNA_READ)
dim(RNA_READ)
dim(t_RNA_READ)

library(dplyr)
###left-join
New_RNA_READ<-left_join(sampleID_Edited,t_RNA_READ,by=c("mi_RNA_format"="Row_names"))
New_RNA_READ<-t(New_RNA_READ)
New_RNA_READ2<-New_RNA_READ[-c(1:3),]
colnames(New_RNA_READ2)<-New_RNA_READ2[1,]
New_RNA_READ2<-New_RNA_READ2[-1,]


t_miRNA_READ<-as.data.frame(t(miRNA_READ))
t_miRNA_READ<-cbind(Row_names=rownames(t_miRNA_READ),t_miRNA_READ)

New_miRNA_READ<-left_join(sampleID_Edited,t_miRNA_READ,by=c("mi_RNA_format"="Row_names"))
New_miRNA_READ<-t(New_miRNA_READ)
New_miRNA_READ2<-New_miRNA_READ[-c(1:3),]
colnames(New_miRNA_READ2)<-New_miRNA_READ2[1,]
New_miRNA_READ2<-New_miRNA_READ2[-1,]

t_cna_READ<-as.data.frame(t(cna_READ))
t_cna_READ<-cbind(Row_names=rownames(t_cna_READ),t_cna_READ)
t_cna_READ<-as.matrix(t_cna_READ)
New_cna_READ<-left_join(sampleID_Edited,as.data.frame(t_cna_READ),by=c("methyl_cna_format"="Row_names"))
New_cna_READ<-t(New_cna_READ)
New_cna_READ2<-New_cna_READ[-c(1:3),]
colnames(New_cna_READ2)<-New_cna_READ2[1,]
New_cna_READ2<-New_cna_READ2[-1,]

t_methyl_READ<-as.data.frame(t(methyl_READ))
t_methyl_READ<-cbind(Row_names=rownames(t_methyl_READ),t_methyl_READ)
t_methyl_READ<-as.matrix(t_methyl_READ)
New_methyl_READ<-left_join(sampleID_Edited,as.data.frame(t_methyl_READ),by=c("methyl_cna_format"="Row_names"))
New_methyl_READ<-t(New_methyl_READ)
New_methyl_READ2<-New_methyl_READ[-c(1:3),]
colnames(New_methyl_READ2)<-New_methyl_READ2[1,]
New_methyl_READ2<-New_methyl_READ2[-1,]


Common_colnames<-Reduce(intersect,list(colnames(New_RNA_READ2),colnames(New_miRNA_READ2),colnames(New_cna_READ2),colnames(New_methyl_READ2)))
length(Common_colnames)
#537
View(Common_colnames)

###We have not convert geneID into symbol yet for RNA dataset
library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(New_RNA_READ2)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
New_RNA_READ2<-subset(New_RNA_READ2,rownames(New_RNA_READ2) %in% geneID_RNA$GENEID)
rownames(New_RNA_READ2)<-geneID_RNA$SYMBOL
New_RNA_READ2<-New_RNA_READ2[duplicated(rownames(New_RNA_READ2))==FALSE,]
dim(New_RNA_READ2) ## this results large number of genes. therefore we filtered ONLY genes significant from GDCRNATools analysis in READ data

##We filter ONLY differentially expressed genes (dePC) and lncRNAs (deLNC) from our previous study :)
dePC_READ <- read.delim("Z:/GENETICS/DNA_Methylation/READ/dePC.txt", row.names=1)
deLNC_READ <- read.delim("Z:/GENETICS/DNA_Methylation/READ/deLNC.txt", row.names=1)

New_RNA_READ2<-subset(New_RNA_READ2,rownames(New_RNA_READ2) %in% c(dePC_READ$symbol,deLNC_READ$symbol))
dim(New_RNA_READ2)

###We observed that there is a large number of data lines in RNA, and target interaction files. Therefore we removed some genes that do not report for all 5 file types (except miRNA expression)

#We concerned on gene-lncRNA pairs that are highly correlated from MC_HC_ceOutput_READ.Then we converted their ensemble ID to symbol

#Then we identified common genes in both lncRNA-based and cancerin analysis's input
New_cna_READ2<-as.data.frame(New_cna_READ2)
New_cna_READ2<-subset(New_cna_READ2,rownames(New_cna_READ2) %in% rownames(New_RNA_READ2))
New_methyl_READ2<-subset(New_methyl_READ2,rownames(New_methyl_READ2) %in% rownames(New_RNA_READ2))
miRNA.target.interactions<-subset(miRNA_target,miRNA_target$target %in% c(dePC_READ$symbol,deLNC_READ$symbol))
tf.target.interactions<-subset(tf_target,tf_target$target %in% c(dePC_READ$symbol,deLNC_READ$symbol))

write.table(miRNA.target.interactions,file="Z:/GENETICS/DNA_Methylation/READ/miRNA.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(tf.target.interactions,file="Z:/GENETICS/DNA_Methylation/READ/tf.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(New_methyl_READ2,file="Z:/GENETICS/DNA_Methylation/READ/methyl.txt",sep = "\t",quote = FALSE)
write.table(New_cna_READ2,file="Z:/GENETICS/DNA_Methylation/READ/cna.txt",sep = "\t",quote = FALSE)
write.table(New_miRNA_READ2,file="Z:/GENETICS/DNA_Methylation/READ/miRNA.txt",sep = "\t",quote = FALSE)
write.table(New_RNA_READ2,file="Z:/GENETICS/DNA_Methylation/READ/RNA.txt",sep = "\t",quote = FALSE)

#we had an issue with duplicated rownames in RNA expression data

New_RNA_READ2<-New_RNA_READ2[duplicated(rownames(New_RNA_READ2))==FALSE,]

## We had an issue when we install "WGCNA" R package in HPC. Therefore we calculate correlation and p-value in desktop for "pair.dt" dataset got from run_Cancerin script

