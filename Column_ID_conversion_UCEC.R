RNA_UCEC <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/RNA_UCEC.txt", row.names=1)
miRNA_UCEC <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/miRNA_UCEC.txt", row.names=1)
cna_UCEC <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/cna_UCEC.txt", row.names=1)
methyl_UCEC <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/methyl_UCEC.txt", row.names=1)

colnames(RNA_UCEC)<-substr(colnames(RNA_UCEC),1,nchar(colnames(RNA_UCEC))-3)
colnames(miRNA_UCEC)<-substr(colnames(miRNA_UCEC),1,nchar(colnames(miRNA_UCEC))-3)
Intersected<-Reduce(intersect,list(colnames(RNA_UCEC),colnames(miRNA_UCEC),colnames(cna_UCEC),colnames(methyl_UCEC)))

RNA_UCEC<-RNA_UCEC[,colnames(RNA_UCEC) %in% Intersected]
miRNA_UCEC<-miRNA_UCEC[1:158,colnames(miRNA_UCEC) %in% Intersected]
cna_UCEC<-cna_UCEC[,colnames(cna_UCEC) %in% Intersected]
methyl_UCEC<-methyl_UCEC[,colnames(methyl_UCEC) %in% Intersected]

RNA<-RNA_UCEC
miRNA<-miRNA_UCEC
cna<-cna_UCEC
methyl<-methyl_UCEC


##We have separate file that converts column names into required format of cancerin. But we removed columns with "Normal" data as Lasso models only consider "Tumor" samples
sampleID_Edited <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/sampleID_Edited.txt")
sampleID_Edited<-sampleID_Edited[-grep("Normal",sampleID_Edited$Required_format),]
dim(sampleID_Edited)
#1091 4 ; earlier it was 1204 with 113 normal samples
library(textshape)
# we should tally columnnames with sampleID
t_RNA_UCEC<-as.data.frame(t(RNA_UCEC))
View(t(RNA_UCEC))
t_RNA_UCEC<-cbind(Row_names=rownames(t_RNA_UCEC),t_RNA_UCEC)
dim(RNA_UCEC)
dim(t_RNA_UCEC)

library(dplyr)
###left-join
New_RNA_UCEC<-left_join(sampleID_Edited,t_RNA_UCEC,by=c("mi_RNA_format"="Row_names"))
New_RNA_UCEC<-t(New_RNA_UCEC)
New_RNA_UCEC2<-New_RNA_UCEC[-c(1:3),]
colnames(New_RNA_UCEC2)<-New_RNA_UCEC2[1,]
New_RNA_UCEC2<-New_RNA_UCEC2[-1,]


t_miRNA_UCEC<-as.data.frame(t(miRNA_UCEC))
t_miRNA_UCEC<-cbind(Row_names=rownames(t_miRNA_UCEC),t_miRNA_UCEC)

New_miRNA_UCEC<-left_join(sampleID_Edited,t_miRNA_UCEC,by=c("mi_RNA_format"="Row_names"))
New_miRNA_UCEC<-t(New_miRNA_UCEC)
New_miRNA_UCEC2<-New_miRNA_UCEC[-c(1:3),]
colnames(New_miRNA_UCEC2)<-New_miRNA_UCEC2[1,]
New_miRNA_UCEC2<-New_miRNA_UCEC2[-1,]

t_cna_UCEC<-as.data.frame(t(cna_UCEC))
t_cna_UCEC<-cbind(Row_names=rownames(t_cna_UCEC),t_cna_UCEC)
t_cna_UCEC<-as.matrix(t_cna_UCEC)
New_cna_UCEC<-left_join(sampleID_Edited,as.data.frame(t_cna_UCEC),by=c("methyl_cna_format"="Row_names"))
New_cna_UCEC<-t(New_cna_UCEC)
New_cna_UCEC2<-New_cna_UCEC[-c(1:3),]
colnames(New_cna_UCEC2)<-New_cna_UCEC2[1,]
New_cna_UCEC2<-New_cna_UCEC2[-1,]

t_methyl_UCEC<-as.data.frame(t(methyl_UCEC))
t_methyl_UCEC<-cbind(Row_names=rownames(t_methyl_UCEC),t_methyl_UCEC)
t_methyl_UCEC<-as.matrix(t_methyl_UCEC)
New_methyl_UCEC<-left_join(sampleID_Edited,as.data.frame(t_methyl_UCEC),by=c("methyl_cna_format"="Row_names"))
New_methyl_UCEC<-t(New_methyl_UCEC)
New_methyl_UCEC2<-New_methyl_UCEC[-c(1:3),]
colnames(New_methyl_UCEC2)<-New_methyl_UCEC2[1,]
New_methyl_UCEC2<-New_methyl_UCEC2[-1,]


Common_colnames<-Reduce(intersect,list(colnames(New_RNA_UCEC2),colnames(New_miRNA_UCEC2),colnames(New_cna_UCEC2),colnames(New_methyl_UCEC2)))
length(Common_colnames)
#537
View(Common_colnames)

###We have not convert geneID into symbol yet for RNA dataset
library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(New_RNA_UCEC2)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
New_RNA_UCEC2<-subset(New_RNA_UCEC2,rownames(New_RNA_UCEC2) %in% geneID_RNA$GENEID)
rownames(New_RNA_UCEC2)<-geneID_RNA$SYMBOL
New_RNA_UCEC2<-New_RNA_UCEC2[duplicated(rownames(New_RNA_UCEC2))==FALSE,]
dim(New_RNA_UCEC2) ## this results large number of genes. therefore we filtered ONLY genes significant from GDCRNATools analysis in UCEC data

##We filter ONLY differentially expressed genes (dePC) and lncRNAs (deLNC) from our previous study :)
dePC_UCEC <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/dePC.txt", row.names=1)
deLNC_UCEC <- read.delim("Z:/GENETICS/DNA_Methylation/UCEC/deLNC.txt", row.names=1)

New_RNA_UCEC2<-subset(New_RNA_UCEC2,rownames(New_RNA_UCEC2) %in% c(dePC_UCEC$symbol,deLNC_UCEC$symbol))
dim(New_RNA_UCEC2)

###We observed that there is a large number of data lines in RNA, and target interaction files. Therefore we removed some genes that do not report for all 5 file types (except miRNA expression)

#We concerned on gene-lncRNA pairs that are highly correlated from MC_HC_ceOutput_UCEC.Then we converted their ensemble ID to symbol

#Then we identified common genes in both lncRNA-based and cancerin analysis's input
New_cna_UCEC2<-as.data.frame(New_cna_UCEC2)
New_cna_UCEC2<-subset(New_cna_UCEC2,rownames(New_cna_UCEC2) %in% rownames(New_RNA_UCEC2))
New_methyl_UCEC2<-subset(New_methyl_UCEC2,rownames(New_methyl_UCEC2) %in% rownames(New_RNA_UCEC2))
miRNA.target.interactions<-subset(miRNA_target,miRNA_target$target %in% c(dePC_UCEC$symbol,deLNC_UCEC$symbol))
tf.target.interactions<-subset(tf_target,tf_target$target %in% c(dePC_UCEC$symbol,deLNC_UCEC$symbol))

write.table(miRNA.target.interactions,file="Z:/GENETICS/DNA_Methylation/UCEC/miRNA.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(tf.target.interactions_PRAD,file="Z:/GENETICS/DNA_Methylation/PRAD/tf.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(New_methyl_UCEC2,file="Z:/GENETICS/DNA_Methylation/UCEC/methyl.txt",sep = "\t",quote = FALSE)
write.table(New_cna_UCEC2,file="Z:/GENETICS/DNA_Methylation/UCEC/cna.txt",sep = "\t",quote = FALSE)
write.table(New_miRNA_UCEC2,file="Z:/GENETICS/DNA_Methylation/UCEC/miRNA.txt",sep = "\t",quote = FALSE)
write.table(New_RNA_UCEC2,file="Z:/GENETICS/DNA_Methylation/UCEC/RNA.txt",sep = "\t",quote = FALSE)

#we had an issue with duplicated rownames in RNA expression data

New_RNA_UCEC2<-New_RNA_UCEC2[duplicated(rownames(New_RNA_UCEC2))==FALSE,]

## We had an issue when we install "WGCNA" R package in HPC. Therefore we calculate correlation and p-value in desktop for "pair.dt" dataset got from run_Cancerin script

