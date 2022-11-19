RNA_BRCA <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/RNA_BRCA.txt", row.names=1)
miRNA_BRCA <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/miRNA_BRCA.txt", row.names=1)
cna_BRCA <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/cna_BRCA.txt", row.names=1)
methyl_BRCA <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/methyl_BRCA.txt", row.names=1)

colnames(RNA_BRCA)<-substr(colnames(RNA_BRCA),1,nchar(colnames(RNA_BRCA))-3)
colnames(miRNA_BRCA)<-substr(colnames(miRNA_BRCA),1,nchar(colnames(miRNA_BRCA))-3)
Intersected<-Reduce(intersect,list(colnames(RNA_BRCA),colnames(miRNA_BRCA),colnames(cna_BRCA),colnames(methyl_BRCA)))

RNA_BRCA<-RNA_BRCA[,colnames(RNA_BRCA) %in% Intersected]
miRNA_BRCA<-miRNA_BRCA[1:158,colnames(miRNA_BRCA) %in% Intersected]
cna_BRCA<-cna_BRCA[,colnames(cna_BRCA) %in% Intersected]
methyl_BRCA<-methyl_BRCA[,colnames(methyl_BRCA) %in% Intersected]

library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(RNA_BRCA)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
RNA_BRCA<-subset(RNA_BRCA,rownames(RNA_BRCA) %in% geneID_RNA$GENEID)
rownames(RNA_BRCA)<-geneID_RNA$SYMBOL
###Now we add "lnRNA" prefix for gene names in RNA dataset as it requires in cancerin pipeline
write.table(RNA_BRCA,file = "BRCA/RNA_BRCA_filtered.txt",sep="\t",quote=FALSE,row.names = TRUE)
RNA_BRCA_withdeLNC<- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/RNA_BRCA_withdeLNC.txt", row.names=1)
RNA_BRCA_withdeLNC<-RNA_BRCA_withdeLNC[,colnames(RNA_BRCA_withdeLNC) %in% Intersected]


###Let's finalise input for pipeline
RNA<-RNA_BRCA_withdeLNC
miRNA<-miRNA_BRCA
cna<-cna_BRCA
methyl<-methyl_BRCA
##We have separate file that converts column names into required format of cancerin. But we removed columns with "Normal" data as Lasso models only consider "Tumor" samples
sampleID_Edited <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/sampleID_Edited.txt")
sampleID_Edited<-sampleID_Edited[-grep("Normal",sampleID_Edited$Required_format),]
dim(sampleID_Edited)
#1091 4 ; earlier it was 1204 with 113 normal samples
library(textshape)
# we should tally columnnames with sampleID
t_RNA_BRCA<-as.data.frame(t(RNA_BRCA))
View(t(RNA_BRCA))
t_RNA_BRCA<-cbind(Row_names=rownames(t_RNA_BRCA),t_RNA_BRCA)
dim(RNA_BRCA)
dim(t_RNA_BRCA)

library(dplyr)
###left-join
New_RNA_BRCA<-left_join(sampleID_Edited,t_RNA_BRCA,by=c("mi_RNA_format"="Row_names"))
New_RNA_BRCA<-t(New_RNA_BRCA)
New_RNA_BRCA2<-New_RNA_BRCA[-c(1:3),]
colnames(New_RNA_BRCA2)<-New_RNA_BRCA2[1,]
New_RNA_BRCA2<-New_RNA_BRCA2[-1,]


t_miRNA_BRCA<-as.data.frame(t(miRNA_BRCA))
t_miRNA_BRCA<-cbind(Row_names=rownames(t_miRNA_BRCA),t_miRNA_BRCA)

New_miRNA_BRCA<-left_join(sampleID_Edited,t_miRNA_BRCA,by=c("mi_RNA_format"="Row_names"))
New_miRNA_BRCA<-t(New_miRNA_BRCA)
New_miRNA_BRCA2<-New_miRNA_BRCA[-c(1:3),]
colnames(New_miRNA_BRCA2)<-New_miRNA_BRCA2[1,]
New_miRNA_BRCA2<-New_miRNA_BRCA2[-1,]

t_cna_BRCA<-as.data.frame(t(cna_BRCA))
t_cna_BRCA<-cbind(Row_names=rownames(t_cna_BRCA),t_cna_BRCA)
t_cna_BRCA<-as.matrix(t_cna_BRCA)
New_cna_BRCA<-left_join(sampleID_Edited,as.data.frame(t_cna_BRCA),by=c("methyl_cna_format"="Row_names"))
New_cna_BRCA<-t(New_cna_BRCA)
New_cna_BRCA2<-New_cna_BRCA[-c(1:3),]
colnames(New_cna_BRCA2)<-New_cna_BRCA2[1,]
New_cna_BRCA2<-New_cna_BRCA2[-1,]

t_methyl_BRCA<-as.data.frame(t(methyl_BRCA))
t_methyl_BRCA<-cbind(Row_names=rownames(t_methyl_BRCA),t_methyl_BRCA)
t_methyl_BRCA<-as.matrix(t_methyl_BRCA)
New_methyl_BRCA<-left_join(sampleID_Edited,as.data.frame(t_methyl_BRCA),by=c("methyl_cna_format"="Row_names"))
New_methyl_BRCA<-t(New_methyl_BRCA)
New_methyl_BRCA2<-New_methyl_BRCA[-c(1:3),]
colnames(New_methyl_BRCA2)<-New_methyl_BRCA2[1,]
New_methyl_BRCA2<-New_methyl_BRCA2[-1,]


Common_colnames<-Reduce(intersect,list(colnames(New_RNA_BRCA2),colnames(New_miRNA_BRCA2),colnames(New_cna_BRCA2),colnames(New_methyl_BRCA2)))
length(Common_colnames)
#537
View(Common_colnames)

###We have not convert geneID into symbol yet for RNA dataset
library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(New_RNA_BRCA2)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
New_RNA_BRCA2<-subset(New_RNA_BRCA2,rownames(New_RNA_BRCA2) %in% geneID_RNA$GENEID)
rownames(New_RNA_BRCA2)<-geneID_RNA$SYMBOL
New_RNA_BRCA2<-New_RNA_BRCA2[duplicated(rownames(New_RNA_BRCA2))==FALSE,]
dim(New_RNA_BRCA2) ## this results large number of genes. therefore we filtered ONLY genes significant from GDCRNATools analysis in BRCA data

##We filter ONLY differentially expressed genes (dePC) and lncRNAs (deLNC) from our previous study :)
dePC_BRCA <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/dePC.txt", row.names=1)
deLNC_BRCA <- read.delim("Z:/GENETICS/DNA_Methylation/BRCA/deLNC.txt", row.names=1)

New_RNA_BRCA2<-subset(New_RNA_BRCA2,rownames(New_RNA_BRCA2) %in% c(dePC_BRCA$symbol,deLNC_BRCA$symbol))
dim(New_RNA_BRCA2)

###We observed that there is a large number of data lines in RNA, and target interaction files. Therefore we removed some genes that do not report for all 5 file types (except miRNA expression)

#We concerned on gene-lncRNA pairs that are highly correlated from MC_HC_ceOutput_BRCA.Then we converted their ensemble ID to symbol

#Then we identified common genes in both lncRNA-based and cancerin analysis's input
New_cna_BRCA2<-as.data.frame(New_cna_BRCA2)
New_cna_BRCA2<-subset(New_cna_BRCA2,rownames(New_cna_BRCA2) %in% rownames(New_RNA_BRCA2))
New_methyl_BRCA2<-subset(New_methyl_BRCA2,rownames(New_methyl_BRCA2) %in% rownames(New_RNA_BRCA2))
miRNA.target.interactions<-subset(miRNA_target,miRNA_target$target %in% c(dePC_BRCA$symbol,deLNC_BRCA$symbol))
tf.target.interactions<-subset(tf_target,tf_target$target %in% c(dePC_BRCA$symbol,deLNC_BRCA$symbol))

write.table(miRNA.target.interactions,file="Z:/GENETICS/DNA_Methylation/BRCA/miRNA.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(tf.target.interactions,file="Z:/GENETICS/DNA_Methylation/BRCA/tf.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(New_methyl_BRCA2,file="Z:/GENETICS/DNA_Methylation/BRCA/methyl.txt",sep = "\t",quote = FALSE)
write.table(New_cna_BRCA2,file="Z:/GENETICS/DNA_Methylation/BRCA/cna.txt",sep = "\t",quote = FALSE)
write.table(New_miRNA_BRCA2,file="Z:/GENETICS/DNA_Methylation/BRCA/miRNA.txt",sep = "\t",quote = FALSE)
write.table(New_RNA_BRCA2,file="Z:/GENETICS/DNA_Methylation/BRCA/RNA.txt",sep = "\t",quote = FALSE)

#we had an issue with duplicated rownames in RNA expression data

New_RNA_BRCA2<-New_RNA_BRCA2[duplicated(rownames(New_RNA_BRCA2))==FALSE,]

## We had an issue when we install "WGCNA" R package in HPC. Therefore we calculate correlation and p-value in desktop for "pair.dt" dataset got from run_Cancerin script

