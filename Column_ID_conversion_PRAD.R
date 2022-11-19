RNA_PRAD <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/RNA_PRAD.txt", row.names=1)
miRNA_PRAD <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/miRNA_PRAD.txt", row.names=1)
cna_PRAD <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/cna_PRAD.txt", row.names=1)
methyl_PRAD <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/methyl_PRAD.txt", row.names=1)

names(RNA_PRAD)<-substr(names(RNA_PRAD),1,nchar(names(RNA_PRAD))-3)
names(miRNA_PRAD)<-substr(names(miRNA_PRAD),1,nchar(names(miRNA_PRAD))-3)
Intersected<-Reduce(intersect,list(names(RNA_PRAD),names(miRNA_PRAD),names(cna_PRAD),names(methyl_PRAD)))

RNA_PRAD<-RNA_PRAD[,colnames(RNA_PRAD) %in% Intersected]
miRNA_PRAD<-miRNA_PRAD[1:60,colnames(miRNA_PRAD) %in% Intersected]
cna_PRAD<-cna_PRAD[,colnames(cna_PRAD) %in% Intersected]
methyl_PRAD<-methyl_PRAD[,colnames(methyl_PRAD) %in% Intersected]

library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(RNA_PRAD)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
RNA_PRAD<-subset(RNA_PRAD,rownames(RNA_PRAD) %in% geneID_RNA$GENEID)
rownames(RNA_PRAD)<-geneID_RNA$SYMBOL

##We have separate file that converts column names into required format of cancerin. But we removed columns with "Normal" data as Lasso models only consider "Tumor" samples
sampleID_Edited <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/sampleID_Edited.txt")
sampleID_Edited<-sampleID_Edited[-grep("Normal",sampleID_Edited$Required_format),]
dim(sampleID_Edited)

library(textshape)
View(t(RNA_PRAD))
# we should tally columnnames with sampleID
t_RNA_PRAD<-as.data.frame(t(RNA_PRAD))
t_RNA_PRAD<-cbind(Row_names=rownames(t_RNA_PRAD),t_RNA_PRAD)
dim(RNA_PRAD)
dim(t_RNA_PRAD)

###left-join
New_RNA_PRAD<-left_join(sampleID_Edited,t_RNA_PRAD,by=c("mi_RNA_format"="Row_names"))
New_RNA_PRAD<-t(New_RNA_PRAD)
New_RNA_PRAD2<-New_RNA_PRAD[-c(1:3),]
colnames(New_RNA_PRAD2)<-New_RNA_PRAD2[1,]
New_RNA_PRAD2<-New_RNA_PRAD2[-1,]


t_miRNA_PRAD<-as.data.frame(t(miRNA_PRAD))
t_miRNA_PRAD<-cbind(Row_names=rownames(t_miRNA_PRAD),t_miRNA_PRAD)

New_miRNA_PRAD<-left_join(sampleID_Edited,t_miRNA_PRAD,by=c("mi_RNA_format"="Row_names"))
New_miRNA_PRAD<-t(New_miRNA_PRAD)
New_miRNA_PRAD2<-New_miRNA_PRAD[-c(1:3),]
colnames(New_miRNA_PRAD2)<-New_miRNA_PRAD2[1,]
New_miRNA_PRAD2<-New_miRNA_PRAD2[-1,]

t_cna_PRAD<-as.data.frame(t(cna_PRAD))
t_cna_PRAD<-cbind(Row_names=rownames(t_cna_PRAD),t_cna_PRAD)
t_cna_PRAD<-as.matrix(t_cna_PRAD)
New_cna_PRAD<-left_join(sampleID_Edited,as.data.frame(t_cna_PRAD),by=c("methyl_cna_format"="Row_names"))
New_cna_PRAD<-t(New_cna_PRAD)
New_cna_PRAD2<-New_cna_PRAD[-c(1:3),]
colnames(New_cna_PRAD2)<-New_cna_PRAD2[1,]
New_cna_PRAD2<-New_cna_PRAD2[-1,]

t_methyl_PRAD<-as.data.frame(t(methyl_PRAD))
t_methyl_PRAD<-cbind(Row_names=rownames(t_methyl_PRAD),t_methyl_PRAD)
t_methyl_PRAD<-as.matrix(t_methyl_PRAD)
New_methyl_PRAD<-left_join(sampleID_Edited,as.data.frame(t_methyl_PRAD),by=c("methyl_cna_format"="Row_names"))
New_methyl_PRAD<-t(New_methyl_PRAD)
New_methyl_PRAD2<-New_methyl_PRAD[-c(1:3),]
colnames(New_methyl_PRAD2)<-New_methyl_PRAD2[1,]
New_methyl_PRAD2<-New_methyl_PRAD2[-1,]


###We have not convert geneID into symbol yet for RNA dataset
library(EnsDb.Hsapiens.v79)
Genes_RNA<-rownames(New_RNA_PRAD2)
geneID_RNA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes_RNA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID_RNA)
New_RNA_PRAD2<-subset(New_RNA_PRAD2,rownames(New_RNA_PRAD2) %in% geneID_RNA$GENEID)
rownames(New_RNA_PRAD2)<-geneID_RNA$SYMBOL
New_RNA_PRAD2<-New_RNA_PRAD2[duplicated(rownames(New_RNA_PRAD2))==FALSE,]
dim(New_RNA_PRAD2) ## this results large number of genes. therefore we filtered ONLY genes significant from GDCRNATools analysis in PRAD data

dePC_PRAD <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/dePC.txt", row.names=1)
deLNC_PRAD <- read.delim("Z:/GENETICS/DNA_Methylation/PRAD/deLNC.txt", row.names=1)
demiR_PRAD<-read.delim("Z:/GENETICS/DNA_Methylation/PRAD/demiR.txt", row.names=1)

New_RNA_PRAD2<-subset(New_RNA_PRAD2,rownames(New_RNA_PRAD2) %in% c(dePC_PRAD$symbol,deLNC_PRAD$symbol))
dim(New_RNA_PRAD2)

###We observed that there is a large number of data lines in RNA, and target interaction files. Therefore we removed some genes that do not report for all 5 file types (except miRNA expression)

#We concerned on gene-lncRNA pairs that are highly correlated from MC_HC_ceOutput_PRAD.Then we converted their ensemble ID to symbol

#Then we identified common genes in both lncRNA-based and cancerin analysis's input
New_cna_PRAD2<-as.data.frame(New_cna_PRAD2)
New_cna_PRAD2<-subset(New_cna_PRAD2,rownames(New_cna_PRAD2) %in% rownames(New_RNA_PRAD2))
New_methyl_PRAD2<-subset(New_methyl_PRAD2,rownames(New_methyl_PRAD2) %in% rownames(New_RNA_PRAD2))
miRNA.target.interactions<-subset(miRNA_target,miRNA_target$target %in% c(dePC_PRAD$symbol,deLNC_PRAD$symbol))
tf.target.interactions<-subset(tf_target,tf_target$target %in% c(dePC_PRAD$symbol,deLNC_PRAD$symbol))

write.table(miRNA.target.interactions,file="Z:/GENETICS/DNA_Methylation/PRAD/miRNA.target.interactions.txt",quote = FALSE,row.names = FALSE)
write.table(tf.target.interactions,file="Z:/GENETICS/DNA_Methylation/PRAD/tf.target.interactions.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(New_methyl_PRAD2,file="Z:/GENETICS/DNA_Methylation/PRAD/methyl.txt",sep = "\t",quote = FALSE)
write.table(New_cna_PRAD2,file="Z:/GENETICS/DNA_Methylation/PRAD/cna.txt",sep = "\t",quote = FALSE)
write.table(New_miRNA_PRAD2,file="Z:/GENETICS/DNA_Methylation/PRAD/miRNA.txt",sep = "\t",quote = FALSE)
write.table(New_RNA_PRAD2,file="Z:/GENETICS/DNA_Methylation/PRAD/RNA.txt",sep = "\t",quote = FALSE)

