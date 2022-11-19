#RNA_COLCA <- read.delim("Z:/GENETICS/DNA_Methylation/COLCA/RNA_COLCA.txt", row.names=1)
miRNA_COLCA <- read.delim("Z:/GENETICS/DNA_Methylation/COLCA/miRNA_COLCA.txt", row.names=1)
cna_COLCA <- read.delim("Z:/GENETICS/DNA_Methylation/COLCA/cna_COADREAD.txt", row.names=1)
methyl_COLCA <- read.delim("Z:/GENETICS/DNA_Methylation/COLCA/methyl_COADREAD.txt", row.names=1)

colnames(RNA_COLCA)<-substr(colnames(RNA_COLCA),1,nchar(colnames(RNA_COLCA))-3)
colnames(miRNA_COLCA)<-substr(colnames(miRNA_COLCA),1,nchar(colnames(miRNA_COLCA))-3)
RNA_COLCA<-RNA
Intersected<-Reduce(intersect,list(colnames(RNA_COLCA),colnames(miRNA_COLCA),colnames(cna_COLCA),colnames(methyl_COLCA)))

RNA_COLCA<-RNA_COLCA[,colnames(RNA_COLCA) %in% Intersected]
miRNA_COLCA<-miRNA_COLCA[1:339,colnames(miRNA_COLCA) %in% Intersected]
cna_COLCA<-cna_COLCA[,colnames(cna_COLCA) %in% Intersected]
methyl_COLCA<-methyl_COLCA[,colnames(methyl_COLCA) %in% Intersected]

RNA<-RNA_COLCA
miRNA<-miRNA_COLCA
cna<-cna_COLCA
methyl<-methyl_COLCA
