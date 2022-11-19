#Here we used "Dorothea" R package to identify TF-gene interaction
library(dorothea)
View(dorothea_hs) #A table reporting signed human TF-target interactions. This database covers in total 1395 TFs targeting 20,244 
#genes with 486,676 unique interactions. In addition, each TF is accompanied with
#an emperical confidence level that was derived from the number of supporting evidences for this
#TF/interaction. The range is from A (high quality) to E (low quality).

dim(dorothea_hs)
#454505      4

#Here we require six data files:
  #1.rna expression data 2. miRNA expression data 3. methylation data 4. copy number alteration data 5. TF-target and 6. miRNA-target
#1. rna 2, miRNA expression data
#we extracted RNA and miRNA expression data from our previous ceRNA analysis

RNA_PRAD <- read.delim("Z:/GENETICS/ceRNA/ceRNA_GDCRNATools/rnaExpr_PRAD.txt", row.names=1)
miRNA_PRAD <- read.delim("Z:/GENETICS/ceRNA/ceRNA_GDCRNATools/mirExpr_PRAD.txt", row.names=1)
 
# 3. methylation (Gene level, HM450K,Beta values, Illumina HM450K platform) and 4. copy number alteration data (gene level, GISTIC2 log ratio format) were downloaded
# from "LinkedOmics" website , http://www.linkedomics.org/data_download/TCGA-PRAD/
# In, methyl and cna data we couldnot notice sample type (01/11). Then we add that prefix for each column of methyl and cna data
col1<-paste(colnames(methyl_PRAD),"01",sep = ".")
colnames(methyl_PRAD)<-col1
col2<-paste(colnames(cna_PRAD),"01",sep = ".")
colnames(cna_PRAD)<-col2
 
# 5. TF-target interactions were extrcated from "dorothea_hs" in dorothea R package
tf.target.interactions_PRAD<-dorothea_hs[,c(1,3)]
tf.target.interactions_PRAD<-as.data.frame(tf.target.interactions_PRAD)
#We had original tf.target.intercations dataset comes with Cancerin pipeline sampledata. we tried to merge both datasets together for a unique rows
r_bind<-rbind(tf.target.interactions,tf.target.interactions_PRAD) ###we had 555317 rows
dim(r_bind)
r_bind<-setkey(r_bind,NULL)
r_bind<-unique(r_bind)
dim(r_bind) ###now 549863 rows


# 6. miRNA-target interactions were gathered from "multimiR" R package, here we can use "validated"/"predicted" options in "table" option
library(multiMiR)
mir<-get_multimir(url = "http://multimir.org/",mirna = rownames(miRNA_PRAD),table = "predicted",use.tibble = FALSE,limit = 100,legacy.out = TRUE) ###putting legacy.out=TRUE will result "S3" object
attributes(mir)
#$names
#[1] "validated"    "predicted"    "disease.drug" "queries"      "summary"     

#$class
#[1] "mmquery"

#$tables
#[1] "validated"

#$org
#[1] "hsa"

#$predicted.cutoff
#numeric(0)

#$predicted.cutoff.type
#[1] "p"

#$predicted.site
#[1] "conserved"

####Here we take first attribute "validated". Below I have shown head of miRNA-target results
head(mir[[1]][,c(3,4)])
#mature_mirna_id target_symbol
#1  hsa-miR-23a-3p      C21orf33
#2  hsa-miR-26a-5p         SMAD1
#3  hsa-miR-26a-5p         SMAD1
#4  hsa-miR-23a-3p        CXCL12
#5  hsa-miR-23a-3p        CXCL12
#6  hsa-miR-23a-3p        CXCL12

miRNA.target.interactions_PRAD<-mir[[2]][,c(3,4)]
col3<-colnames(miRNA.target.interactions)
colnames(miRNA.target.interactions_PRAD)<-col3
View(miRNA.target.interactions_PRAD)

#We found intersetced column list for all four data types
#Intersected_colnames<-Reduce(intersect,list(colnames(RNA_PRAD),colnames(miRNA_PRAD),colnames(methyl_PRAD),colnames(cna_PRAD)))

#Then we filtered only shared colum list for all four datasets

RNA_PRAD_common<-RNA_PRAD[,colnames(RNA_PRAD) %in% Intersected_colnames]
miRNA_PRAD_common<-miRNA_PRAD[,colnames(miRNA_PRAD) %in% Intersected_colnames]
methyl_PRAD_common<-methyl_PRAD[,colnames(methyl_PRAD) %in% Intersected_colnames]
cna_PRAD_common<-cna_PRAD[,colnames(cna_PRAD) %in% Intersected_colnames]

#These ensembl ID should be converted into symbol

library(EnsDb.Hsapiens.v79)
geneID<-ensembldb::select(EnsDb.Hsapiens.v79,keys =Genes,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
dim(geneID)
#60433 2 here some genes are not included in conversion. dim(RNA) was 60483 2

RNA<-subset(RNA,rownames(RNA) %in% geneID$GENEID)
dim(RNA) #now dim is 60433 485
#Now we have to tally extrcated rownames of RNA with symbol dataset
#Let's put rownames into columns
RNA<-setDT(RNA,keep.rownames=TRUE) 
head(RNA)
#Now we take intersection (left_join of two datasets)
library(dplyr)
RNA<-dplyr::left_join(geneID,RNA,by=c("GENEID"="rn"))
RNA<-RNA[,-2]
dim(RNA)
install.packages("textshape")
library(textshape)
column_to_rownames(RNA,"SYMBOL")#we could not convert as repeated values were in symbol
RNA2<-RNA[!duplicated(RNA$SYMBOL),] ###removed duplicates
RNA3<-column_to_rownames(RNA2,"SYMBOL")


