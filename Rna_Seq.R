setwd("E:/cancer genomics/cancer genomics")
if(!requireNamespace("BioManager",quitely=TRUE))
  install.packages("BioManager")
  BiocManager::install("3.18")

BiocManager::version()

install.packages("R.utils")
library(R.utils)#for coping files, readlines,checking permission


BiocManager::install("Rsubread")
#This package is used for Ngs data analysis where it can have functions that align the reads to
#reference genome,quantification,tools for differential gene expression,visualization
library(Rsubread) 
# 
search()
sessionInfo()
getRversion()
R.version

install.packages("data.table")
library(data.table)

BiocManager::install("RUVSeq")
library(RUVSeq)

search()
BiocManager::install("DESeq2")
library(DESeq2)

#this file tells the conidtion whether normal or cancer for each patiet id
data=read.table("TCGA-column-data (1).csv",header=T,sep=",",row.names = 1,check.names = FALSE)
head(data)
dim(data)
# count of genes in each patient

count_data=read.table("TCGA-count-data (1).csv",header=T,sep=",",row.names = 1,check.names = FALSE)
head(count_data)
View(count_data)
dim(count_data)
all(row.names(data)==colnames(count_data))
?read.table
class(count_data)
class(data)
is.na(data)
(is.na(count_data))
which(is.na(count_data),arr.ind = TRUE) # tells which row

sum(is.na(count_data)) # count total number of missing values
library(zoo)
install.packages("zoo")
count_data=t(na.aggregate(t(count_data)))
#t- transposses
which(is.na(count_data),arr.ind = TRUE)

des<-DESeqDataSetFromMatrix(countData = round(count_data),colData = data,design=~condition)
colnames(des)
row.names(des)
class(des)
#to check any non inte_ger values
which(is.integer(count_data))
any(count_data!=is.integer(count_data))

des$condition
data
library(dplyr)
group_by(data,condition)%>%count() # number in each clss
str(data)
#count_per_class <- data %>%group_by(condition) %>% summarize(count = n())

des$condition<-relevel(des$condition,ref='normal')
#print(count_per_class)

des<-DESeq(des)
result<-results(des)
summary(result)

class(result)
 
#  subsetting the upregulated and downregulated genes 
up_reg<-subset(result,padj<0.05 & log2FoldChange>1)
write.csv(up_reg,"upregulated_gene.csv")

down_reg<-subset(result,padj<0.05 & log2FoldChange< -1)
write.csv(up_reg,"downregulated_gene.csv")

de<-subset(result,padj<0.05 & log2FoldChange>1|padj<0.05 & log2FoldChange< -1)
write.csv(de,"de.csv")

re=read.csv("de.csv")
head(re)

write.table(result, file="Volcano", row.names=F, sep="\t") 


library(EnhancedVolcano)
BiocManager::install("EnhancedVolcano")



EnhancedVolcano(result,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL,
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )
dim(result)
EnhancedVolcano(result,
               lab = rownames(result),
               x = 'log2FoldChange',
               y = 'pvalue',
               selectLab = c('ENSG00000287048.1','ENSG00000280032.1'),
               title = NULL,
               cutoffLineType = 'twodash',
               cutoffLineWidth = 0.8,
               xlim = c(-8,8),
               xlab = bquote(~Log[2]~ 'fold change'),
               ylim = c(0,12),
               ylab = bquote(~-Log[10]~italic(P)),
               pCutoff = 0.05,
               FCcutoff = 1.0,
               #transcriptPointSize = 0.5,
               #transcriptLabSize = 4.0,
               colAlpha = 1,
               shape = 19,
               subtitle = NULL,
               legendPosition = 'top',
               legendLabSize = 12,
               legendIconSize = 4.0,
               gridlines.major = FALSE,
               gridlines.minor = FALSE,
               drawConnectors = TRUE,
               widthConnectors = 1.0,
               colConnectors = 'black',
               border = 'full' )




# patients wh
head(up_reg)
View(up_reg)
row.names(up_reg)
all(row.names(data)==colnames(count_data))

row.names(up_reg) #genes names
row.names(count_data)
colnames(count_data)
column<-colData(des)



dim(up_reg)
deseq_results_up <- count_data[rownames(count_data) %in% rownames(up_reg), ]
View(deseq_results_up)
dim(deseq_results_up)
deseq_results_down<-count_data[rownames(count_data) %in% rownames(down_reg),]
dim(down_reg)
dim(deseq_results_down)
