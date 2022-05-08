
#Author: ReyWu

library(tidyverse)
library(DESeq2)
library(dplyr)
library(psych) 
library(formattable)
library(RColorBrewer)
library(pheatmap)
library(gplots)

#---------------------------------------------------raw data------------------------------------------------
rawCounts <- read.delim("/Users/rey/Desktop/bf591_R/FinalProject_reywu/data/raw-counts.tsv")
rawCounts <- rawCounts[!(is.na(rawCounts$Gene.Name) | rawCounts$Gene.Name==""), ]
rawCounts <- rawCounts[1:20526,]
write.csv(rawCounts, paste("data/rawCount_table.csv"), row.names = FALSE)
#------------------------------------------------------------------------------------------------------------


#----------------------------------------------------norm data------------------------------------------------
sampleData <- read.delim("/Users/rey/Desktop/bf591_R/FinalProject_reywu/data/experiment-design.tsv")
sampleData_1 <- read.delim("/Users/rey/Desktop/bf591_R/FinalProject_reywu/data/experiment-design.tsv")

#Convert count data to a matrix of appropriate form that DEseq2 can read
rawCounts <- rawCounts %>%
  select(!Gene.Name) %>%
  data.frame(row.names = 1)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
sampleData <- sampleData %>%
  column_to_rownames(var="Run") %>%
  select(c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual."))%>%
  dplyr::rename(tissueType="Sample.Characteristic.biopsy.site.",
         individualID = "Sample.Characteristic.individual.") 
sampleData$individualID <- factor(sampleData$individualID)
sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal", "primary tumor", "colorectal cancer metastatic in the liver"))
write.csv(sampleData, paste("data/meta_info.csv"))

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Create the DEseq2DataSet object
dds<- DESeqDataSetFromMatrix(
  countData=rawCounts, 
  colData=sampleData, 
  design= ~ individualID + tissueType
  )

#Create Normalized count matrix 
dds_norm_factor <- estimateSizeFactors(dds)
deseq_norm_counts <- as.data.frame(counts(dds_norm_factor,normalized=TRUE))
#rownames(deseq_norm_counts) <- rownames(rawCounts)
deseq_norm_counts <- mutate(deseq_norm_counts, across(where(is.numeric), ~ round(., digits = 0)))
write.csv(deseq_norm_counts, paste("data/deseq_norm_counts.csv"))
#-----------------------------------------------------------------------------------------------------------------


#------------------------------------------------------DE data----------------------------------------------------
#DE result by deseq2
dds_Data <- DESeq(dds)
res <- results(dds_Data, contrast=c("tissueType", "primary tumor", "normal"))
res_table <- as.data.frame(res)
write.csv(res_table, paste("data/deseq_DE_table.csv"))
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------fgsea result table------------------------------------------------
library(fgsea)
res <- read.csv("/Users/rey/Desktop/bf591_R/FinalProject_reywu/data/deseq_DE_table.csv")
genes <- select(rawCounts, c(Gene.ID, Gene.Name))
res<-inner_join(res, genes, by=c("X"="Gene.ID"))

res2 <- res %>% 
  dplyr::select(Gene.Name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene.Name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

pathways <- fgsea::gmtPathways('data/c2.cp.v7.5.1.symbols.gmt')

fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)
fgseaRes <- select(fgseaRes, -leadingEdge, -ES, -nMoreExtreme) 
write.csv(fgseaRes, paste("data/fgsea_Result.csv"),row.names = FALSE)
#----------------------------------------------------------------------------------------------------------------
