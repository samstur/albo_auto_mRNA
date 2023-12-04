
########################## Installations ###########################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install("apeglm")
BiocManager::install("GOplot")
BiocManager::install("mygene")
install.packages("RColorBrewer")
install.packages("ggolot2")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("viridis")
install.packages("vsn")
install.packages("ggthemes")
install.packages("VennDiagram")
BiocManager::install("genefilter")
BiocManager::install("ggrepel")


########################## Load Libraries ###########################
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(GOplot) 
library("apeglm")
library(mygene)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(DESeq2)
library(vsn)
library(VennDiagram)
library(genefilter)
library(ggrepel)

########################## Input HTSeq data files ###########################
#Choose directory with htseq-count data
directory<-("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Auto_Biting_RNASeq/data")
#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- list.files(directory)
head(sampleFiles)
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleNames <- substr(sampleNames, 1, nchar(sampleNames)-4) # this keeps only the treatment plus the replicate
head(sampleNames)
sampleConditions <- substr(sampleFiles, 1, 1)#to get conditions I'm pulling the first letter, which is either A (auto) or M (Manassas, anautogenous)
head(sampleConditions)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions) ## condition is either A for auto or M for Manassas (non-autogenous)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)
View(sampleTable)

########################## Make the DESeq dataset from this HTSeq count data ###########################

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds


########################## Pre-filtering ###########################
#DESeq recommends a pre-filtering step to reduce memory size and increase speed. 
#They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

########################## Re-level the reference condition ###########################
#Relevel the condition for what to compare to. 
## (Default is first condition alphabetically)
dds$condition <- relevel(dds$condition, ref = "A") # setting reference level as the autogenous group
head(dds$condition)

########################## Look at normalized data ###########################
# this isn't used downstream in the DE analysis because the DESeq() function does the normalization automatically behind the scenes
# instead we can use these normalized counts for downstream visualization if we want
dds_counts <- estimateSizeFactors(dds)
sizeFactors(dds_counts) # tells us the normalization factors for each sample
dds_counts <- counts(dds_counts, normalized = TRUE) #save this to an excel file to look at normalized counts for visualization

########################## Quality control then PCA Visualization ##########################
cds <- estimateSizeFactors(dds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds, blind=TRUE)
theme_set(theme_bw())
meanSdPlot(assay(vsd))
plotDispEsts(cds)
nudge <- position_nudge(y = 10)
PCA_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE,ntop=500)
View(PCA_data)
percentVar <- round(100 * attr(PCA_data, "percentVar"))
nudge <- position_nudge(y = 4)
pca_plot <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = condition, position_nudge(y=10))) +
  geom_point(size =3, position = position_jitter(w=0.05, h=0.05)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("dark blue", "dark red"),labels=c("Autogenous","Manassas 
Anautogenous"))+
  labs(color = "Population")
 # geom_text_repel()
# By default plotPCA() uses the top 500 most variable genes. 
# You can change this by adding the ntop= argument and specifying how many of the genes you want the function to consider.
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Auto_Biting_RNASeq/R_output_MvA")
ggsave("MvA_pcaplot.png",plot=pca_plot,dpi=600,units='in',width=6,height=5)


########################## Actually run DESeq Analysis ###########################
data <- DESeq(dds)

## Plot dispersion estimates
#   You expect your data to generally scatter around the curve, with the dispersion
#   decreasing with increasing mean expression levels. If you see a cloud or 
#   different shapes, then you might want to explore your data more to see if 
#   you have contamination (mitochondrial, etc.) or outlier samples. 
plotDispEsts(data)
#   our plot looks okay because we see the decrease of dispersion values as you get higher mean values
# https://github.com/hbctraining/DGE_workshop_salmon_online/blob/master/lessons/04b_DGE_DESeq2_analysis.md

########################## Hypothesis Testing ###########################
# DESeq2 uses a negative binomial distribution to model RNA-seq counts since they exhibit overdispersion
# (variance > mean). 
# It is common to shrink the log fold change estimate for better visualization and ranking of genes.
#Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(data)
res_LFC <- lfcShrink(data, coef="condition_M_vs_A", type="apeglm") ## plug in the output from resultsNames(dds) line as the coef
res_LFC
head(res_LFC)
# Order the table by smallest p value
resOrdered <- res_LFC[order(res_LFC$padj),]
summary(resOrdered)
head(resOrdered)
View(resOrdered)

# Write out a table of all genes for KEGG enrichment
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Auto_Biting_RNASeq/R_output_MvA")
write.csv(resOrdered, 
          file="./MvA_allgenes.csv", row.names = T)

# Make a volcano plot
## visualizer each results set with a volcano plot
VP <- EnhancedVolcano(res_LFC,
                      lab = rownames(res_LFC),
                      x = 'log2FoldChange',
                      y = 'padj',
                      ylab = "-Log10(p-adjusted)",
                      selectLab = NA,
                      #drawConnectors = TRUE,
                      xlim = c(-2.5, 2.5),
                      ylim = c(0,25),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      pointSize = 2.0,
                      labSize = 5.0)
VP

ggsave("volcano_MvA.png",plot=VP,dpi=600,units='in',width=8,height=6)

sig_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
View(sig_res)
length(unique(sig_res$gene)) #793 


# adding gene names to this deg
genes <- queryMany(sig_res$gene, scopes="symbol", fields=c("name"))
colnames(genes)<-c("gene","id","score","gene name")
genes <- as.data.frame(subset(genes, select = c("gene", "gene name")))
View(genes)
length(unique(genes$gene)) #793
length(genes$gene) #802
genes <- (genes[!duplicated(genes), ]) ##getting rid of duplicate gene name rows for Trnad-guc
sig_res.gene <- as.data.frame(merge(genes,sig_res,by="gene"))
View(sig_res.gene) 
length(unique(sig_res.gene$gene)) #793
length(sig_res.gene$gene) #793

# Write out a table of these significant differentially expressed genes
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Auto_Biting_RNASeq/R_output_MvA")
write.csv(sig_res.gene, 
          file="./MvA_sig.gene.csv", row.names = F)
