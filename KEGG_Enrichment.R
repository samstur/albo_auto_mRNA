
# Load libraries
install.packages("KEGGREST")
library(KEGGREST)

#### original workflow from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
# Set some variables
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Auto_Biting_RNASeq/R_output_MvA")
genes_with_pvalues <- "./MvA_allgenes.csv" # with at least columns "gene", "padj"

# Read in data
all_genes_list <- read.csv(genes_with_pvalues, header = T)
View(all_genes_list)

# Make a list of all the ncbi to gene ids using the organism code for albopictus (aalb)
convs <- keggConv("ncbi-geneid", "aalb") 
head(convs)

all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$gene)
head(all_genes_list)
# Find the matching ncbi id in the conversion list and take the name associated 
#   (the kegg id) and assign that in the kegg_id column of the gene_list dataframe
all_genes_list$kegg_id = names(convs)[match(all_genes_list$ncbi_geneid, as.character(convs))]
View(all_genes_list)

# Get the pathways list from KEGG
pathways.list <- keggList("pathway", "aalb")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
View(genes.by.pathway)

geneList <- all_genes_list$padj 
head(geneList)
geneLFClist <- all_genes_list$log2FoldChange 

names(geneList) <- sub("ncbi-geneid:","", all_genes_list$ncbi_geneid)
head(geneList)

names(geneLFClist) <- sub("ncbi-geneid:","", all_genes_list$ncbi_geneid)
head(geneLFClist)

pathway_pval <- data.frame()
for (pathway in 1:length(genes.by.pathway)){
  pathway.genes <- genes.by.pathway[pathway] ## removed extra set of brackets around pathway
  if (!is.na(pathway.genes)){
    list.genes.in.pathway <- intersect(names(geneList), pathway.genes[[1]]) ## added extra brackets here
    list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
    scores.in.pathway <- geneList[list.genes.in.pathway]
    scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
    mean.score.in.pathway <- mean(scores.in.pathway,na.rm=T) ##added line
    mean.score.not.in.pathway <- mean(scores.not.in.pathway,na.rm=T) ##added line
    lfc.in.pathway <- geneLFClist[list.genes.in.pathway]
    mean.lfc.in.pathway <- mean(lfc.in.pathway,na.rm=T) ##added line
    if (length(scores.in.pathway) > 0){
      p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
    } else{
      p.value <- NA
    }
    new_row <- c(names(genes.by.pathway[pathway]), # pathway id
                 p.value, # pvalue
                 length(list.genes.in.pathway), #number of genes from the full DESeq2 output list in the pathway 
                 #(including DESeq2 output genes that didn't meet pvalue and logFC thresholds)
                 sum(abs(geneLFClist[list.genes.in.pathway])>1&geneList[list.genes.in.pathway]<0.05), #number of sig DE genes in pathway 
                 sum(geneLFClist[list.genes.in.pathway]>1&geneList[list.genes.in.pathway]<0.05), # number of sig up DE genes in pathway
                 sum(geneLFClist[list.genes.in.pathway]< -1&geneList[list.genes.in.pathway]<0.05),
                 mean.score.in.pathway, ## new line
                 mean.score.not.in.pathway, ## new line
                 mean.lfc.in.pathway) ## new line
    pathway_pval <- rbind(pathway_pval, new_row)
  }
}


head(pathway_pval)


colnames(pathway_pval) <- c("pathwayCode", "pval", "annotated", "DEG", "up", "down","mean p of genes in pathway",
                            "mean p of genes not in pathway","mean LFC of genes in pathway")   ### changed line
pathway_pval <- pathway_pval[complete.cases(pathway_pval),]
head(pathway_pval$pathwayCode)
pathway_pval$pathwayName = pathways.list[match(pathway_pval$pathwayCode, sub("path:","", names(pathways.list)))]

head(pathway_pval)
pathway_pval$pval <- as.numeric(pathway_pval$pval)

pathway_pval <- pathway_pval[order(pathway_pval$pval),]
head(pathway_pval)
View(pathway_pval)
# Write out a csv with these data
write.csv(pathway_pval, 
          file="./MvA_keggPathwayEnrichment_allgenes_DEGlog2fc1_updated.csv", row.names = F)




# Address the magnitude and direction of log2 fold-change for specific pathways
# loop version pull out all sig DE genes in each 
# enriched pathway
library("tidyverse")
# first load in the pathways that we're interested in. I've put all the ones that came up as sig enriched
#   from KEGG analysis above and that had >= 2 DE genes 
list_interesting_pathways <- c("aalb04142","aalb00040")
all_genes_list$no_loc_geneID <- sub("LOC","", all_genes_list$gene) # makes new column in all_genes_list that has the geneID without the LOC at the start
head(all_genes_list)
df <- data.frame(gene="NA",baseMean=NA,log2FoldChange=NA,lfcSE=NA,pvalue=NA,padj=NA,ncbi_geneid=NA,
                 kegg_id="NA",no_loc_geneID=NA,pathway_code=NA) ##makes blank dataframe to add to
for (pathway in list_interesting_pathways){ # for each pathway, pulls out all genes in both the pathway and the DE gene list
  pathway_genes <- genes.by.pathway[[pathway]]
  pathway_genes_data <- all_genes_list[all_genes_list$no_loc_geneID %in% pathway_genes,]
  pathway_genes_sig <- pathway_genes_data %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) # only keeps ones that are sig DE
  pathway_genes_sig$pathway_code <- pathway
  newdata = data.frame(pathway_genes_sig)
  df <- rbind(df,newdata)
}
head(df)
Finaldf<-df[!(df$gene=="NA"),]
head(Finaldf)
## adding pathway description
Finaldf$pathwayName = pathways.list[match(Finaldf$pathway_code, sub("path:","", names(pathways.list)))]

View(Finaldf)

## adding gene names

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Auto_Biting_RNASeq/R_output_MvA")
gene_df <- read.csv("MvA_sig.gene.csv", header = T)
head(gene_df)
match(Finaldf$gene, gene_df$gene)
Finaldf$gene_name <- gene_df$gene.name[match(Finaldf$gene, gene_df$gene)]
head(Finaldf)
# reordering column positions
col_order <- c("pathway_code", "pathwayName", "gene","gene_name","baseMean","log2FoldChange",
               "lfcSE","pvalue","padj","ncbi_geneid","kegg_id","no_loc_geneID")
Finaldf_2 <- Finaldf[, col_order]
Finaldf_3 <- Finaldf_2[1:9]
View(Finaldf_3)

write.csv(Finaldf_3, 
          file="./MvA_DEgenes_in_EnrichedKEGGPathways.csv", row.names = F)
