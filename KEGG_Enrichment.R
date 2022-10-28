# KEGG analysis
# The pathway enrichment analysis for this I modified from "https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html"

# BiocManager::install("KEGGREST")

# Load libraries
install.packages("KEGGREST")
library(KEGGREST)

# Set some variables
setwd("~/Georgetown/Armbruster Lab")
genes_with_pvalues <- "./Auto_Biting_RNASeq/AvM_LFCshrink_padj.txt" # with at least columns "gene", "padj"

keggGeneID_output <- "./Auto_Biting_RNASeq/AvM_keggID.txt"

# Read in data
gene_list <- read.csv(genes_with_pvalues, header = T)
View(gene_list)

# Make a list of all the ncbi to albopictus gene ids
convs <- keggConv("ncbi-geneid", "aalb")
head(convs)
#convs2 <- keggConv("aalb", "uniprot")
#convs3 <- keggConv("aalb", "ncbi-proteinid")

# Convert gene ids from list to something relevant to our work. Note that all our gene ids begin with "LOC". For NCBI LOC1234 is equivalent to GeneID = 1234. The "LOC"+GeneID is when orthologs have not yet been determined and a published symbol is not available. 
gene_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", gene_list$gene)
head(gene_list)

# Find the matching ncbi id in the conversion list and take the name associated (the kegg id) and assign that in the kegg_id column of the gene_list dataframe
gene_list$kegg_id = names(convs)[match(gene_list$ncbi_geneid, as.character(convs))]
#gene_list$uniprot_id = names(convs2)[match(gene_list$kegg_id, as.character(convs2))]
#gene_list$protein_id = names(convs3)[match(gene_list$kegg_id, as.character(convs3))] # getting the protein id does nothing for us downstream as these are not searchable in the fasta file. I think this is likely due to protein "versions" (basically all the protein ids in the fasta end with a ".#" following so I think they are alternative proteins for the genes)
head(gene_list)
# If you want to write out the KEGG ID list and do this online
write.table(gene_list$kegg_id, 
            file=keggGeneID_output, col.names = F, row.names = F, quote = F)
###
# You can input this KEGG list at "https://www.kegg.jp/kegg/tool/map_pathway1.html" to find enriched pathways
###


########
######
####
#Trying to automate this process with KEGG pathway enrichment
#all_genes_list <- read.csv(file="../misc/DESeq_results_pharatelarvae.csv")

all_genes_list <- read.csv(file="AVM_allgenes.csv")
head(all_genes_list)
#all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$geneID)
all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$gene)
head(all_genes_list)
all_genes_list$kegg_id = names(convs)[match(all_genes_list$ncbi_geneid, as.character(convs))]


# Get the pathways list from KEGG
pathways.list <- keggList("pathway", "aalb")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
head(pathway.codes)
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

#geneList <- all_genes_list$X11d_padj #CHANGE THIS LINE
geneList <- all_genes_list$padj #CHANGE THIS LINE
head(geneList)
#geneLFClist <- all_genes_list$X11d_Log2FoldChange #CHANGE THIS LINE
geneLFClist <- all_genes_list$log2FoldChange 

#names(geneList) <- sub("aalb:","", all_genes_list$kegg_id) # get rid of the beginning "aalb:" since the gene list we brought from kegg doesn't have this
names(geneList) <- sub("aalb:","", all_genes_list$kegg_id) # get rid of the beginning "aalb:" since the gene list we brought from kegg doesn't have this
head(geneList)
View(all_genes_list)
View(geneList)
names(geneLFClist) <- sub("aalb:","", all_genes_list$kegg_id)
head(geneLFClist)

#genes.by.pathway_40d <- genes.by.pathway[-c(99, 120)]

pathway_pval <- data.frame()

for (pathway in 1:length(genes.by.pathway)){
  pathway.genes <- genes.by.pathway[[pathway]]
  if (!is.na(pathway.genes)){
    list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
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
                 length(list.genes.in.pathway), #number of genes from the DE list in the pathway 
                 #(including DE genes that didn't meet pvalue and logFC thresholds)
                 sum(abs(geneLFClist[list.genes.in.pathway])>0.58&geneList[list.genes.in.pathway]<0.05), #number of sig DE genes in pathway 
                 sum(geneLFClist[list.genes.in.pathway]>0.58&geneList[list.genes.in.pathway]<0.05), # number of sig up DE genes in pathway
                 sum(geneLFClist[list.genes.in.pathway]< -0.58&geneList[list.genes.in.pathway]<0.05),
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
          file="./AvM_keggPathwayEnrichment_allgenes.csv", row.names = F)




# Address the magnitude and direction of log2 fold-change for specific pathways
# loop version of Sarah's original code to pull out all sig DE genes in each 
# enriched pathway
library("tidyverse")
# first load in the pathways that we're interested in. I've put all the ones that came up as sig enriched
#   from KEGG analysis above
list_interesting_pathways <- c("aalb04142","aalb00500",
                               "aalb00250","aalb00270","aalb00520","aalb00910",
                               "aalb00051","aalb00040","aalb00220", "aalb04141",
                               "aalb04213", "aalb00785", "aalb00531", "aalb00620",
                               "aalb00290","aalb00790","aalb00980","aalb00240",
                               "aalb00010","aalb00563")
all_genes_list$no_loc_geneID <- sub("LOC","", all_genes_list$gene) # makes new column in all_genes_list that has the geneID without the LOC at the start

df <- data.frame(geneID="NA",baseMean=NA,log2FoldChange=NA,lfcSE=NA,pvalue=NA,padj=NA,ncbi_geneid=NA,
                 kegg_id="NA",no_loc_geneID=NA,pathway_code=NA) ##makes blank dataframe to add to
for (pathway in list_interesting_pathways){ # for each pathway, pulls out all genes in both the pathway and the DE gene list
  pathway_genes <- genes.by.pathway[[pathway]]
  pathway_genes_data <- all_genes_list[all_genes_list$no_loc_geneID %in% pathway_genes,]
  pathway_genes_sig <- pathway_genes_data %>%
    filter(padj < 0.05, abs(log2FoldChange) > 0.58) # only keeps ones that are sig DE
  pathway_genes_sig$pathway_code <- pathway
  newdata = data.frame(pathway_genes_sig)
  df <- rbind(df,newdata)
}
head(df)
Finaldf<-df[!(df$geneID=="NA"),]
head(Finaldf)
## adding pathway description
Finaldf$pathwayName = pathways.list[match(Finaldf$pathway_code, sub("path:","", names(pathways.list)))]

## adding gene names
genes <- "./Albo_DESeq_results_with_genenames.csv" # with at least columns for LOC id and gene name
gene_df <- read.csv(genes, header = T)
head(gene_df)
head(match(Finaldf$geneID, gene_df$geneID))
Finaldf$gene_name = gene_df$name[match(Finaldf$geneID, gene_df$geneID)]

# reordering column positions
col_order <- c("pathway_code", "pathwayName", "geneID","gene_name","baseMean","log2FoldChange",
               "lfcSE","pvalue","padj","ncbi_geneid","kegg_id","no_loc_geneID")
Finaldf_2 <- Finaldf[, col_order]

View(Finaldf_2)

write.csv(Finaldf_2, 
          file="./AvM_DEgenes_in_EnrichedKEGGPathways.csv", row.names = F)



## same but for B v NB data
file <- read.csv("BVN_KEGG_pathway.csv", header = T)
head(file)
list_interesting_pathways <- file$pathwayCode
all_genes_list$no_loc_geneID <- sub("LOC","", all_genes_list$gene) # makes new column in all_genes_list that has the geneID without the LOC at the start

df <- data.frame(geneID="NA",baseMean=NA,log2FoldChange=NA,lfcSE=NA,pvalue=NA,padj=NA,ncbi_geneid=NA,
                 kegg_id="NA",no_loc_geneID=NA,pathway_code=NA) ##makes blank dataframe to add to
for (pathway in list_interesting_pathways){ # for each pathway, pulls out all genes in both the pathway and the DE gene list
  pathway_genes <- genes.by.pathway[[pathway]]
  pathway_genes_data <- all_genes_list[all_genes_list$no_loc_geneID %in% pathway_genes,]
  pathway_genes_sig <- pathway_genes_data %>%
    filter(padj < 0.05, abs(log2FoldChange) > 0.58) # only keeps ones that are sig DE
  pathway_genes_sig$pathway_code <- pathway
  newdata = data.frame(pathway_genes_sig)
  df <- rbind(df,newdata)
}
head(df)
Finaldf<-df[!(df$geneID=="NA"),]
head(Finaldf)
## adding pathway description
Finaldf$pathwayName = pathways.list[match(Finaldf$pathway_code, sub("path:","", names(pathways.list)))]

## adding gene names
genes <- "./Albo_DESeq_results_with_genenames.csv" # with at least columns for LOC id and gene name
gene_df <- read.csv(genes, header = T)
head(gene_df)
head(match(Finaldf$geneID, gene_df$geneID))
Finaldf$gene_name = gene_df$name[match(Finaldf$geneID, gene_df$geneID)]

# reordering column positions
col_order <- c("pathway_code", "pathwayName", "geneID","gene_name","baseMean","log2FoldChange",
               "lfcSE","pvalue","padj","ncbi_geneid","kegg_id","no_loc_geneID")
Finaldf_2 <- Finaldf[, col_order]

View(Finaldf_2)

write.csv(Finaldf_2, 
          file="./BvN_DEgenes_in_EnrichedKEGGPathways.csv", row.names = F)










## sarah's orginal code to do this pathway by pathway
interesting_pathway <- "aalb00220"

pathway_genes <- genes.by.pathway[[interesting_pathway]]
head(pathway_genes)

head(all_genes_list)
all_genes_list$no_loc_geneID <- sub("LOC","", all_genes_list$gene)

pathway_genes_data <- all_genes_list[all_genes_list$no_loc_geneID %in% pathway_genes,]

library("tidyverse")
pathway_genes_sig <- pathway_genes_data %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)
View(pathway_genes_sig)
pathway_genes_sig$pathway_code <- interesting_pathway
