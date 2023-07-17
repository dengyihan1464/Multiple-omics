library(DESeq2)
library(clusterProfiler)
library(pheatmap)
library(enrichplot)

# load data-------------------------------------
dataDir <- '/Sshare/home/lanLessonPublic/RNA-seq/5.featureCounts/'
counts_data <- read.table(file = paste0(dataDir, 'Counts.txt'), header = T, sep = '\t', stringsAsFactors = F)
Geneid <- counts_data$Geneid
counts_data <- counts_data[, -grep('Geneid', colnames(counts_data))]
colnames(counts_data) <- c('HBR_Rep1', 'HBR_Rep2', 'HBR_Rep3', 'UHR_Rep1', 'UHR_Rep2', 'UHR_Rep3')
metaData <- data.frame(id = colnames(counts_data), group = c(rep('HBR', 3), rep('UHR', 3)))
rownames(counts_data) <- Geneid
zero_index <- apply(counts_data, 1, FUN = function(x){
  return(length(which(x == 0)))
})

zero_index <- which(zero_index == ncol(counts_data))
counts_data <- counts_data[-zero_index, ]

# DEGs----------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts_data, 
                              colData=metaData, 
                              design=~group)
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy = T)) 
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- res[which(res$padj<0.05), ]
res <- res[order(res$log2FoldChange, decreasing = T), ]

# plot
plot_data <- rownames(res)[which(abs(res$log2FoldChange)>2)]
plot_data <- as.matrix(counts_data[which(rownames(counts_data) %in% plot_data), ])
pheatmap(plot_data, angle_col = 45, scale = 'row', show_rownames = F)

# GSEA GO----------------------------------------------------------------
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
gse <- gseGO(geneList=geneList,
             ont ="ALL",
             keyType = "ENSEMBL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = 'org.Hs.eg.db',
             pAdjustMethod = "none")

p1 <- dotplot(gse, showCategory=30) + ggtitle("dotplot for GSEA")
print(p1)
gseaplot2(gse, geneSetID=c(15, 25, 37))
gse <- as.data.frame(gse)
gse <- gse[order(gse$NES, decreasing = T), ]

# enrichGO---------------------------------------------------------------------
geneList <- rownames(res)[which(abs(res$log2FoldChange)>2)]
enrich <- enrichGO(gene=geneList,
                   ont ="ALL",
                   keyType = "ENSEMBL",
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   OrgDb = 'org.Hs.eg.db',
                   pAdjustMethod = "none")
p2 <- barplot(enrich, showCategory=20)
print(p2)
