#Volcano plot with EnhancedVolcano (https://github.com/kevinblighe/EnhancedVolcano)
#  BiocManager::install('EnhancedVolcano')
#with results<- a dataframe containing the gene ids, log2FoldChange, adjusted-pvalues (padj)
#label_volcano<-a dataframe with the gene ids to label.
#Sec<-a dataframe with the gene ids of the genes encoding secreted proteins
#keyval.shape<-a numeric to assign a specific point shape to each category of interest (Sec/DEGs for instance)


#Here I make a volcano plot where I (1) label only the genes encoding effectors and (2) modify the size and point shape of the gene encoding secreted proteins. 
#shape 9 is a diamond, shape 16 a normal plain dot 
keyvals.shape_sec <- ifelse(
rownames(results) %in% unique(Sec$Gene), 9,
16)
names(keyvals.shape_sec)[keyvals.shape_sec == 9] <- 'Sec'
names(keyvals.shape_sec)[keyvals.shape_sec == 16] <- 'DEGs'

#I already have loaded in R the df Sec, results (can use the DEseq2 contrast alternatively) and label_volcano.

EnhancedVolcano(results, lab=results$Gene, selectLab = unique((label_volcano$Gene)), x = 'log2FoldChange', y = 'padj',xlim = c(-2.5, 2.5),xlab = bquote(~Log[2]~ 'fold change'),ylab = bquote(~-Log[10]~adjusted~italic(P)),FCcutoff = 0.5,pCutoff = 0.001,colAlpha = 0.6,pointSize = c(ifelse(results$Gene%in%Sec$Gene, 3, 1)),gridlines.major = FALSE,gridlines.minor = FALSE,vline=0.5,boxedLabels = FALSE,labCol = "black",labFace = "bold",title = "Regulated genes",shapeCustom = keyvals.shape_sec)

setwd("C:/Users/pjt6/OneDrive - University of St Andrews/Bat_FOXP2_overexpression_Dec_2021_final")
library(EnhancedVolcano)
library(magrittr)
in_data = read.table("data.txt", header=T,stringsAsFactors=F)
genes = in_data$gene
padj = in_data$padj
# simple  log2FoldChange	padj
 EnhancedVolcano(in_data,
    lab = rownames(in_data),
    x = 'log2FoldChange',
    y = 'pvalue', 
    xlim = c(-5,5), 
    ylim = c(0,4),
    title = 'GFP_striatum_vs_KI_striatum', 
    pCutoff = 0.05, 
    FCcutoff = 1.0, 
    pointSize = 1.5, labSize = 2.0,
    legendLabSize = 10,
    colAlpha = 1)

    
EnhancedVolcano(in_data, lab=in_data$gene, selectLab = unique((label_volcano$Gene)), x = 'log2FoldChange', y = 'padj',xlim = c(-2.5, 2.5),xlab = bquote(~Log[2]~ 'fold change'),ylab = bquote(~-Log[10]~adjusted~italic(P)),FCcutoff = 0.5,pCutoff = 0.001,colAlpha = 0.6,pointSize = c(ifelse(results$Gene%in%Sec$Gene, 3, 1)),gridlines.major = FALSE,gridlines.minor = FALSE,vline=0.5,boxedLabels = FALSE,labCol = "black",labFace = "bold",title = "Regulated genes",shapeCustom = keyvals.shape_sec)
