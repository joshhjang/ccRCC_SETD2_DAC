#R version 4.0.2 (2020-06-22)
library(ggplot2) #3.3.5
library(RColorBrewer) #1.1-2
library(plyr) #1.8.6
library("ComplexHeatmap") #2.6.2
library("circlize") #0.4.13
library(reshape2) #1.4.4
library(DESeq2)#1.30.1
library(matrixStats) #0.60.1

mypalette <- brewer.pal(12,"Paired")
mypalette2 <- brewer.pal(12,"Set3")
mypalette3 <- brewer.pal(8,"Pastel1")
mypalette4 <- brewer.pal(8,"Set2")
mypalette5 <- brewer.pal(11,"Spectral")
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


###meta data###
setwd("~/Desktop/Projects/SETD2_Minmin_Gangning/Submission_code")
meta <- read.table("Cell_Sample_info_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)
meta <- meta[order(meta$Type,meta$Cell,meta$Day,meta$Rep),]
meta1 <- meta
meta1$Cell <- factor(meta1$Cell, levels = c("ACHN","O786","CAKI","A498"))
meta1 <- meta1[order(meta1$Cell,meta1$Type,meta1$Day,meta1$Rep),]#order cell type based on SETD2 mutation status

############# DESeq2: check differentially expressed genes ###############
setwd("./DESeq2")
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))#from Stringtie prepDE.py
colData <- read.csv("SETD2_Sample_Metadata_forDESeq2.txt", sep="\t", row.names=1)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ CellDayType)
keep <- rowSums(counts(dds)) >= 1 #filter lowly expressed genes
dds <- dds[keep,]
dds <- DESeq(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("Cell","Day","Rep","CellDayType")])

#Make Sample-to-sample Distance Heatmap (Fig S1A)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( c(brewer.pal(5, "Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#Make PCA plot of RNAseq data (Fig S2). Change colors manually on illustrator.
plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}

a<-plotPCA(vsd, intgroup=c("CellDayType"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a


### call DEGs for each timepoint comparison (LFC >1, adj pvalue cutoff < 0.01)  ###
res_O786EV_D0D5 <- results(dds,contrast=c("CellDayType","O786D0EV","O786D5EV"),lfcThreshold=1,alpha=0.01)
res_O786EV_D0D15 <- results(dds,contrast=c("CellDayType","O786D0EV","O786D15EV"),lfcThreshold=1,alpha=0.01)
res_O786EV_D0D23 <- results(dds,contrast=c("CellDayType","O786D0EV","O786D23EV"),lfcThreshold=1,alpha=0.01)

res_O786KO_D0D5 <- results(dds,contrast=c("CellDayType","O786D0KO","O786D5KO"),lfcThreshold=1,alpha=0.01)
res_O786KO_D0D15 <- results(dds,contrast=c("CellDayType","O786D0KO","O786D15KO"),lfcThreshold=1,alpha=0.01)
res_O786KO_D0D23 <- results(dds,contrast=c("CellDayType","O786D0KO","O786D23KO"),lfcThreshold=1,alpha=0.01)

res_O786EV_KO_D0 <- results(dds,contrast=c("CellDayType","O786D0EV","O786D0KO"),lfcThreshold=1,alpha=0.01)
res_O786EV_KO_D5 <- results(dds,contrast=c("CellDayType","O786D5EV","O786D5KO"),lfcThreshold=1,alpha=0.01)
res_O786EV_KO_D15 <- results(dds,contrast=c("CellDayType","O786D15EV","O786D15KO"),lfcThreshold=1,alpha=0.01)
res_O786EV_KO_D23 <- results(dds,contrast=c("CellDayType","O786D23EV","O786D23KO"),lfcThreshold=1,alpha=0.01)

res_A498_D0D5 <- results(dds,contrast=c("CellDayType","A498D0Mut","A498D5Mut"),lfcThreshold=1,alpha=0.01)
res_A498_D0D15 <- results(dds,contrast=c("CellDayType","A498D0Mut","A498D15Mut"),lfcThreshold=1,alpha=0.01)

res_ACHN_D0D5 <- results(dds,contrast=c("CellDayType","ACHND0WT","ACHND5WT"),lfcThreshold=1,alpha=0.01)
res_ACHN_D0D14 <- results(dds,contrast=c("CellDayType","ACHND0WT","ACHND14WT"),lfcThreshold=1,alpha=0.01)
res_ACHN_D0D22 <- results(dds,contrast=c("CellDayType","ACHND0WT","ACHND22WT"),lfcThreshold=1,alpha=0.01)

res_CAKI_D0D12 <- results(dds,contrast=c("CellDayType","CAKID0Mut","CAKID12Mut"),lfcThreshold=1,alpha=0.01)
res_CAKI_D0D25 <- results(dds,contrast=c("CellDayType","CAKID0Mut","CAKID25Mut"),lfcThreshold=1,alpha=0.01)



######  fGSEA Hallmark and Reactome pathways ###
setwd("../fGSEA")
library(fgsea) #1.16.0
library(tibble) #0.3.5
library(tidyr) #1.1.3
library(dplyr) #1.0.7

Ensb_genelist<-read.table("ENSEMBL_gene_transcript_ID_Genename.txt", sep = "\t", header =T, stringsAsFactors = F) #table containing Gene name associated with ENSEMBL label
colnames(Ensb_genelist)[4] <- "GeneID"
ens2symbol <- as_tibble(Ensb_genelist)
colnames(ens2symbol)[4] <- "GeneID"

pathways.hallmark <- gmtPathways("../GSEA gene set lists/h.all.v7.1.symbols.gmt") #Hallmark
pathways.reactome <- gmtPathways("../GSEA gene set lists/c2.cp.reactome.v7.1.symbols.gmt")# Reactome


#iterate through each pairwise comparison and save gsea results into new variable for Pathway or Reactome analysis
DEG_comp <- res_O786EV_KO_D0
DEG_comp <- res_O786EV_KO_D5
DEG_comp <- res_O786EV_KO_D15
DEG_comp <- res_O786EV_KO_D23

DEG_comp <- res_O786KO_D0D5
DEG_comp <- res_O786KO_D0D15
DEG_comp <- res_O786KO_D0D23

DEG_comp <- res_O786EV_D0D5
DEG_comp <- res_O786EV_D0D15
DEG_comp <- res_O786EV_D0D23

DEG_comp <- res_A498_D0D5
DEG_comp <- res_A498_D0D15

DEG_comp <- res_CAKI_D0D12
DEG_comp <- res_CAKI_D0D25

DEG_comp <- res_ACHN_D0D5
DEG_comp <- res_ACHN_D0D14
DEG_comp <- res_ACHN_D0D22

####Run fGSEA for Hallmark enrichment####
df <- data.frame(DEG_comp)
df$row <- rownames(DEG_comp)
df$Gene.name<- apply(df,1,function(x) {unlist(strsplit(x[7],"\\|"))[2]})
df2 <- df %>% 
  dplyr::select(Gene.name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene.name) %>% dplyr::summarize(stat= mean(stat))

ranks <- deframe(df2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks) #For Hallmark analysis

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Take top 10 positive or negative NES after filtering for hallmark pathways that are padj <0.05 
dataplot <- fgseaResTidy[fgseaResTidy$padj <0.05,]
dataplot_pos <- dataplot[dataplot$NES > 0,]
dataplot_pos <- dataplot_pos[order(dataplot_pos$padj,-(dataplot_pos$NES)),]
dataplot_pos <- dataplot_pos[c(1:10),]
dataplot_neg <- dataplot[dataplot$NES < 0,]
dataplot_neg <- dataplot_neg[order(dataplot_neg$padj,dataplot_neg$NES),]
dataplot_neg <- dataplot_neg[c(1:10),]
dataplot_top10 <- rbind(dataplot_pos,dataplot_neg)
dataplot_top10 <-dataplot_top10 %>% filter(!is.na(pathway))

#store GSEA hallmark top 15 NES results
res_O786EV_KO_D0_HALLMARK <-data.frame(dataplot_top10)
res_O786EV_KO_D5_HALLMARK <-data.frame(dataplot_top10)
plotEnrichment(pathways.hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],stats=ranks) #plot for EV vs KO Day5 comparison (Fig S2C)
plotEnrichment(pathways.hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],stats=ranks) #plot for EV vs KO Day5 comparison (Fig S2C)
res_O786EV_KO_D15_HALLMARK <-data.frame(dataplot_top10)
res_O786EV_KO_D23_HALLMARK <-data.frame(dataplot_top10)

res_O786KO_D0D5_HALLMARK <-data.frame(dataplot_top10)
res_O786KO_D0D15_HALLMARK <-data.frame(dataplot_top10)
res_O786KO_D0D23_HALLMARK <-data.frame(dataplot_top10)

res_O786EV_D0D5_HALLMARK <-data.frame(dataplot_top10)
res_O786EV_D0D15_HALLMARK <-data.frame(dataplot_top10)
res_O786EV_D0D23_HALLMARK <-data.frame(dataplot_top10)

res_A498_D0D5_HALLMARK <-data.frame(dataplot_top10)
res_A498_D0D15_HALLMARK <-data.frame(dataplot_top10)

res_CAKI_D0D12_HALLMARK <-data.frame(dataplot_top10)
res_CAKI_D0D25_HALLMARK <-data.frame(dataplot_top10)

res_ACHN_D0D5_HALLMARK <-data.frame(dataplot_top10)
res_ACHN_D0D14_HALLMARK <-data.frame(dataplot_top10)
res_ACHN_D0D22_HALLMARK <-data.frame(dataplot_top10)


### #Run fGSEA for Reactome enrichment ###
df <- data.frame(DEG_comp)
df$row <- rownames(DEG_comp)
df$Gene.name<- apply(df,1,function(x) {unlist(strsplit(x[7],"\\|"))[2]})
df2 <- df %>% 
  dplyr::select(Gene.name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene.name) %>% dplyr::summarize(stat= mean(stat))

ranks <- deframe(df2)

fgseaRes <- fgseaMultilevel(pathways=pathways.reactome, stats=ranks) #For Reactome analysis

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Take top 10 positive or negative NES after filtering for hallmark pathways that are padj <0.01
dataplot <- fgseaResTidy[fgseaResTidy$padj <0.01,]
dataplot_pos <- dataplot[dataplot$NES > 0,]
dataplot_pos <- dataplot_pos[order(dataplot_pos$padj,-(dataplot_pos$NES)),]
dataplot_pos <- dataplot_pos[c(1:10),]
dataplot_neg <- dataplot[dataplot$NES < 0,]
dataplot_neg <- dataplot_neg[order(dataplot_neg$padj,dataplot_neg$NES),]
dataplot_neg <- dataplot_neg[c(1:10),]
dataplot_top10 <- rbind(dataplot_pos,dataplot_neg)
dataplot_top10 <-dataplot_top10 %>% filter(!is.na(pathway))

#store GSEA Reactome top 5 NES results
res_O786EV_KO_D0_Reactome <-data.frame(dataplot_top10)
plotEnrichment(pathways.reactome[["REACTOME_MRNA_SPLICING"]],stats=ranks) ##(Fig S4B)
res_O786EV_KO_D5_Reactome <-data.frame(dataplot_top10)
res_O786EV_KO_D15_Reactome <-data.frame(dataplot_top10)
res_O786EV_KO_D23_Reactome <-data.frame(dataplot_top10)

res_O786KO_D0D5_Reactome <-data.frame(dataplot_top10)
res_O786KO_D0D15_Reactome <-data.frame(dataplot_top10)
res_O786KO_D0D23_Reactome <-data.frame(dataplot_top10)

res_O786EV_D0D5_Reactome <-data.frame(dataplot_top10)
res_O786EV_D0D15_Reactome <-data.frame(dataplot_top10)
res_O786EV_D0D23_Reactome <-data.frame(dataplot_top10)

res_A498_D0D5_Reactome <-data.frame(dataplot_top10)
res_A498_D0D15_Reactome <-data.frame(dataplot_top10)

res_CAKI_D0D12_Reactome <-data.frame(dataplot_top10)
res_CAKI_D0D25_Reactome <-data.frame(dataplot_top10)

res_ACHN_D0D5_Reactome <-data.frame(dataplot_top10)
res_ACHN_D0D14_Reactome <-data.frame(dataplot_top10)
res_ACHN_D0D22_Reactome <-data.frame(dataplot_top10)


##### Generate combined heatmap of GSEA pathway enrichment (Fig S2A) ####
#786-O EV vs KO GSEA Hallmark Enrichment 
df.R<-data.frame(unique(c(res_O786EV_KO_D0_HALLMARK$pathway,
                          res_O786EV_KO_D5_HALLMARK$pathway,
                          res_O786EV_KO_D15_HALLMARK$pathway,
                          res_O786EV_KO_D23_HALLMARK$pathway)))
colnames(df.R) <-"pathway"
df.R <-na.omit(df.R)
df.R<- merge(df.R,res_O786EV_KO_D0_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_KO_D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_KO_D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_KO_D23_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","EV_KO_D0","EV_KO_D5","EV_KO_D15","EV_KO_D23")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:5])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:5])
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)

#786-O EV across timepoint comparison & 786-O KO across timepoint, GSEA Hallmark Enrichment (Fig 3B)
df.R<-data.frame(unique(c(res_O786EV_D0D5_HALLMARK$pathway,
                          res_O786EV_D0D15_HALLMARK$pathway,
                          res_O786EV_D0D23_HALLMARK$pathway,
                          res_O786KO_D0D5_HALLMARK$pathway,
                          res_O786KO_D0D15_HALLMARK$pathway,
                          res_O786KO_D0D23_HALLMARK$pathway)))
df.R <-na.omit(df.R)
colnames(df.R) <-"pathway"
df.R<- merge(df.R,res_O786EV_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_D0D23_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D23_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","EV_D0_D5","EV_D0_D15","EV_D0_D23","KO_D0_D5","KO_D0_D15","KO_D0_D23")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:7])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:7])
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)

#A498 GSEA Hallmark Enrichment across timepoint comparison (Fig S2A)
df.R<-data.frame(unique(c(res_A498_D0D5_HALLMARK[order(res_A498_D0D5_HALLMARK$NES),]$pathway,
                          res_A498_D0D15_HALLMARK[order(res_A498_D0D15_HALLMARK$NES),]$pathway)))
colnames(df.R) <-"pathway"
df.R <-na.omit(df.R)
df.R<- merge(df.R,res_A498_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_A498_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","A498_D0_D5","A498_D0_D15")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:3])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:3])
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)

#Caki-1 GSEA Hallmark Enrichment across timepoint comparison (Fig S2A)
df.R<-data.frame(unique(c(res_CAKI_D0D12_HALLMARK[order(res_CAKI_D0D12_HALLMARK$NES),]$pathway,
                          res_CAKI_D0D25_HALLMARK[order(res_CAKI_D0D25_HALLMARK$NES),]$pathway)))
colnames(df.R) <-"pathway"
df.R <-na.omit(df.R)
df.R<- merge(df.R,res_CAKI_D0D12_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_CAKI_D0D25_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","CAKI-1_D0_D12","CAKI-1_D0_D25")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:3])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:3])
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)

#ACHN GSEA Hallmark Enrichment across timepoint comparison (Fig S2A)
df.R<-data.frame(unique(c(res_ACHN_D0D5_HALLMARK[order(res_ACHN_D0D5_HALLMARK$NES),]$pathway,
                          res_ACHN_D0D14_HALLMARK[order(res_ACHN_D0D14_HALLMARK$NES),]$pathway,
                          res_ACHN_D0D22_HALLMARK[order(res_ACHN_D0D22_HALLMARK$NES),]$pathway)))
colnames(df.R) <-"pathway"
df.R <-na.omit(df.R)
df.R<- merge(df.R,res_ACHN_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ACHN_D0D14_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ACHN_D0D22_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","ACHN_D0_D5","ACHN_D0_D14","ACHN_D0_D22")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})  

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:4])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:4])
pheatmap(df.OG2,clustering_method = "single", color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)


#####786-O EV vs KO GSEA Reactome Enrichment (Figure S4A)
df.R<-data.frame(unique(c(res_O786EV_KO_D0_Reactome$pathway,
                          res_O786EV_KO_D5_Reactome$pathway,
                          res_O786EV_KO_D15_Reactome$pathway,
                          res_O786EV_KO_D23_Reactome$pathway)))
colnames(df.R) <-"pathway"
df.R <-na.omit(df.R)
df.R<- merge(df.R,res_O786EV_KO_D0_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_KO_D5_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_KO_D15_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_KO_D23_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","EV_KO_D0","EV_KO_D5","EV_KO_D15","EV_KO_D23")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"TOME_"))[2]})

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:5])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:5])
pheatmap(df.OG2,clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)

#786-O EV across timepoint comparison & 786-O KO across timepoint, GSEA Reactome Enrichment (Figure S4B)
df.R<-data.frame(unique(c(res_O786EV_D0D5_Reactome$pathway,
                          res_O786EV_D0D15_Reactome$pathway,
                          res_O786EV_D0D23_Reactome$pathway,
                          res_O786KO_D0D5_Reactome$pathway,
                          res_O786KO_D0D15_Reactome$pathway,
                          res_O786KO_D0D23_Reactome$pathway)))
df.R <-na.omit(df.R)
colnames(df.R) <-"pathway"
df.R<- merge(df.R,res_O786EV_D0D5_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_D0D15_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_D0D23_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D5_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D15_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D23_Reactome[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","EV_D0_D5","EV_D0_D15","EV_D0_D23","KO_D0_D5","KO_D0_D15","KO_D0_D23")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"TOME_"))[2]})

name <- df.R$pathway_simple
df.OG2 <- data.matrix(df.R[,2:7])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:7])
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)



###########Save normalized count values for Zscore heatmaps ###########
cdds<-data.frame(counts(dds, normalized = TRUE))
cdds$GeneID <- rownames(cdds)
cdds$Gene.name <- apply(cdds, 1, function(x) {unlist(strsplit(x[37],"\\|"))[2]})
colnames(cdds)[37:38] <- c("Gene","GeneID")
cdds<-cdds[,c("Gene","GeneID",meta1$Sample)]



#################### Gene pathway Z-score heatmaps ##########################
##Viral mimicry gene Zscore heatmap (Fig 3A)
VMG <- c("DDX58","DHX58","IFI16","IFI27","IFI30","IFI6","IFIH1","IFIT1","IFIT2","IFITM1","IFITM3","IRF7","IRF9","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OASL","STAT1","IFI35","IFIT3","OAS3","NFKBIA","NMI","B2M","STING1")

SETD2_VM <-cdds[which(cdds$GeneID %in% VMG==T),]
SETD2_VM <-SETD2_VM[order(SETD2_VM$GeneID),]

SETD2_VM1<-log2(SETD2_VM[,c(3:38)]+0.01)#Zscore
SETD2_VM1_Zscore<- (SETD2_VM1-rowMeans(SETD2_VM1))/(rowSds(as.matrix(SETD2_VM1)))[row(SETD2_VM1)]
name <- SETD2_VM[,c("GeneID")]
df.OG2 <- data.matrix(SETD2_VM1_Zscore)
row.names(df.OG2) <- name
pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =T, cluster_cols=F)

##IFNy Gene expression (Fig S2B)
IFNy <- read.table("GSEA_Hallmark_interferon_gamma_response.txt", sep = "\t", header = F, stringsAsFactors = F)
colnames(IFNy) <-"GeneID"

SETD2_IFNy <-cdds[which(cdds$GeneID %in% IFNy$GeneID==T),]
SETD2_IFNy <-SETD2_IFNy[order(SETD2_IFNy$GeneID),]

SETD2_IFNy1<-log2(SETD2_IFNy[,c(3:38)]+0.01)#Zscore
SETD2_IFNy1_Zscore<- (SETD2_IFNy1-rowMeans(SETD2_IFNy1))/(rowSds(as.matrix(SETD2_IFNy1)))[row(SETD2_IFNy1)]
name <- SETD2_IFNy[,c("GeneID")]
df.OG2 <- data.matrix(SETD2_IFNy1_Zscore)
row.names(df.OG2) <- name
pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =T, cluster_cols=F)

##mRNA splicing gene Zscore heatmap, only 786-O (Fig S4D)
splicing <- read.table("../fGSEA/REACTOME_mRNA_Splicing.txt", sep = "\t")
SETD2_SPL <-cdds[which(cdds$GeneID %in% splicing$V1==T),]
SETD2_SPL <-SETD2_SPL[order(SETD2_SPL$GeneID),]

SETD2_SPL1<-log2(SETD2_SPL[,c(11:26)]+0.01)#Zscore
SETD2_SPL1_Zscore<- (SETD2_SPL1-rowMeans(SETD2_SPL1))/(rowSds(as.matrix(SETD2_SPL1)))[row(SETD2_SPL1)]
name <- SETD2_SPL[,c("GeneID")]
df.OG2 <- data.matrix(SETD2_SPL1_Zscore)
row.names(df.OG2) <- name
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =T, cluster_cols=F)




###################### Transposable Element Analysis ################################
setwd("../TEexpression")
TE_cpm <- read.table("SETD2_TE_expression_all_sample_CPM.bed", sep ="\t", header = T, stringsAsFactors = F)
TE_cpm <- TE_cpm[,c(colnames(TE_cpm[1:6]),meta1$Sample)] 

TE_dic <- read.table("GRCh38_GENCODE_rmsk_TE_dictionary.txt",sep ="\t", header = F, stringsAsFactors = F)
colnames(TE_dic) <- c("TE_Subfamily","Geneid","TE_Family","TE_Class")
TE_dic <- TE_dic[!duplicated(TE_dic),]

TE_cpm <- merge(TE_cpm, TE_dic, by = "Geneid", x.all =T)
TE_cpm <- TE_cpm[TE_cpm$TE_Class == "LTR" | TE_cpm$TE_Class == "LINE" | TE_cpm$TE_Class == "SINE" | TE_cpm$TE_Class == "DNA" | TE_cpm$TE_Class == "SINE" ,] #only LINE, SINE, LTR, DNA


TE_annot <- read.delim("GRCh38_GENCODE_rmsk_TE_GTF_parsed.V37annotate.bed", sep ="\t", header =T, stringsAsFactors = F)
TE_annot <- TE_annot[,c(2,3,4,8)]
TE_annot$Anno2 <- apply(TE_annot, 1, function(x) {unlist(strsplit(x[4], " \\("))[1]})
TE_annot$Start <- TE_annot$Start-1 #make into 0-based
TE_cpm_anno_all<- merge(TE_cpm,TE_annot, by=c("Chr","Start","End")) 
TE_cpm_anno_all <-TE_cpm_anno_all[,c(1:6,43,44,45,47,7:42)]


filter_1cpm <- which(apply(TE_cpm_anno_all[,11:46], 1, function (x) (max(x)>=1))) ## 1CPM filter
TE_cpm_anno <- TE_cpm_anno_all[filter_1cpm,] 
TE_cpm_anno <-TE_cpm_anno[TE_cpm_anno$Anno2 == "intron" | TE_cpm_anno$Anno2 == "Intergenic",] #only look at intergenic or intronic TEs
#Visualize with heatmaps to verify replicability of biological replicates before averaging biological replicates

###Average CPM for biological Reps####
TE_cpm_anno$ACHNWA0_ave <- rowMeans(TE_cpm_anno[,c(11:12)])
TE_cpm_anno$ACHNWA5_ave <- rowMeans(TE_cpm_anno[,c(13:14)])
TE_cpm_anno$ACHNWA14_ave <- rowMeans(TE_cpm_anno[,c(15:16)])
TE_cpm_anno$ACHNWA22_ave <- rowMeans(TE_cpm_anno[,c(17:18)])
TE_cpm_anno$OEVD0_ave <- rowMeans(TE_cpm_anno[,c(19:20)])
TE_cpm_anno$OEVD5_ave <- rowMeans(TE_cpm_anno[,c(21:22)])
TE_cpm_anno$OEVD15_ave <- rowMeans(TE_cpm_anno[,c(23:24)])
TE_cpm_anno$OEVD23_ave <- rowMeans(TE_cpm_anno[,c(25:26)])
TE_cpm_anno$SK0D0_ave <- rowMeans(TE_cpm_anno[,c(27:28)])
TE_cpm_anno$SK0D5_ave <- rowMeans(TE_cpm_anno[,c(29:30)])
TE_cpm_anno$SK0D15_ave <- rowMeans(TE_cpm_anno[,c(31:32)])
TE_cpm_anno$SK0D23_ave <- rowMeans(TE_cpm_anno[,c(33:34)])
TE_cpm_anno$CAKIWA0_ave <- rowMeans(TE_cpm_anno[,c(35:36)])
TE_cpm_anno$CAKIWA12_ave <- rowMeans(TE_cpm_anno[,c(37:38)])
TE_cpm_anno$CAKIWA25_ave  <- rowMeans(TE_cpm_anno[,c(39:40)])
TE_cpm_anno$A498MD0_ave<- rowMeans(TE_cpm_anno[,c(41:42)])
TE_cpm_anno$A498MA5_ave <- rowMeans(TE_cpm_anno[,c(43:44)])
TE_cpm_anno$A498MA15_ave <- rowMeans(TE_cpm_anno[,c(45:46)])

TE_cpm_anno_ave_ACHN <- TE_cpm_anno[,c(1:10,47:50)]
TE_cpm_anno_ave_O786 <- TE_cpm_anno[,c(1:10,51:58)]
TE_cpm_anno_ave_CAKI <- TE_cpm_anno[,c(1:10,59:61)]
TE_cpm_anno_ave_A498 <- TE_cpm_anno[,c(1:10,62:64)]

filter_1cpm <- which(apply(TE_cpm_anno_ave_ACHN[,11:14], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_ACHN <- TE_cpm_anno_ave_ACHN[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_O786[,11:18], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_O786<- TE_cpm_anno_ave_O786[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_CAKI[,11:13], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_CAKI <- TE_cpm_anno_ave_CAKI[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_A498[,11:13], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_A498 <- TE_cpm_anno_ave_A498[filter_1cpm,]

TE_cpm_anno_ave <- TE_cpm_anno[,c(1:10,47:64)]
filter_1cpm <- which(apply(TE_cpm_anno_ave[,11:28], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave <- TE_cpm_anno_ave[filter_1cpm,] 


##Calculate TE expression fold change
write.table(TE_cpm_anno_ave[,11:28]+0.01, "SETD2_TE_averageCPM_V37anno_Days.txt", sep = "\t", col.names = F, row.names = F, quote = F)
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1, $4/$1, $5/$5, $6/$5, $7/$5, $8/$5, $9/$9, $10/$9, $11/$9, $12/$9, $13/$13, $14/$13, $15/$13, $16/$16, $17/$16, $18/$16}' SETD2_TE_averageCPM_V37anno_Days.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange.txt

TE_cpm_anno_ave_FC <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC <- cbind(TE_cpm_anno_ave[,c(1:10)],TE_cpm_anno_ave_FC)
colnames(TE_cpm_anno_ave_FC) <- colnames(TE_cpm_anno_ave)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC[,11:28], 1, function (x) (max(x)>=2))) 

TE_cpm_anno_ave_FC2 <- TE_cpm_anno_ave_FC[filter_2FC,]


write.table(TE_cpm_anno_ave_ACHN[,11:14]+0.01, "SETD2_TE_averageCPM_V37anno_Days_ACHN.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_O786[,11:18]+0.01, "SETD2_TE_averageCPM_V37anno_Days_O786.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_CAKI[,11:13]+0.01, "SETD2_TE_averageCPM_V37anno_Days_CAKI.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_A498[,11:13]+0.01, "SETD2_TE_averageCPM_V37anno_Days_A498.txt", sep = "\t", col.names = F, row.names = F, quote = F)

#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1}' SETD2_TE_averageCPM_V37anno_Days_ACHN.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_ACHN.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1,$5/$5,$6/$5,$7/$5,$8/$5}' SETD2_TE_averageCPM_V37anno_Days_O786.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_O786.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1}' SETD2_TE_averageCPM_V37anno_Days_CAKI.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_CAKI.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1}' SETD2_TE_averageCPM_V37anno_Days_A498.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_A498.txt

TE_cpm_anno_ave_FC_ACHN <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_ACHN.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_ACHN <- cbind(TE_cpm_anno_ave_ACHN[,c(1:10)],TE_cpm_anno_ave_FC_ACHN)
colnames(TE_cpm_anno_ave_FC_ACHN) <- colnames(TE_cpm_anno_ave_ACHN)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_ACHN[,11:14], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_ACHN <- TE_cpm_anno_ave_FC_ACHN[filter_2FC,]

TE_cpm_anno_ave_FC2_both_ACHN<-TE_cpm_anno_ave_FC2_ACHN #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_ACHN<-TE_cpm_anno_ave_FC2_both_ACHN[order(TE_cpm_anno_ave_FC2_both_ACHN$TE_Class,TE_cpm_anno_ave_FC2_both_ACHN$TE_Family,TE_cpm_anno_ave_FC2_both_ACHN$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_ACHN <-TE_cpm_anno_ave_FC2_ACHN[TE_cpm_anno_ave_FC2_ACHN$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_ACHN<-TE_cpm_anno_ave_FC2_intron_ACHN[order(TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class,TE_cpm_anno_ave_FC2_intron_ACHN$TE_Family,TE_cpm_anno_ave_FC2_intron_ACHN$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_ACHN <-TE_cpm_anno_ave_FC2_ACHN[TE_cpm_anno_ave_FC2_ACHN$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_ACHN<-TE_cpm_anno_ave_FC2_intergenic_ACHN[order(TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class,TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Family,TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Subfamily),] #sort by TE name



TE_cpm_anno_ave_FC_O786 <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_O786.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_O786 <- cbind(TE_cpm_anno_ave_O786[,c(1:10)],TE_cpm_anno_ave_FC_O786)
colnames(TE_cpm_anno_ave_FC_O786) <- colnames(TE_cpm_anno_ave_O786)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_O786[,11:18], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_O786 <- TE_cpm_anno_ave_FC_O786[filter_2FC,]

TE_cpm_anno_ave_FC2_both_O786<-TE_cpm_anno_ave_FC2_O786 #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_O786<-TE_cpm_anno_ave_FC2_both_O786[order(TE_cpm_anno_ave_FC2_both_O786$TE_Class,TE_cpm_anno_ave_FC2_both_O786$TE_Family,TE_cpm_anno_ave_FC2_both_O786$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_O786 <-TE_cpm_anno_ave_FC2_O786[TE_cpm_anno_ave_FC2_O786$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_O786<-TE_cpm_anno_ave_FC2_intron_O786[order(TE_cpm_anno_ave_FC2_intron_O786$TE_Class,TE_cpm_anno_ave_FC2_intron_O786$TE_Family,TE_cpm_anno_ave_FC2_intron_O786$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_O786 <-TE_cpm_anno_ave_FC2_O786[TE_cpm_anno_ave_FC2_O786$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_O786<-TE_cpm_anno_ave_FC2_intergenic_O786[order(TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class,TE_cpm_anno_ave_FC2_intergenic_O786$TE_Family,TE_cpm_anno_ave_FC2_intergenic_O786$TE_Subfamily),] #sort by TE name



TE_cpm_anno_ave_FC_CAKI <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_CAKI.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_CAKI <- cbind(TE_cpm_anno_ave_CAKI[,c(1:10)],TE_cpm_anno_ave_FC_CAKI)
colnames(TE_cpm_anno_ave_FC_CAKI) <- colnames(TE_cpm_anno_ave_CAKI)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_CAKI[,11:13], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_CAKI <- TE_cpm_anno_ave_FC_CAKI[filter_2FC,]

TE_cpm_anno_ave_FC2_both_CAKI<-TE_cpm_anno_ave_FC2_CAKI #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_CAKI<-TE_cpm_anno_ave_FC2_both_CAKI[order(TE_cpm_anno_ave_FC2_both_CAKI$TE_Class,TE_cpm_anno_ave_FC2_both_CAKI$TE_Family,TE_cpm_anno_ave_FC2_both_CAKI$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_CAKI <-TE_cpm_anno_ave_FC2_CAKI[TE_cpm_anno_ave_FC2_CAKI$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_CAKI<-TE_cpm_anno_ave_FC2_intron_CAKI[order(TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class,TE_cpm_anno_ave_FC2_intron_CAKI$TE_Family,TE_cpm_anno_ave_FC2_intron_CAKI$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_CAKI <-TE_cpm_anno_ave_FC2_CAKI[TE_cpm_anno_ave_FC2_CAKI$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_CAKI<-TE_cpm_anno_ave_FC2_intergenic_CAKI[order(TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class,TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Family,TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Subfamily),] #sort by TE name


TE_cpm_anno_ave_FC_A498 <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_A498.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_A498 <- cbind(TE_cpm_anno_ave_A498[,c(1:10)],TE_cpm_anno_ave_FC_A498)
colnames(TE_cpm_anno_ave_FC_A498) <- colnames(TE_cpm_anno_ave_A498)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_A498[,11:13], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_A498 <- TE_cpm_anno_ave_FC_A498[filter_2FC,]

TE_cpm_anno_ave_FC2_both_A498<-TE_cpm_anno_ave_FC2_A498 #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_A498<-TE_cpm_anno_ave_FC2_both_A498[order(TE_cpm_anno_ave_FC2_both_A498$TE_Class,TE_cpm_anno_ave_FC2_both_A498$TE_Family,TE_cpm_anno_ave_FC2_both_A498$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_A498 <-TE_cpm_anno_ave_FC2_A498[TE_cpm_anno_ave_FC2_A498$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_A498<-TE_cpm_anno_ave_FC2_intron_A498[order(TE_cpm_anno_ave_FC2_intron_A498$TE_Class,TE_cpm_anno_ave_FC2_intron_A498$TE_Family,TE_cpm_anno_ave_FC2_intron_A498$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_A498 <-TE_cpm_anno_ave_FC2_A498[TE_cpm_anno_ave_FC2_A498$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_A498<-TE_cpm_anno_ave_FC2_intergenic_A498[order(TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class,TE_cpm_anno_ave_FC2_intergenic_A498$TE_Family,TE_cpm_anno_ave_FC2_intergenic_A498$TE_Subfamily),] #sort by TE name



save.image("SETD2_TE_expression.RData")


#####Bar plot of number of TEs that are >2FC (Fig 3C)#######
bothTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22"),c(nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$SK0D23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A498[TE_cpm_anno_ave_FC2_both_A498$A498MA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A498[TE_cpm_anno_ave_FC2_both_A498$A498MA15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI[TE_cpm_anno_ave_FC2_both_CAKI$CAKIWA12_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI[TE_cpm_anno_ave_FC2_both_CAKI$CAKIWA25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_ACHN[TE_cpm_anno_ave_FC2_both_ACHN$ACHNWA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_ACHN[TE_cpm_anno_ave_FC2_both_ACHN$ACHNWA14_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_ACHN[TE_cpm_anno_ave_FC2_both_ACHN$ACHNWA22_ave >=2,])))
colnames(bothTE_FC2) <-c("Sample","Count")
bothTE_FC2<-merge(meta[,c(1:6)],bothTE_FC2, by="Sample")
bothTE_FC2$Day <-factor(bothTE_FC2$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2$Cell <-factor(bothTE_FC2$Cell, levels = c("ACHN","O786","CAKI","A498"))
bothTE_FC2$CellType <- paste(bothTE_FC2$Cell,bothTE_FC2$Type,sep="_")

p <- ggplot(bothTE_FC2, aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic+intron TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+ scale_y_continuous(limits = c(0, max(bothTE_FC2$Count)))+scale_fill_manual(values = mypalette[c(1,5,7,9,10)])+scale_color_manual(values= mypalette[c(1,5,7,9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p+ facet_wrap( ~ Cell, scales="free",ncol=4)+coord_cartesian(ylim = c(500,2000))


#Intergenic TE Bar plot of Number of TEs >2FC (Fig. 5A)
intergenicTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22"),c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA12_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA14_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA22_ave >=2,])))
colnames(intergenicTE_FC2) <-c("Sample","Count")
intergenicTE_FC2<-merge(meta[,c(1:6)],intergenicTE_FC2, by="Sample")
intergenicTE_FC2$Day <-factor(intergenicTE_FC2$Day, levels = c("5","12","14","15","22","23","25"))
intergenicTE_FC2$Cell <-factor(intergenicTE_FC2$Cell, levels = c("ACHN","O786","CAKI","A498"))
intergenicTE_FC2$CellType <- paste(intergenicTE_FC2$Cell,intergenicTE_FC2$Type,sep="_")

p <- ggplot(intergenicTE_FC2, aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+ scale_y_continuous(limits = c(0, max(intergenicTE_FC2$Count)))+scale_fill_manual(values = mypalette[c(1,5,7,9,10)])+scale_color_manual(values= mypalette[c(1,5,7,9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p1<-p+ facet_wrap( ~ Cell, scales="free",ncol=4)+coord_cartesian(ylim = c(0,1500))


#Intron TE Bar plot of Number of TEs >2FC (Fig. 5A)
intronTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22"),c(nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA12_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA14_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA22_ave >=2,])))
colnames(intronTE_FC2) <-c("Sample","Count")
intronTE_FC2<-merge(meta[,c(1:6)],intronTE_FC2, by="Sample")
intronTE_FC2$Day <-factor(intronTE_FC2$Day, levels = c("5","12","14","15","22","23","25"))
intronTE_FC2$Cell <-factor(intronTE_FC2$Cell, levels = c("ACHN","O786","CAKI","A498"))
intronTE_FC2$CellType <- paste(intronTE_FC2$Cell,intronTE_FC2$Type,sep="_")

p <- ggplot(intronTE_FC2, aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intron TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+ scale_y_continuous(limits = c(0, max(intronTE_FC2$Count)))+scale_fill_manual(values = mypalette[c(1,5,7,9,10)])+scale_color_manual(values= mypalette[c(1,5,7,9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p2<-p+ facet_wrap( ~ Cell, scales="free",ncol=4)+coord_cartesian(ylim = c(0,1500))


###    Make Z-score of TEs >2FC, 786-O only (Fig 3D)   ###############
O786_aveCPM_with_FC <- merge(TE_cpm_anno_ave_O786,TE_cpm_anno_ave_FC2_O786, by= c(colnames(TE_cpm_anno_ave)[1:10]))

O786_aveCPM_with_FC_Zscore<- (log2(O786_aveCPM_with_FC[,c(11:18)]+0.01)-rowMeans(log2(O786_aveCPM_with_FC[,c(11:18)]+0.01)))/(rowSds(as.matrix(log2(O786_aveCPM_with_FC[,c(11:18)]+0.01))))[row(log2(O786_aveCPM_with_FC[,c(11:18)]+0.01))]
name <- rownames(O786_aveCPM_with_FC_Zscore)
df.OG2 <- data.matrix(O786_aveCPM_with_FC_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(O786_aveCPM_with_FC_Zscore)
htname_786O =  Heatmap(df.OG2, column_title = "Timepoint",name= "O786 TEs\n(Zscore)",col = colorRamp2(c(-1.5,-0.8,0, 0.8, 1.5), c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000")), 
                       cluster_rows = T, cluster_columns = FALSE,show_row_names = F)
htname_786O




#######Identify Activated TEs (>1cpm) in O786 EV and KO COMPARISION###############
TE_cpm_anno_ave_O786 <- TE_cpm_anno_ave[,c(1:10,15:22)]
filter_1cpm <- which(apply(TE_cpm_anno_ave_O786[,11:18], 1, function (x) (max(x)>=1))) #some TEs are not 1cpm in 786-O, but might be included due to all cell line filter
TE_cpm_anno_ave_O786 <- TE_cpm_anno_ave_O786[filter_1cpm,] 

TE_cpm_anno_ave_O786_intron <-TE_cpm_anno_ave_O786[TE_cpm_anno_ave_O786$Anno2 == "intron",]
TE_cpm_anno_ave_O786_intron<-TE_cpm_anno_ave_O786_intron[order(TE_cpm_anno_ave_O786_intron$TE_Class,TE_cpm_anno_ave_O786_intron$TE_Family,TE_cpm_anno_ave_O786_intron$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_O786_intergenic <-TE_cpm_anno_ave_O786[TE_cpm_anno_ave_O786$Anno2 == "Intergenic",]
TE_cpm_anno_ave_O786_intergenic<-TE_cpm_anno_ave_O786_intergenic[order(TE_cpm_anno_ave_O786_intergenic$TE_Class,TE_cpm_anno_ave_O786_intergenic$TE_Family,TE_cpm_anno_ave_O786_intergenic$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_O786_both <-TE_cpm_anno_ave_O786[TE_cpm_anno_ave_O786$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786$Anno2 == "intron",]
TE_cpm_anno_ave_O786_both<-TE_cpm_anno_ave_O786_both[order(TE_cpm_anno_ave_O786_both$TE_Class,TE_cpm_anno_ave_O786_both$TE_Family,TE_cpm_anno_ave_O786_both$TE_Subfamily),] #sort by TE name


##identify TE >2FC
write.table(TE_cpm_anno_ave_O786[,11:18]+0.001, "SETD2_TE_averageTPM_Only786.txt", sep = "\t", col.names = F, row.names = F, quote = F) ### MAKE SURE TO REMAKE ONLY USING 786 CPM CUTOFF >1
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1, $4/$1, $5/$5, $6/$5, $7/$5, $8/$5}' SETD2_TE_averageTPM_Only786.txt > SETD2_TE_averageTPM_Only786_FoldChange.txt
#make fold change comparing EV to KO: awk -v OFS="\t" '{print $5/$1, $6/$2, $7/$3, $8/$4}' SETD2_TE_averageTPM_Only786.txt > SETD2_TE_averageTPM_Only786_FoldChange_CompareEVvsKO.txt
TE_cpm_anno_ave_O786_FC_EVvKO <- read.table("SETD2_TE_averageTPM_Only786_FoldChange_CompareEVvsKO.txt", sep = "\t", header =F, stringsAsFactors = F)

TE_cpm_anno_ave_O786_FC <- read.table("SETD2_TE_averageTPM_Only786_FoldChange.txt", sep = "\t", header =F, stringsAsFactors = F)
TE_cpm_anno_ave_O786_FC <- cbind(TE_cpm_anno_ave_O786[,1:10],TE_cpm_anno_ave_O786_FC,TE_cpm_anno_ave_O786_FC_EVvKO )
colnames(TE_cpm_anno_ave_O786_FC) <-c(colnames(TE_cpm_anno_ave_O786[,1:10]),paste0(colnames(TE_cpm_anno_ave_O786[,11:18]),"_FC"),"EVvKO_D0","EVvKO_D5","EVvKO_D15","EVvKO_D23")
filter_2FC <- which(apply(TE_cpm_anno_ave_O786_FC[,11:18], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_O786_FC2 <- TE_cpm_anno_ave_O786_FC[filter_2FC,] 

TE_cpm_anno_ave_O786_FC2_intron <-TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]
TE_cpm_anno_ave_O786_FC2_intron<-TE_cpm_anno_ave_O786_FC2_intron[order(TE_cpm_anno_ave_O786_FC2_intron$TE_Class,TE_cpm_anno_ave_O786_FC2_intron$TE_Family,TE_cpm_anno_ave_O786_FC2_intron$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_O786_FC2_intergenic <-TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]
TE_cpm_anno_ave_O786_FC2_intergenic<-TE_cpm_anno_ave_O786_FC2_intergenic[order(TE_cpm_anno_ave_O786_FC2_intergenic$TE_Class,TE_cpm_anno_ave_O786_FC2_intergenic$TE_Family,TE_cpm_anno_ave_O786_FC2_intergenic$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_O786_alldata_FC <- merge(TE_cpm_anno_ave_O786,TE_cpm_anno_ave_O786_FC[c(4,11:22)],by="Geneid")
filter_2FC <- which(apply(TE_cpm_anno_ave_O786_alldata_FC[,19:26], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_O786_alldata_FC2<- TE_cpm_anno_ave_O786_alldata_FC[filter_2FC,]
TE_cpm_anno_ave_O786_alldata_FC2_intron <- TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "intron",]

TE_cpm_anno_ave_O786_alldata_FC2_allKO <- TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC > 2,] 
TE_cpm_anno_ave_O786_alldata_FC2_allKO2 <- TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >1,] #686 intronic & 284 intergenic
TE_cpm_anno_ave_O786_alldata_FC2_allKOup <- TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 > 2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >1,] #344 intronic & 207 intergenic
TE_cpm_anno_ave_O786_alldata_FC2_allKOonly <- TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 > 2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2$OEVD5_ave <1 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >1,] #298 intronic & 184 intergenic


## Make plot comparing the CPM expression of induced INTRON TEs (>2FC) in EV vs KO (Fig. 5B) ####
TE_cpm_anno_ave_O786_alldata_FC2_intron_allTE <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$OEVD5_ave_FC > 2 | TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2,] 
TE_cpm_anno_ave_O786_alldata_FC2_intron_allKO <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2,] 
TE_cpm_anno_ave_O786_alldata_FC2_intron_allKO2 <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave >1,] #686 intronic 
TE_cpm_anno_ave_O786_alldata_FC2_intron_allKOup <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$EVvKO_D5 > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave >1,] #344 intronic 
TE_cpm_anno_ave_O786_alldata_FC2_intron_allKOonly <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$EVvKO_D5 > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$OEVD5_ave <1 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave >1,] #298 intronic

p<-ggplot(TE_cpm_anno_ave_O786_alldata_FC2_intron_allTE, aes(x=log2(OEVD5_ave) , y=log2(SK0D5_ave)) ) +
  geom_point(color="grey") + 
  theme_bw()
p<-p + geom_point(data = TE_cpm_anno_ave_O786_alldata_FC2_intron_allKOup, aes(x=log2(OEVD5_ave) , y=log2(SK0D5_ave)), color = "#ffcc00")	
p<-p + geom_point(data = TE_cpm_anno_ave_O786_alldata_FC2_intron_allKOonly, aes(x=log2(OEVD5_ave) , y=log2(SK0D5_ave)), color = "#cc3300")
p+geom_abline(slope=1, intercept=0)+ggtitle("EV intron TE vs KO intron TE (KO D5 > 2 FC)")+coord_cartesian(xlim=c(-6,6),ylim=c(-6,6))



#########Perform De Novo transcript assembly using Stringtie and identify which intronic TEs (induced >2FC) overlap with denovo transcript isoforms#######
### SPLICING DENOVO TRANSCRIPTS detection vis Stringtie #####
stringtie "${line}.bam" -p 16 -l ${line} -o ./stringtie_V37/GTF_denovo/"${line}_denovoAssembly2.gtf" --fr
stringtie --merge -c 2 -f 0.05 -i -p 16 --rf -o O786_stringtie_GTF_denovo_merged.gtf -G /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/genomes/hg38/gencode.v37.annotation.gtf O786_assembly_GTF_list.txt
gffcompare -r /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/genomes/hg38/gencode.v37.annotation.gtf -o /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/totalRNAseq/aligned/HISAT2_stranded/stringtie_V37/GTF_denovo/O786_gffcompare_V37_denovo O786_stringtie_GTF_denovo_merged.gtf

#Overlap induced intronic TEs with de novo transcript generated by Stringtie
setwd("../DeNovoSplicingIsoform")
tmap_denovo <- read.delim("gffcompare_V37_denovo.O786_stringtie_GTF_denovo_merged.gtf.tmap", sep = "\t", header =T, stringsAsFactors = F) #Categorize type of splicing using gffcompare
TE_denovo <- read.delim("TE_cpm_anno_ave_O786_FC2_Stringtie_Denovo_transcript_overlap_V37.bed", sep = "\t", header =F, stringsAsFactors = F)
TE_denovo <- TE_denovo[,c(1:10,31:40)]
colnames(TE_denovo) <-c("Chr","Start","End","Strand","TE_Subfamily","TE_Family","TE_Class","Geneid","Length","Anno2","chr","start","end","tool","anno","V6","strand","V8","GeneInfo","overlap")

TE_denovo2 <- merge(TE_denovo,TE_cpm_anno_ave_O786_alldata_FC2[,c(1,11:30)])
TE_denovo2 <- TE_denovo2[,c(1:10,21:40,11:20)]

TE_denovo_allKO <- TE_denovo2[TE_denovo2$SK0D5_ave_FC >= 2 & TE_denovo2$SK0D5_ave >=1,]
TE_denovo_allKO <- TE_denovo2[TE_denovo2$SK0D5_ave_FC >= 2 & TE_denovo2$EVvKO_D5 >= 2 & TE_denovo2$SK0D5_ave >=1,]
TE_denovo_allKO <- TE_denovo2[TE_denovo2$SK0D5_ave_FC >= 2 & TE_denovo2$EVvKO_D5 >= 2 & TE_denovo2$SK0D5_ave >=1 & TE_denovo2$OEVD5_ave <1,]

#iterate with different filters to get quantify what % of TEs overlap denovo or canonical exons transcript (Fig. S5A)
length(unique(TE_denovo_allKO[TE_denovo_allKO$Anno2 == "Intergenic",1])) #284/207/184
length(unique(TE_denovo_allKO[TE_denovo_allKO$Anno2 == "intron",1])) #683/344/298

TE_denovo_allKO_exon <- TE_denovo_allKO[(TE_denovo_allKO$Anno2 == "intron" | TE_denovo_allKO$Anno2 == "Intergenic") & TE_denovo_allKO$anno == "exon",]
TE_denovo_allKO_exon$qry_id <- apply(TE_denovo_allKO_exon ,1,function(x) {gsub(";","",unlist(strsplit(x[39], " "))[2])})
TE_denovo_allKO_exon$exon <- apply(TE_denovo_allKO_exon ,1,function(x) {gsub(";","",unlist(strsplit(x[39], " "))[6])})
TE_denovo_allKO_exon <- TE_denovo_allKO_exon[,c(1:34,37,41,42)]

TE_denovo_allKO_exon <- merge(TE_denovo_allKO_exon,tmap_denovo[,c("class_code","qry_gene_id","qry_id","num_exons")], by = "qry_id")
TE_denovo_allKO_exon <-TE_denovo_allKO_exon[!duplicated(TE_denovo_allKO_exon),]
TE_denovo_allKO_exon <- TE_denovo_allKO_exon[which(TE_denovo_allKO_exon$Strand == TE_denovo_allKO_exon$strand),] 
length(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "Intergenic",]$Geneid)) #199/143/121 or 67.7%/69.1%/65.8%
length(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron",]$Geneid)) #347/196/162 or 50.8%/57.0%/54.4%

t<- TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron",]
t$ExonLength <- t$end-t$start
t_filter_longest_exon <- unique(t[,c("Geneid","chr","start","end","ExonLength")]) #since lot of denovo exons overlap, annotate each TE to the longest exon detected

require(data.table) 
t_filter_longest_exon <- as.data.table(t_filter_longest_exon)
t_filter_longest_exon<-data.frame(t_filter_longest_exon[t_filter_longest_exon[, .I[which.max(ExonLength)], by=Geneid]$V1])##Filter longest exon containing TE
t_filter_longest_exon$ExonName <- paste(t_filter_longest_exon$chr,t_filter_longest_exon$start,t_filter_longest_exon$end,sep="_")
t_filter_longest_exon_location <- t_filter_longest_exon[,c("chr","start","end")]
t_filter_longest_exon_location <- t_filter_longest_exon_location[!duplicated(t_filter_longest_exon_location),]
tcount<-data.frame(table(t_filter_longest_exon$ExonName))

t_filter_longest_exon_with_info <- merge(t_filter_longest_exon[,c(1:6)],TE_denovo_allKO_exon,by = c("Geneid","chr","start","end"))
t_filter_longest_exon_with_info_only_exons <- unique(t_filter_longest_exon_with_info[,c("chr","start","end","ExonName","class_code")])
table(t_filter_longest_exon_with_info_only_exons$class_code) #Type of exon detected with induced intron TEs
##gffcompare exon categories
= : Exact match
c : containted in reference(intron compatible)
k : containment of reference (reverse containment)
m : retained intron(s), all introns matched or retained
n : retained introns, not all introns matched/covered
j : multi-exon with at least one junction match
e : single exon transfrag partially covering an intron, possible pre-mRNA fragment
o : other same strand overlap with reference exons
s : intron match on the opposite strand (likely a mapping error)
x : exonic overlap on the opposite strand ( like o or e but on the opposite strand) 
i : fully contained within a reference intron
y : contains a reference withint its intron(s)
p : possible polymerase run-on (no actual overlap)
r : repeat (at least 50% bases soft-masked)
u : none of the above (unknown, intergenic)
#

t_filter_longest_exon_with_info_only_exons2 <- t_filter_longest_exon_with_info_only_exons %>%
  dplyr::group_by(ExonName) %>%
  dplyr::summarise(class_code = paste(class_code, collapse = ","))
ttt1<-data.frame(t_filter_longest_exon_with_info_only_exons2)
table(t_filter_longest_exon_with_info_only_exons2$class_code)

=       =,j     =,j,k     =,k,j         i         j     j,=,n       j,i       j,k     j,k,= j,n,k,=,m         k       k,=     k,=,j       k,j     k,j,i         m       m,j       m,k         n       n,k         o 
26         1         2         2        40        51         1         1         6         1         1        40         2         2         4         1         9         1         1        10         1         3 
x         y 
16         1 

##Combine transcript class - Group 0 (=), Group 1 exon discrepancy (j,k,o), intron retention Group 2 (m,n), Group 3 novel transcript (i,x,y). Priority: Group1 > Group 2 >Group 3 > Group 0 
group 0: 26
group 1: 1+2+2+51+1+1+6+1+1+40+2+2+4+1+1+1+1+3 = 121
group 2: 9+10 = 19
group 3: 40 + 16+1 = 57

#make figure for frequency of different groups of exons (Fig 5C)
denovo <- data.frame(c(26,121,19,57),c(rep("allKO",4)),c(rep(c("g0","g1","g2","g3"),1)))
colnames(denovo) <- c("Count","Type","Group")
figurecolor <- c("#44bec7","#fa3c4c","#ffc300","#d696bb")

p <- ggplot(denovo, aes(x=Type, y=Count, fill=Group)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("De novo Exon type")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Type", y = "TE count")+scale_fill_manual(values = mypalette[c(1,5,7,9,10)])+scale_color_manual(values= figurecolor)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))







############################ HISTONE CHIP-seq ANALYSIS ############################
#R version 4.0.2 (2020-06-22)
library(ggplot2) #3.3.5
library(RColorBrewer) #1.1-2
library(plyr) #1.8.6
library("pheatmap") #1.0.12
library("ComplexHeatmap") #2.6.2
library("circlize") #0.4.13
library(reshape2) #1.4.4
library(matrixStats) #0.60.1
library(reshape2) #1.4.4
library(dplyr) #1.0.7
library(DiffBind) #3.0.15
library(data.table) #1.14.0
mypalette <- brewer.pal(12,"Paired")
mypalette2 <- brewer.pal(12,"Set3")
mypalette5 <- brewer.pal(9,"Pastel1")
mypalette6 <- brewer.pal(9,"Purples")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

####Diffbind
setwd("../ChIPseq/Diffbind")
samples<-read.csv("Diffbind_Sample_info_IDRpeaks.csv",stringsAsFactors =F)
samples_K4me3 <- samples[samples$Factor == "K4me3",] #narrow
samples_K36me3 <- samples[samples$Factor == "K36me3",] #broad
samples_K9me3 <- samples[samples$Factor == "K9me3",] #narrow
samples_K27me3 <- samples[samples$Factor == "K27me3",] #broad
samples_K27ac <- samples[samples$Factor == "K27ac",] #narrow

dba.samples_K4me3 <- dba(sampleSheet=samples_K4me3) %>%
  dba.blacklist(blacklist=DBA_BLACKLIST_HG38, greylist=F) %>%
  dba.count(summits=F) %>%
  dba.normalize(normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL,background=F) %>%
  dba.contrast(categories=DBA_TISSUE,minMembers = 2) %>%
  dba.analyze(bGreylist = F)

dba.samples_K36me3 <- dba(sampleSheet=samples_K36me3) %>%
  dba.blacklist(blacklist=DBA_BLACKLIST_HG38, greylist=F) %>%
  dba.count(summits=F) %>%
  dba.normalize(normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL,background=F) %>%
  dba.contrast(categories=DBA_TISSUE,minMembers = 2) %>%
  dba.analyze(bGreylist = F)

dba.samples_K9me3 <- dba(sampleSheet=samples_K9me3) %>%
  dba.blacklist(blacklist=DBA_BLACKLIST_HG38, greylist=F) %>%
  dba.count(summits=F) %>%
  dba.normalize(normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL,background=F) %>%
  dba.contrast(categories=DBA_TISSUE,minMembers = 2) %>%
  dba.analyze(bGreylist = F)

dba.samples_K27me3 <- dba(sampleSheet=samples_K27me3) %>%
  dba.blacklist(blacklist=DBA_BLACKLIST_HG38, greylist=F) %>%
  dba.count(summits=F) %>%
  dba.normalize(normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL,background=F) %>%
  dba.contrast(categories=DBA_TISSUE,minMembers = 2) %>%
  dba.analyze(bGreylist = F)

dba.samples_K27ac <- dba(sampleSheet=samples_K27ac) %>%
  dba.blacklist(blacklist=DBA_BLACKLIST_HG38, greylist=F) %>%
  dba.count(summits=F) %>%
  dba.normalize(normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL,background=F) %>%
  dba.contrast(categories=DBA_TISSUE,minMembers = 2) %>%
  dba.analyze(bGreylist = F)


##plot FRiP (Fig. S3B)
info <- dba.show(dba.samples_K36me3)
info <-rbind(info,dba.show(dba.samples_K9me3))
info <-rbind(info,dba.show(dba.samples_K4me3))
info <-rbind(info,dba.show(dba.samples_K27me3))
info <-rbind(info,dba.show(dba.samples_K27ac))
info$CellDay <-paste0(info$Condition,info$Treatment)
info$Factor <- factor(info$Factor ,levels = c("K4me3","K27ac","K27me3","K9me3","K36me3"))
info$CellDay <-factor(info$CellDay,levels=c("EV0","EV5","EV23","KO0","KO5","KO23"))

mypalette3 <- c("#2A4E69", "#4C86B5","#ADCBE2","#A91E22","#EF3B43","#F1BDBC")

pp <- ggplot(info, aes(x=Factor, y=FRiP, fill=CellDay)) +geom_bar(stat="identity",colour = "black", position=position_dodge(),alpha=0.6)+ggtitle("Frequency of Reads in Peak (FRiP)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Histone Marks", y = "FRiP")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp



##plot PCA from Diffbind (Fig. S3C)
dba.plotPCA(dba.samples_K4me3,score = "DBA_SCORE_RPKM_MINUS")
dba.plotPCA(dba.samples_K36me3,score = "DBA_SCORE_RPKM_MINUS")
dba.plotPCA(dba.samples_K27me3,score = "DBA_SCORE_RPKM_MINUS")
dba.plotPCA(dba.samples_K9me3,score = "DBA_SCORE_RPKM_MINUS")
dba.plotPCA(dba.samples_K27ac,score = "DBA_SCORE_RPKM_MINUS")



#### RPKM of each compiled peak by diff bind
DBA <-  dba.count(dba.samples_K36me3,peaks=NULL,score=DBA_SCORE_RPKM_MINUS)
K36me3_rpkm <- data.frame(dba.peakset(DBA, bRetrieve=TRUE))
DBA <-  dba.count(dba.samples_K27me3,peaks=NULL,score=DBA_SCORE_RPKM_MINUS)
K27me3_rpkm <- data.frame(dba.peakset(DBA, bRetrieve=TRUE))
DBA <-  dba.count(dba.samples_K9me3,peaks=NULL,score=DBA_SCORE_RPKM_MINUS)
K9me3_rpkm <- data.frame(dba.peakset(DBA, bRetrieve=TRUE))
DBA <-  dba.count(dba.samples_K4me3,peaks=NULL,score=DBA_SCORE_RPKM_MINUS)
K4me3_rpkm <- data.frame(dba.peakset(DBA, bRetrieve=TRUE))
DBA <-  dba.count(dba.samples_K27ac,peaks=NULL,score=DBA_SCORE_RPKM_MINUS)
K27ac_rpkm <- data.frame(dba.peakset(DBA, bRetrieve=TRUE))

write.table(K36me3_rpkm, "H3K36me3_Peak_RPKM_allsample.bed", sep ="\t", col.names = F, row.names = F, quote=F)
write.table(K27me3_rpkm, "H3K27me3_Peak_RPKM_allsample.bed", sep ="\t", col.names = F, row.names = F, quote=F)
write.table(K4me3_rpkm, "H3K4me3_Peak_RPKM_allsample.bed", sep ="\t", col.names = F, row.names = F, quote=F)
write.table(K9me3_rpkm, "H3K9me3_Peak_RPKM_allsample.bed", sep ="\t", col.names = F, row.names = F, quote=F)
write.table(K27ac_rpkm, "H3K27ac_Peak_RPKM_allsample.bed", sep ="\t", col.names = F, row.names = F, quote=F)

#Annotate all peaks using Homer (Fig S3D)
for i in *_Peak_RPKM_allsample.bed; do annotatePeaks.pl $i hg38 -genomeOntology ./Homer_annotate_V37/${i/.bed/}"_genomeOntology" -gtf /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/genomes/hg38/gencode.v37.annotation.gtf > ./Homer_annotate_V37/${i/.bed/}"_annotatePeaks_V37.bed";done
#
K36peak_anno <- read.delim("H3K36me3_Peak_RPKM_allsample_annotatePeaks_V37.bed", sep = "\t", header = T, stringsAsFactors = F)
K9peak_anno <- read.delim("H3K9me3_Peak_RPKM_allsample_annotatePeaks_V37.bed", sep = "\t", header = T, stringsAsFactors = F)
K4peak_anno <- read.delim("H3K4me3_Peak_RPKM_allsample_annotatePeaks_V37.bed", sep = "\t", header = T, stringsAsFactors = F)
K27peak_anno <- read.delim("H3K27me3_Peak_RPKM_allsample_annotatePeaks_V37.bed", sep = "\t", header = T, stringsAsFactors = F)
K27acpeak_anno <- read.delim("H3K27ac_Peak_RPKM_allsample_annotatePeaks_V37.bed", sep = "\t", header = T, stringsAsFactors = F)

K36peak_anno$Anno2 <-apply(K36peak_anno, 1, function(x) {unlist(strsplit(x[8], " \\("))[1]})
K9peak_anno$Anno2 <-apply(K9peak_anno, 1, function(x) {unlist(strsplit(x[8], " \\("))[1]})
K4peak_anno$Anno2 <-apply(K4peak_anno, 1, function(x) {unlist(strsplit(x[8], " \\("))[1]})
K27peak_anno$Anno2 <-apply(K27peak_anno, 1, function(x) {unlist(strsplit(x[8], " \\("))[1]})
K27acpeak_anno$Anno2 <-apply(K27acpeak_anno, 1, function(x) {unlist(strsplit(x[8], " \\("))[1]})

K36peak_anno$Start <- K36peak_anno$Start-1
K9peak_anno$Start <- K9peak_anno$Start-1
K4peak_anno$Start <- K4peak_anno$Start-1
K27peak_anno$Start <- K27peak_anno$Start-1
K27acpeak_anno$Start <- K27acpeak_anno$Start-1

colnames(K36me3_rpkm)[1:3] <-c("Chr","Start","End")
colnames(K9me3_rpkm)[1:3] <-c("Chr","Start","End")
colnames(K4me3_rpkm)[1:3] <-c("Chr","Start","End")
colnames(K27me3_rpkm)[1:3] <-c("Chr","Start","End")
colnames(K27ac_rpkm)[1:3] <-c("Chr","Start","End")

Peak_K36me3_anno <- merge(K36me3_rpkm,K36peak_anno[,c(2:4,20,19,8)],by =c("Chr","Start","End"))
Peak_K9me3_anno <- merge(K9me3_rpkm,K9peak_anno[,c(2:4,20,19,8)],by =c("Chr","Start","End"))
Peak_K4me3_anno <- merge(K4me3_rpkm,K4peak_anno[,c(2:4,20,19,8)],by =c("Chr","Start","End"))
Peak_K27me3_anno <- merge(K27me3_rpkm,K27peak_anno[,c(2:4,20,19,8)],by =c("Chr","Start","End"))
Peak_K27ac_anno <- merge(K27ac_rpkm,K27acpeak_anno[,c(2:4,20,19,8)],by =c("Chr","Start","End"))

Peak_K36me3_anno <-na.omit(Peak_K36me3_anno)
Peak_K9me3_anno <-na.omit(Peak_K9me3_anno)
Peak_K4me3_anno <-na.omit(Peak_K4me3_anno)
Peak_K27me3_anno <-na.omit(Peak_K27me3_anno)
Peak_K27ac_anno <-na.omit(Peak_K27ac_anno)

t_K36peak <- data.frame(table(Peak_K36me3_anno$Anno2))
t_K36peak$histone <- "K36me3"
t_K9peak <- data.frame(table(Peak_K9me3_anno$Anno2))
t_K9peak$histone <- "K9me3"
t_K4peak <- data.frame(table(Peak_K4me3_anno$Anno2))
t_K4peak$histone <- "K4me3"
t_K27peak <-data.frame( table(Peak_K27me3_anno$Anno2))
t_K27peak$histone <- "K27me3"
t_K27acpeak <- data.frame(table(Peak_K27ac_anno$Anno2))
t_K27acpeak$histone <- "K27ac"

t_allpeaks <-rbind(t_K36peak,t_K9peak,t_K4peak,t_K27peak,t_K27acpeak) 
t_allpeaks$histone <- factor(t_allpeaks$histone, levels = c("K4me3","K27ac","K27me3","K9me3","K36me3"))
t_allpeaks$Var1 <- factor(t_allpeaks$Var1, levels = c("Intergenic","promoter-TSS","5' UTR","exon","intron","3' UTR","TTS","non-coding"))

p2<-ggplot(t_allpeaks, aes(x=histone, y=Freq, fill=Var1)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("Consensus Histone Peaks Annotation v37\n(identified by Diffbind)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Histone", y = "Frequency")+scale_fill_manual(values = mypalette5)+scale_color_manual(values=mypalette5)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p2 


#### MAKE BIGWIG FROM BAM AND PLOT GENEBODY HISTONE LEVELS ####
cd ./ChIPseq/Rep1/aligned
for i in *.sorted.rmdup.bam; do bamCoverage --bam $i -o ./ChIPseq/bigwig_RPGC/${i/.sorted.rmdup.bam/_RPGC.bw} -p 16 --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX chrM chrY --extendReads;done
cd ./ChIPseq/Rep2/aligned
for i in *.sorted.rmdup.bam; do bamCoverage --bam $i -o ./ChIPseq/bigwig_RPGC/${i/.sorted.rmdup.bam/_RPGC.bw} -p 16 --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX chrM chrY --extendReads;done

#combine biological replicates for visualization
cd ./ChIPseq/bigwig_RPGC
for i in *_Rep1_RPGC.bw; do bigWigMerge $i ${i/Rep1/Rep2} ${i/_Rep1_RPGC.bw/_Combined.bedGraph};done
for i in *_Combined.bedGraph; do sort -k1,1 -k2,2n $i > ${i/.bed/.sorted.bed};done
for i in *_Combined.sorted.bedGraph; do bedGraphToBigWig $i /genomes/hg38/hg38.chrom.sizes ${i/.sorted.bedGraph/.bw};done

#deeptools heatmap and linegraph of histone signal across gene body (Fig 4B) . H3K9me3 example.
computeMatrix scale-regions --missingDataAsZero -bs 100 -p 12 -R /genomes/hg38/gencode.v37.AutosomalTranscriptsOnly.forDeepTools.gtf -S EV0K9_Combined_RPGC.bw EV5K9_Combined_RPGC.bw EV23K9_Combined_RPGC.bw SK0K9_Combined_RPGC.bw SK5K9_Combined_RPGC.bw SK23K9_Combined_RPGC.bw -m 5000 -b 3000 -a 3000 -out H3K9me3_GeneBody_O786_RPGC_combined.tab.gz 
plotHeatmap  --heatmapWidth 6 -m H3K9me3_GeneBody_O786_RPGC_combined.tab.gz --perGroup -out H3K9me3_GeneBody_O786_RPGC.pdf --missingDataColor "grey" --colorMap Greens --heatmapHeight 10




#### PERFORM TE centric analysis to quantify histone changes ###########
setwd("../TE_centric")
write.table(TE_cpm_anno_ave_O786_alldata_FC2[,c(2:4,1,5:30)],"SETD2_TE_averageTPM_Only786_FoldChangeGreaterThan2_V37.txt",sep = "\t", col.names = F, row.names =F, quote =F)
write.table(TE_cpm_anno_ave_O786_alldata_FC2[,c(2:4,1)],"TE_FC2_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F)

######### Quantify TE histone marks (RPGC) vis bedtools  #######
sort -k1,1 -k2,2n TE_FC2_O786_location.bed > TE_FC2_O786_location.sorted.bed
module load bbc/bedtools/bedtools-2.29.2 
for i in *K9_*Rep1.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *K4_*Rep1.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *K36_*Rep1.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *K27_*Rep1.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *Ac_*Rep1.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done

for i in *K9_*Rep2.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *K4_*Rep2.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *K36_*Rep2.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *K27_*Rep2.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done
for i in *Ac_*Rep2.sorted.bedgraph ; do bedtools map -a /varidata/research/home/josh.jang/jones-secondary/projects/JJANG/projects/SETD2_ccRCC/ChIPseq/Rcode/RPGC_bedgraph/bedtools_map2/TE_FC2_O786_location.sorted.bed  -b $i -c 4 -o mean > ./bedtools_map2/"TE_FC2_O786_location_"$i;done

bothTE_FC2_Histone <- read.table("TE_FC2_O786_location_Histone_BedtoolsMAP_RPGC.bed", sep = "\t", header = T, stringsAsFactors = F)
bothTE_FC2_Histone <-bothTE_FC2_Histone[,c("chr","start","end","EV0K4_Rep1","EV0K4_Rep2","EV5K4_Rep1","EV5K4_Rep2","EV23K4_Rep1","EV23K4_Rep2","SK0K4_Rep1","SK0K4_Rep2","SK5K4_Rep1","SK5K4_Rep2","SK23K4_Rep1","SK23K4_Rep2","EV0Ac_Rep1","EV0Ac_Rep2","EV5Ac_Rep1","EV5Ac_Rep2","EV23Ac_Rep1","EV23Ac_Rep2","SK0Ac_Rep1","SK0Ac_Rep2","SK5Ac_Rep1","SK5Ac_Rep2","SK23Ac_Rep1","SK23Ac_Rep2","EV0K27_Rep1","EV0K27_Rep2","EV5K27_Rep1","EV5K27_Rep2","EV23K27_Rep1","EV23K27_Rep2","SK0K27_Rep1","SK0K27_Rep2","SK5K27_Rep1","SK5K27_Rep2","SK23K27_Rep1","SK23K27_Rep2","EV0K9_Rep1","EV0K9_Rep2","EV5K9_Rep1","EV5K9_Rep2","EV23K9_Rep1","EV23K9_Rep2","SK0K9_Rep1","SK0K9_Rep2","SK5K9_Rep1","SK5K9_Rep2","SK23K9_Rep1","SK23K9_Rep2","EV0K36_Rep1","EV0K36_Rep2","EV5K36_Rep1","EV5K36_Rep2","EV23K36_Rep1","EV23K36_Rep2","SK0K36_Rep1","SK0K36_Rep2","SK5K36_Rep1","SK5K36_Rep2","SK23K36_Rep1","SK23K36_Rep2")]
colnames(bothTE_FC2_Histone)[1:3] <- c("Chr","Start","End")

TE_cpm_anno_aveCPM_O786_FC2_both_histone <- merge(TE_cpm_anno_ave_O786_alldata_FC2, bothTE_FC2_Histone, by = c("Chr","Start","End"))


TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K4_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K4_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K4_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K4_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K4_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K4_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K4_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K4_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K4_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K4_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K4_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K4_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K4_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K4_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K4_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K4_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K4_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K4_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0Ac_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0Ac_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0Ac_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5Ac_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5Ac_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5Ac_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23Ac_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23Ac_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23Ac_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0Ac_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0Ac_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0Ac_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5Ac_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5Ac_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5Ac_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23Ac_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23Ac_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23Ac_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K27_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K27_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K27_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K27_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K27_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K27_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K27_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K27_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K27_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K27_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K27_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K27_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K27_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K27_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K27_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K27_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K27_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K27_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K9_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K9_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K9_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K9_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K9_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K9_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K9_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K9_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K9_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K9_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K9_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K9_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K9_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K9_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K9_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K9_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K9_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K9_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K36_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K36_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV0K36_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K36_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K36_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV5K36_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K36_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K36_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$EV23K36_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K36_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K36_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK0K36_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K36_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K36_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK5K36_Rep2)/2
TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K36_ave <- (TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K36_Rep1+TE_cpm_anno_aveCPM_O786_FC2_both_histone$SK23K36_Rep2)/2

###all time points comparison but order on D5 TE expression
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(11:18)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2<- cbind(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(1:3)],TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5 <- TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$SK0D5_ave

TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE)]

#TE expression Zscore heatmap (Fig 4C)
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),])
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
colnames(df.OG2) <- colnames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),])
d<-pheatmap(df.OG2,color = colorRampPalette(c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000"))(100),breaks = c(seq(-1.5,1.5,by=0.03)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

#TE expression CPM boxplot (Fig 4C)
TE_aveCPM_m<-melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(11:18)])
colnames(TE_aveCPM_m) <- c("Sample","CPM")
TE_aveCPM_m$Cell <- apply(TE_aveCPM_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
TE_aveCPM_m$Sample <- factor(TE_aveCPM_m$Sample, levels = c("OEVD0_ave","OEVD5_ave","OEVD15_ave","OEVD23_ave","SK0D0_ave","SK0D5_ave","SK0D15_ave","SK0D23_ave"))
TE_aveCPM_m$Cell <- factor(TE_aveCPM_m$Cell, levels = c("OEVD0","OEVD5","OEVD15","OEVD23","SK0D0","SK0D5","SK0D15","SK0D23"))
ppppp<-ggplot(TE_aveCPM_m, aes(x=Sample, y=CPM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("TE Expression Distribution (CPM)")+ coord_cartesian(ylim=c(0, 3))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "CPM")+scale_fill_manual(values = c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+scale_color_manual(values=c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot



### histone Zscore heatmaps (Fig 4C) ###
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),c(91:96)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),c(97:102)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),c(103:108)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),c(109:114)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),c(115:120)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36)]

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-1.5,1.5,by=0.03)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-1.5,1.5,by=0.03)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-1.5,1.5,by=0.03)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-1.5,1.5,by=0.03)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-1.5,1.5,by=0.03)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d



##Boxplot
mypalette3 <- c("#2a4d69", "#4b86b4","#adcbe3","#9c0000","#fd0000","#ff9797")

K36me3_rpgc_m<-melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(115:120)])
colnames(K36me3_rpgc_m) <- c("Sample","RPGC")
K36me3_rpgc_m$Cell <- apply(K36me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K36me3_rpgc_m$Sample <- factor(K36me3_rpgc_m$Sample, levels = c("EV0K36_ave","EV5K36_ave","EV23K36_ave","SK0K36_ave","SK5K36_ave","SK23K36_ave"))
K36me3_rpgc_m$Cell <- factor(K36me3_rpgc_m$Cell, levels = c("EV0K36","EV5K36","EV23K36","SK0K36","SK5K36","SK23K36"))
ppppp<-ggplot(K36me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K36me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 26))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

K9me3_rpgc_m<-melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(109:114)])
colnames(K9me3_rpgc_m) <- c("Sample","RPGC")
K9me3_rpgc_m$Cell <- apply(K9me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K9me3_rpgc_m$Sample <- factor(K9me3_rpgc_m$Sample, levels = c("EV0K9_ave","EV5K9_ave","EV23K9_ave","SK0K9_ave","SK5K9_ave","SK23K9_ave"))
K9me3_rpgc_m$Cell <- factor(K9me3_rpgc_m$Cell, levels = c("EV0K9","EV5K9","EV23K9","SK0K9","SK5K9","SK23K9"))
pppp<-ggplot(K9me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K9me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pppp##boxplot


K27me3_rpgc_m<-melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(103:108)])
colnames(K27me3_rpgc_m) <- c("Sample","RPGC")
K27me3_rpgc_m$Cell <- apply(K27me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K27me3_rpgc_m$Sample <- factor(K27me3_rpgc_m$Sample, levels = c("EV0K27_ave","EV5K27_ave","EV23K27_ave","SK0K27_ave","SK5K27_ave","SK23K27_ave"))
K27me3_rpgc_m$Cell <- factor(K27me3_rpgc_m$Cell, levels = c("EV0K27","EV5K27","EV23K27","SK0K27","SK5K27","SK23K27"))
ppp<-ggplot(K27me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K27me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 4))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppp##boxplot

K4me3_rpgc_m<-melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(91:96)])
colnames(K4me3_rpgc_m) <- c("Sample","RPKM")
K4me3_rpgc_m$Cell <- apply(K4me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K4me3_rpgc_m$Sample <- factor(K4me3_rpgc_m$Sample, levels = c("EV0K4_ave","EV5K4_ave","EV23K4_ave","SK0K4_ave","SK5K4_ave","SK23K4_ave"))
K4me3_rpgc_m$Cell <- factor(K4me3_rpgc_m$Cell, levels = c("EV0K4","EV5K4","EV23K4","SK0K4","SK5K4","SK23K4"))
pp<-ggplot(K4me3_rpgc_m, aes(x=Sample, y=RPKM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K4me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 22))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPKM")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pp##boxplot


K27ac_rpgc_m<-melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(97:102)])
colnames(K27ac_rpgc_m) <- c("Sample","RPGC")
K27ac_rpgc_m$Cell <- apply(K27ac_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K27ac_rpgc_m$Sample <- factor(K27ac_rpgc_m$Sample, levels = c("EV0Ac_ave","EV5Ac_ave","EV23Ac_ave","SK0Ac_ave","SK5Ac_ave","SK23Ac_ave"))
K27ac_rpgc_m$Cell <- factor(K27ac_rpgc_m$Cell, levels = c("EV0Ac","EV5Ac","EV23Ac","SK0Ac","SK5Ac","SK23Ac"))
p<-ggplot(K27ac_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K27ac Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 12))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
p##boxplot

multiplot(ppppp,pppp,ppp,pp,p,cols=5)



####### DNA METHYLATION ANALYSIS ######
############################ PERFORM DNA METHYLATION ANALYSIS ############################
setwd("../../DNAmethylation")
meth <- read.csv("K786O_methylation_Only_Array_data.csv", header = T, stringsAsFactors = F)
colnames(meth) <- c("probeID", "EV_D0","EV_D5", "KO_D0","KO_D5","EV_D15","KO_D15","EV_D39","KO_D39")
anno <- read.delim("EPIC.hg38.manifest.tsv",header =T, sep = "\t",stringsAsFactors = F)
meth_anno <- merge(meth,anno,by=("probeID"))
meth_anno_CpG <- (meth_anno[,c("CpG_chrm","CpG_beg","CpG_end", "EV_D0","EV_D5", "KO_D0","KO_D5","EV_D15","KO_D15","EV_D39","KO_D39")])
meth_anno_CpG<-  meth_anno_CpG[complete.cases(meth_anno_CpG), ]
meth_anno_CpG$CpG_end <- meth_anno_CpG$CpG_end-1
library(data.table)
keys <- colnames(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end")])
X <- as.data.table(meth_anno_CpG)
meth_anno_CpG<-data.frame(X[,lapply(.SD,mean),keys])### Average CpG beta vaues

write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","EV_D0")], "EV_D0_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","EV_D5")], "EV_D5_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","EV_D15")], "EV_D15_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","EV_D39")], "EV_D39_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","KO_D0")], "KO_D0_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","KO_D5")], "KO_D5_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","KO_D15")], "KO_D15_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(meth_anno_CpG[,c("CpG_chrm","CpG_beg","CpG_end","KO_D39")], "KO_D39_CpG_Methylation_Beta_values.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#Use Deeptools to generate Fig 4A.
computeMatrix scale-regions -bs 200 -p 12 -R /genomes/hg38/gencode.v37.AutosomalTranscriptsOnly.forDeepTools.gtf -S EV_D0_CpG_Methylation_Beta_values.sorted.bw EV_D5_CpG_Methylation_Beta_values.sorted.bw EV_D15_CpG_Methylation_Beta_values.sorted.bw EV_D39_CpG_Methylation_Beta_values.sorted.bw KO_D0_CpG_Methylation_Beta_values.sorted.bw KO_D5_CpG_Methylation_Beta_values.sorted.bw KO_D15_CpG_Methylation_Beta_values.sorted.bw KO_D39_CpG_Methylation_Beta_values.sorted.bw -m 5000 -b 3000 -a 3000 -out DNAmethylation_GeneBody_O786_BetaValues.tab.gz 
plotHeatmap  --heatmapWidth 6 -m DNAmethylation_GeneBody_O786_BetaValues.tab.gz --perGroup -out DNAmethylation_GeneBody_O786_BetaValues.pdf --missingDataColor "white" -min 0 -max 1 --colorList "blue,yellow,red" --heatmapHeight 10


## Distribution of probe methylation
meth_anno_CpG <- meth_anno[meth_anno$probeType == "cg" & meth_anno$MASK_general==F,] #filter masked probes
m_meth_anno_CpG <-melt(meth_anno_CpG[,c(2:11)])
m_meth_anno_CpG$variable <- factor(m_meth_anno_CpG$variable, levels = c("EV_D0","EV_D5","EV_D15","EV_D39","KO_D0","KO_D5","KO_D15","KO_D39"))

mypalette6 <- c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797")
ppp<-ggplot(m_meth_anno_CpG, aes(x=variable, y=value)) +
  geom_violin(aes(fill=variable),adjust =1)+geom_boxplot(fill="#e7e7e7",position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA,width = 0.1) +ggtitle("DNA methylation distribution (Beta)")+ coord_cartesian(ylim=c(0, 1))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "Beta")+scale_fill_manual(values = mypalette6)+scale_color_manual(values=mypalette6)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppp #Fig S3A






