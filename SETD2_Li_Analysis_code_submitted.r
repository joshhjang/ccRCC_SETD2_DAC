#R version 4.0.2 
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


### meta data for RNA-seq samples ###
setwd("~/SETD2_project")
meta <- read.table("Cell_Sample_info_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)
meta <- meta[order(meta$Type,meta$Cell,meta$Day,meta$Rep),]
meta1 <- meta
meta1$Cell <- factor(meta1$Cell, levels =c("ACHN","O786","CAKI","A498","CAKI2","A704","P769"))
meta1$Type <- factor(meta1$Type, levels = c("WT","EV","KO","Mut"))
meta1 <- meta1[order(meta1$Cell,meta1$Type,meta1$Day,meta1$Rep),]#order cell type based on SETD2 mutation status


############# DESeq2: check differential expressed genes ###############
setwd("./Rebuttal_GeneExpression")
countData <- as.matrix(read.csv("gene_count_matrix_rebuttal.csv", row.names="gene_id"))#from Stringtie prepDE.py
colData <- read.csv("SETD2_Sample_Metadata_forDESeq2_rebuttal_final.txt", sep="\t", row.names=1)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ CellDayType)
keep <- rowSums(counts(dds)) >= 1 #filter non expressed genes 
dds <- dds[keep,]#filter lowly expressed genes
dds <- DESeq(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("Cell","Day","Rep","CellDayType")])

#Make Sample-to-sample Distance Heatmap (Fig. S2A)
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

cellcolor <- c("#2389da","#74ccf4","#5abcd8",
               "#c89f73","#f1cc8f","#fbe7a1","#d9b380",
               "#a70000","#ff5252","#ffbaba","#ff0000",
               "#317256","#49ab81","#52bf90","#398564", 
               "#ff7400","#ff9a00","#ffc100","#ff0081","#493267","#ff77bc","#e86af0","#ffcae5","#e1c4ff","#ff48a5","#9e379f","#343d46","#65737e","#a7adba","#4f5b66")

pcaData <- plotPCA(vsd, intgroup=c("CellDayType"),returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
#PCA plot for RNA-seq (Fig. S2B)
ggplot(pcaData, aes(x = PC1, y = PC2)) + 
  geom_point(size =7, aes(fill=CellDayType),shape=21) +  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))+scale_fill_manual(values = cellcolor)+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))


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

res_CAKI2WT_D0D5 <- results(dds,contrast=c("CellDayType","CAKI2D0WT","CAKI2D5WT"),lfcThreshold=1,alpha=0.01)
res_CAKI2WT_D0D15 <- results(dds,contrast=c("CellDayType","CAKI2D0WT","CAKI2D15WT"),lfcThreshold=1,alpha=0.01)
res_CAKI2WT_D0D25 <- results(dds,contrast=c("CellDayType","CAKI2D0WT","CAKI2D25WT"),lfcThreshold=1,alpha=0.01)

res_P769_D0D5 <- results(dds,contrast=c("CellDayType","P769D0WT","P769D5WT"),lfcThreshold=1,alpha=0.01)
res_P769_D0D15 <- results(dds,contrast=c("CellDayType","P769D0WT","P769D15WT"),lfcThreshold=1,alpha=0.01)
res_P769_D0D25 <- results(dds,contrast=c("CellDayType","P769D0WT","P769D25WT"),lfcThreshold=1,alpha=0.01)

res_A704_D0D5 <- results(dds,contrast=c("CellDayType","A704D0Mut","A704D5Mut"),lfcThreshold=1,alpha=0.01)
res_A704_D0D15 <- results(dds,contrast=c("CellDayType","A704D0Mut","A704D15Mut"),lfcThreshold=1,alpha=0.01)
res_A704_D0D25 <- results(dds,contrast=c("CellDayType","A704D0Mut","A704D25Mut"),lfcThreshold=1,alpha=0.01)



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

DEG_comp <- res_CAKI2WT_D0D5 
DEG_comp <- res_CAKI2WT_D0D15
DEG_comp <- res_CAKI2WT_D0D25

DEG_comp <- res_P769_D0D5
DEG_comp <- res_P769_D0D15
DEG_comp <- res_P769_D0D25

DEG_comp <- res_A704_D0D5
DEG_comp <- res_A704_D0D15
DEG_comp <- res_A704_D0D25

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

#Take top 15 positive or negative NES after filtering for hallmark pathways that are padj <0.01
dataplot <- fgseaResTidy[fgseaResTidy$padj <0.01,]


#store GSEA hallmark top 15 NES results
res_O786EV_KO_D0_HALLMARK <-data.frame(dataplot)
res_O786EV_KO_D5_HALLMARK <-data.frame(dataplot)
plotEnrichment(pathways.hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],stats=ranks) #plot for EV vs KO Day5 comparison (Fig S4C)
plotEnrichment(pathways.hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],stats=ranks) #plot for EV vs KO Day5 comparison (Fig S4C)
res_O786EV_KO_D15_HALLMARK <-data.frame(dataplot)
res_O786EV_KO_D23_HALLMARK <-data.frame(dataplot)

res_O786KO_D0D5_HALLMARK <-data.frame(dataplot)
res_O786KO_D0D15_HALLMARK <-data.frame(dataplot)
res_O786KO_D0D23_HALLMARK <-data.frame(dataplot)

res_O786EV_D0D5_HALLMARK <-data.frame(dataplot)
res_O786EV_D0D15_HALLMARK <-data.frame(dataplot)
res_O786EV_D0D23_HALLMARK <-data.frame(dataplot)

res_A498_D0D5_HALLMARK <-data.frame(dataplot)
res_A498_D0D15_HALLMARK <-data.frame(dataplot)

res_CAKI_D0D12_HALLMARK <-data.frame(dataplot)
res_CAKI_D0D25_HALLMARK <-data.frame(dataplot)

res_ACHN_D0D5_HALLMARK <-data.frame(dataplot)
res_ACHN_D0D14_HALLMARK <-data.frame(dataplot)
res_ACHN_D0D22_HALLMARK <-data.frame(dataplot)


res_CAKI2WT_D0D5_HALLMARK <-data.frame(dataplot) 
res_CAKI2WT_D0D15_HALLMARK <-data.frame(dataplot)
res_CAKI2WT_D0D25_HALLMARK <-data.frame(dataplot)

res_P769_D0D5_HALLMARK <-data.frame(dataplot)
res_P769_D0D15_HALLMARK <-data.frame(dataplot)
res_P769_D0D25_HALLMARK <-data.frame(dataplot)

res_A704_D0D5_HALLMARK <-data.frame(dataplot)
res_A704_D0D15_HALLMARK <-data.frame(dataplot)
res_A704_D0D25_HALLMARK <-data.frame(dataplot)

### #Run fGSEA for Reactome enrichment ###
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
plotEnrichment(pathways.reactome[["REACTOME_MRNA_SPLICING"]],stats=ranks) ##(Fig S7B)
res_O786EV_KO_D5_Reactome <-data.frame(dataplot_top10)
res_O786EV_KO_D15_Reactome <-data.frame(dataplot_top10)
res_O786EV_KO_D23_Reactome <-data.frame(dataplot_top10)

res_O786KO_D0D5_Reactome <-data.frame(dataplot_top10)
res_O786KO_D0D15_Reactome <-data.frame(dataplot_top10)
res_O786KO_D0D23_Reactome <-data.frame(dataplot_top10)

res_O786EV_D0D5_Reactome <-data.frame(dataplot_top10)
res_O786EV_D0D15_Reactome <-data.frame(dataplot_top10)
res_O786EV_D0D23_Reactome <-data.frame(dataplot_top10)


##### Generate combined heatmap of GSEA pathway enrichment ####
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
#Fig. 2C
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F)



#####Check if IFN-gamma, IFN-alpha, or Inflammatory responses are differentially expressed after DAC treatment across cell types
df.R<-data.frame(c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE"))
colnames(df.R) <-"pathway"
df.R<- merge(df.R,res_O786EV_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786EV_D0D23_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_O786KO_D0D23_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_A498_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_A498_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_CAKI_D0D12_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_CAKI_D0D25_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ACHN_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ACHN_D0D14_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ACHN_D0D22_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_P769_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_P769_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_P769_D0D25_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_A704_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_A704_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_A704_D0D25_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_CAKI2WT_D0D5_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_CAKI2WT_D0D15_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_CAKI2WT_D0D25_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","EV_D0_D5","EV_D0_D15","EV_D0_D23","KO_D0_D5","KO_D0_D15","KO_D0_D23","A498_D0_D5","A498_D0_D15","CAKI1_D0_D12","CAKI1_D0_D25","ACHN_D0_D5","ACHN_D0_D14","ACHN_D0_D22","P769_D0_D5","P769_D0_D15","P769_D0_D25","A704_D0_D5","A704_D0_D15","A704_D0_D25","CAKI2_D0_D5","CAKI2_D0_D15","CAKI2_D0_D25")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE"))
colnames(df.R2) <-"pathway"
df.R2<- merge(df.R2,res_O786EV_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_O786EV_D0D15_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_O786EV_D0D23_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_O786KO_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_O786KO_D0D15_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_O786KO_D0D23_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_A498_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_A498_D0D15_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_CAKI_D0D12_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_CAKI_D0D25_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ACHN_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ACHN_D0D14_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ACHN_D0D22_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_P769_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_P769_D0D15_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_P769_D0D25_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_A704_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_A704_D0D15_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_A704_D0D25_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_CAKI2WT_D0D5_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_CAKI2WT_D0D15_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_CAKI2WT_D0D25_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
colnames(df.R2)<-c("pathway","EV_D0_D5","EV_D0_D15","EV_D0_D23","KO_D0_D5","KO_D0_D15","KO_D0_D23","A498_D0_D5","A498_D0_D15","CAKI1_D0_D12","CAKI1_D0_D25","ACHN_D0_D5","ACHN_D0_D14","ACHN_D0_D22","P769_D0_D5","P769_D0_D15","P769_D0_D25","A704_D0_D5","A704_D0_D15","A704_D0_D25","CAKI2_D0_D5","CAKI2_D0_D15","CAKI2_D0_D25")
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:23)]<- -log(df.R2[,c(2:23)])

df.R_M <- melt(df.R[,c(24,2:23)])
df.R2_M <- melt(df.R2[,c(24,2:23)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("EV_D0_D5","EV_D0_D15","EV_D0_D23","KO_D0_D5","KO_D0_D15","KO_D0_D23","A498_D0_D5","A498_D0_D15","CAKI1_D0_D12","CAKI1_D0_D25","ACHN_D0_D5","ACHN_D0_D14","ACHN_D0_D22","P769_D0_D5","P769_D0_D15","P769_D0_D25","A704_D0_D5","A704_D0_D15","A704_D0_D25","CAKI2_D0_D5","CAKI2_D0_D15","CAKI2_D0_D25"))  
#Fig. S3 
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =c(brewer.pal(n = 5, name = "BrBG")),limits=c(-2, 2)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




#####786-O EV vs KO GSEA Reactome Enrichment (Figure S7A)
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

#786-O EV across timepoint comparison & 786-O KO across timepoint, GSEA Reactome Enrichment present in >=2 samples (Figure S7C)
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

name <- df.R[rowSums(df.R[,2:7] != 0) >= 2,]$pathway_simple
df.OG2 <- data.matrix(df.R[rowSums(df.R[,2:7] != 0) >= 2, 2:7])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(df.R[,2:7])
pheatmap(df.OG2, color = colorRampPalette((brewer.pal(n = 5, name = "BrBG")))(100),breaks = c(seq(-2,2,by=0.04)),cluster_rows =T, cluster_cols=F,clustering_method ="ward.D2")




###########Save normalized count values for Zscore heatmaps ###########
cdds<-data.frame(counts(dds, normalized = TRUE))
cdds$GeneID <- rownames(cdds)
cdds$Gene.name <- apply(cdds, 1, function(x) {unlist(strsplit(x[61],"\\|"))[2]})
colnames(cdds)[61:62] <- c("Gene","GeneID")
cdds<-cdds[,c("Gene","GeneID",meta1$Sample)]

####SETD2 Gene expression (Fig. S2C)####
SETD2_gene <-cdds[which(cdds$GeneID %in% c("SETD2")),]
SETD2_gene <-SETD2_gene[order(SETD2_gene$GeneID),]

SETD2_gene1<-log2(SETD2_gene[,c(3:62)]+0.01)

SETD2info <- data.frame(c(colnames(SETD2_gene1)),c(t(SETD2_gene[,c(3:62)])))
colnames(SETD2info) <- c("Sample","GEX")
SETD2info <- merge(SETD2info,meta1[,c(1:6)],by = "Sample")

SETD2info$Day <-factor(SETD2info$Day, levels = c("0","5","12","14","15","22","23","25"))
SETD2info$Cell <-factor(SETD2info$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A498","A704"))
SETD2info$CellDayType <- paste(SETD2info$Cell,SETD2info$Day,SETD2info$Type,sep="_")
SETD2info$CellType <- paste(SETD2info$Cell,SETD2info$Type,sep="_")
SETD2info$CellType <-factor(SETD2info$CellType, levels = c("A498_Mut","A704_Mut","ACHN_WT","CAKI_Mut","CAKI2_WT","O786_EV","O786_KO","P769_WT"))

cellcolor <- c("#2389da","#74ccf4","#5abcd8","#2389da",
               "#c89f73","#f1cc8f","#fbe7a1","#d9b380", "#c89f73",
               "#a70000","#ff5252","#ffbaba","#ff0000", "#ff5252","#ff7400","#ff9a00","#ffc100","#ff7400",
               "#317256","#49ab81","#52bf90","#398564", "#5cb85c",
               "#ff0081","#493267","#ff77bc","#e86af0","#ffcae5","#e1c4ff","#ff48a5","#9e379f","#ff0081","#493267","#343d46","#65737e","#a7adba","#4f5b66","#65737e")


ggplot(SETD2info)+geom_boxplot(aes(x=Cell, y=log2(GEX), fill= CellType),width=0.6,outlier.shape = NA)+geom_point(aes(x=Cell, y=log2(GEX), fill=CellDayType),color = "black",position = position_jitter(0.2),pch = 21,size=5) +ggtitle("SETD2 Expression (Normalized Counts)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Cell", y = "log2(Normalized Counts)")+scale_fill_manual(values = cellcolor)+scale_color_manual(values= cellcolor)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))

#################### Gene pathway Z-score heatmaps ##########################
##IFNy Gene expression (Fig 3A, S4A)
IFNy <- read.table("../fGSEA/GSEA_Hallmark_interferon_gamma_response.txt", sep = "\t", header = F, stringsAsFactors = F)
colnames(IFNy) <-"GeneID"

SETD2_IFNy <-cdds[which(cdds$GeneID %in% IFNy$GeneID==T),]
SETD2_IFNy <-SETD2_IFNy[order(SETD2_IFNy$GeneID),]

SETD2_IFNy1<-log2(SETD2_IFNy[,c(3:62)]+0.01)#Zscore

SETD2_IFNy1_Zscore<- (SETD2_IFNy1-rowMeans(SETD2_IFNy1))/(rowSds(as.matrix(SETD2_IFNy1)))[row(SETD2_IFNy1)]
name <- SETD2_IFNy[,c("GeneID")]
df.OG2 <- data.matrix(SETD2_IFNy1_Zscore)
row.names(df.OG2) <- name
d<-pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =T, cluster_cols=F)

name <- SETD2_IFNy[row_order(d),c("GeneID")]
df.OG2 <- data.matrix(SETD2_IFNy1[row_order(d),]) #normCounts heatmap, fix roworder based on clustering of 
pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),breaks = c(seq(0,16,by=0.16)),show_rownames = T,cluster_rows =F, cluster_cols=F)


#inflammatory response (Fig 3B, S4B)
inflam <- read.table("../fGSEA/GSEA_Hallmark_inflammatory_response.txt", sep = "\t", header = F, stringsAsFactors = F)

colnames(inflam) <-"GeneID"

SETD2_inflam <-cdds[which(cdds$GeneID %in% inflam$GeneID==T),]
SETD2_inflam <-SETD2_inflam[order(SETD2_inflam$GeneID),]

SETD2_inflam1<-log2(SETD2_inflam[,c(3:62)]+0.01)#Zscore

SETD2_inflam1_Zscore<- (SETD2_inflam1-rowMeans(SETD2_inflam1))/(rowSds(as.matrix(SETD2_inflam1)))[row(SETD2_inflam1)]
name <- SETD2_inflam[,c("GeneID")]
df.OG2 <- data.matrix(SETD2_inflam1_Zscore)
row.names(df.OG2) <- name
d<-pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =T, cluster_cols=F)

name <- SETD2_inflam[row_order(d),c("GeneID")]
df.OG2 <- data.matrix(SETD2_inflam1[row_order(d),]) #normCounts heatmap, fix roworder based on clustering of 
pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),breaks = c(seq(0,16,by=0.16)),show_rownames = T,cluster_rows =F, cluster_cols=F)



####partition into High, medium, low genes in 786-O (for Fig. S6A)
cdds_786O <- cdds[,c("OEVD0","O786WD0","OEVD5","O786WD5","OEVD15","O786WD15","OEVD23","O786WD23","SK0D0","O786MD0","SK0D5","O786MD5","SK0D15","O786MD15","SK0D23","O786MD23","GeneID")]
cdds_786O$OEVD0_ave <- rowMeans(cdds_786O[,c(1:2)])
cdds_786O$OEVD5_ave <- rowMeans(cdds_786O[,c(3:4)])
cdds_786O$OEVD15_ave <- rowMeans(cdds_786O[,c(5:6)])
cdds_786O$OEVD23_ave <- rowMeans(cdds_786O[,c(7:8)])
cdds_786O$SK0D0_ave <- rowMeans(cdds_786O[,c(9:10)])
cdds_786O$SK0D5_ave <- rowMeans(cdds_786O[,c(11:12)])
cdds_786O$SK0D15_ave <- rowMeans(cdds_786O[,c(13:14)])
cdds_786O$SK0D23_ave <- rowMeans(cdds_786O[,c(15:16)])
cdds_786O_ave <- cdds_786O[,c(17:25)]
cdds_786O_ave$Allsample_ave <- rowMeans(cdds_786O_ave[,2:8])
m_cdds_786O_ave <-reshape2::melt(cdds_786O_ave[,2:10])

g<-ggplot(m_cdds_786O_ave, aes(x=log(value+0.001), color = variable)) +
  geom_density(lwd=2,alpha = 0.5)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "log(normalized Counts)", y = "Density")+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
#Fig. S6A
g+geom_vline(xintercept = -5,linetype="dotted")+geom_vline(xintercept=1,linetype="dotted")+geom_vline(xintercept=5,linetype="dotted")


cdds_786O_ave_NoExp <- cdds_786O_ave[log(cdds_786O_ave$Allsample_ave) < -5,] #6691
cdds_786O_ave_low <- cdds_786O_ave[log(cdds_786O_ave$Allsample_ave) > -5 & log(cdds_786O_ave$Allsample_ave) < 1,] #15944
cdds_786O_ave_middle <- cdds_786O_ave[log(cdds_786O_ave$Allsample_ave) > 1 & log(cdds_786O_ave$Allsample_ave) < 5,] #13607
cdds_786O_ave_high <- cdds_786O_ave[log(cdds_786O_ave$Allsample_ave) >=5,]#13215

write.table(cdds_786O_ave_NoExp[,c("GeneID","Allsample_ave")],"cdds_786O_ave_NoExp.txt",col.names = T, row.names = F, sep ="\t", quote = F)
write.table(cdds_786O_ave_low[,c("GeneID","Allsample_ave")],"cdds_786O_ave_low.txt",col.names = T, row.names = F, sep ="\t", quote = F)
write.table(cdds_786O_ave_middle[,c("GeneID","Allsample_ave")],"cdds_786O_ave_middle.txt",col.names = T, row.names = F, sep ="\t", quote = F)
write.table(cdds_786O_ave_high[,c("GeneID","Allsample_ave")],"cdds_786O_ave_high.txt",col.names = T, row.names = F, sep ="\t", quote = F)
#Use python to parse out gtf based on gene expression and use Deeptools with histone data to generate Fig. S6B.




##mRNA splicing gene Zscore heatmap, only 786-O (Fig. S7D)
splicing <- read.table("../fGSEA/REACTOME_mRNA_Splicing.txt", sep = "\t")
SETD2_SPL <-cdds[which(cdds$GeneID %in% splicing$V1==T),]
SETD2_SPL <-SETD2_SPL[order(SETD2_SPL$GeneID),]

SETD2_SPL1<-log2(SETD2_SPL[,c(11:26)]+0.01)#Zscore for 786-O
df.OG2 <- data.matrix(SETD2_SPL1) #normCounts heatmap, fix roworder based on clustering of 
d<-pheatmap(df.OG2,clustering_method = "ward.D2",color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),breaks = c(seq(5,16,by=0.11)),show_rownames = T,cluster_rows =T, cluster_cols=F)


#fold-change 786-O
SETD2_SPL2<-SETD2_SPL[row_order(d),c(11:26)]
SETD2_SPL2$OEVD0_ave <- rowMeans(SETD2_SPL2[,c(1:2)])
SETD2_SPL2$OEVD5_ave <- rowMeans(SETD2_SPL2[,c(3:4)])
SETD2_SPL2$OEVD15_ave <- rowMeans(SETD2_SPL2[,c(5:6)])
SETD2_SPL2$OEVD23_ave <- rowMeans(SETD2_SPL2[,c(7:8)])
SETD2_SPL2$SK0D0_ave <- rowMeans(SETD2_SPL2[,c(9:10)])
SETD2_SPL2$SK0D5_ave <- rowMeans(SETD2_SPL2[,c(11:12)])
SETD2_SPL2$SK0D15_ave <- rowMeans(SETD2_SPL2[,c(13:14)])
SETD2_SPL2$SK0D23_ave <- rowMeans(SETD2_SPL2[,c(15:16)])
write.table(SETD2_SPL2[,17:24]+0.01, "SETD2_splicing_normCountsAve_V37anno_Days.txt", sep = "\t", col.names = F, row.names = F, quote = F)
#make fold change: awk -v OFS="\t" '{print $5/$1,$1/$1, $2/$1, $3/$1, $4/$1, $5/$5, $6/$5, $7/$5, $8/$5}' SETD2_splicing_normCountsAve_V37anno_Days.txt > SETD2_splicing_normCountsAve_V37anno_FoldChange.txt
SETD2_SPL_FC <- read.table("SETD2_splicing_normCountsAve_V37anno_FoldChange.txt",sep ="\t", header = F, stringsAsFactors = F)
SETD2_SPL_FC <- log2(SETD2_SPL_FC)
pheatmap(data.matrix(SETD2_SPL_FC),color = colorRampPalette(rev(brewer.pal(n = 9, name = "PiYG")))(100),breaks = c(seq(-1,1,by=0.02)),show_rownames = F,cluster_rows =F, cluster_cols=F)




###################### Transposable Element Analysis ################################
setwd("../TEexpression")
TE_dic <- read.table("GRCh38_GENCODE_rmsk_TE_dictionary.txt",sep ="\t", header = F, stringsAsFactors = F)
colnames(TE_dic) <- c("TE_Subfamily","Geneid","TE_Family","TE_Class")
TE_dic <- TE_dic[!duplicated(TE_dic),]

TE_annot <- read.delim("GRCh38_GENCODE_rmsk_TE_GTF_parsed.V37annotate.bed", sep ="\t", header =T, stringsAsFactors = F)
TE_annot <- TE_annot[,c(2,3,4,8)]
TE_annot$Anno2 <- apply(TE_annot, 1, function(x) {unlist(strsplit(x[4], " \\("))[1]})
TE_annot$Start <- TE_annot$Start-1 #make into 0-based


setwd("../Rebuttal_GeneExpression")
TE_cpm <- read.table("SETD2_TE_expression_all_sample_CPM_rebuttal.bed", sep ="\t", header = T, stringsAsFactors = F)
TE_cpm <- TE_cpm[,c(colnames(TE_cpm[1:6]),meta1$Sample)] 

TE_cpm <- merge(TE_cpm, TE_dic, by = "Geneid", x.all =T)
TE_cpm <- TE_cpm[TE_cpm$TE_Class == "LTR" | TE_cpm$TE_Class == "LINE" | TE_cpm$TE_Class == "SINE" | TE_cpm$TE_Class == "DNA" | TE_cpm$TE_Class == "SINE" ,] #only LINE, SINE, LTR, DNA

TE_cpm_anno_all<- merge(TE_cpm,TE_annot, by=c("Chr","Start","End")) 
TE_cpm_anno_all <-TE_cpm_anno_all[,c(1:6,67,68,69,71,7:66)]


filter_1cpm <- which(apply(TE_cpm_anno_all[,11:70], 1, function (x) (max(x)>=1))) ## 1CPM filter
TE_cpm_anno <- TE_cpm_anno_all[filter_1cpm,] 
TE_cpm_anno <-TE_cpm_anno[TE_cpm_anno$Anno2 == "intron" | TE_cpm_anno$Anno2 == "Intergenic",] #only look at intergenic or intronic TEs
#Visualize with heatmaps to verify replicability of biological replicates before averaging biological replicates

TE_cpm_anno_intron <-TE_cpm_anno[TE_cpm_anno$Anno2 == "intron",]
TE_cpm_anno_intron<-TE_cpm_anno_intron[order(TE_cpm_anno_intron$TE_Class,TE_cpm_anno_intron$TE_Family,TE_cpm_anno_intron$TE_Subfamily),] #sort by TE name

TE_cpm_anno_intron_LTR <- TE_cpm_anno_intron[TE_cpm_anno_intron$TE_Class == "LTR",] 
TE_cpm_anno_intron_SINE<- TE_cpm_anno_intron[TE_cpm_anno_intron$TE_Class == "SINE",] 
TE_cpm_anno_intron_LINE<- TE_cpm_anno_intron[TE_cpm_anno_intron$TE_Class == "LINE",] 

TE_cpm_anno_intergenic <-TE_cpm_anno[TE_cpm_anno$Anno2 == "Intergenic",]
TE_cpm_anno_intergenic<-TE_cpm_anno_intergenic[order(TE_cpm_anno_intergenic$TE_Class,TE_cpm_anno_intergenic$TE_Family,TE_cpm_anno_intergenic$TE_Subfamily),] #sort by TE name

TE_cpm_anno_intergenic_LTR <- TE_cpm_anno_intergenic[TE_cpm_anno_intergenic$TE_Class == "LTR",] 
TE_cpm_anno_intergenic_SINE<- TE_cpm_anno_intergenic[TE_cpm_anno_intergenic$TE_Class == "SINE",] 
TE_cpm_anno_intergenic_LINE<- TE_cpm_anno_intergenic[TE_cpm_anno_intergenic$TE_Class == "LINE",] 


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
TE_cpm_anno$CaSCRD0_ave <- rowMeans(TE_cpm_anno[,c(47:48)])
TE_cpm_anno$CaSCRD5_ave <- rowMeans(TE_cpm_anno[,c(49:50)])
TE_cpm_anno$CaSCRD15_ave <- rowMeans(TE_cpm_anno[,c(51:52)])
TE_cpm_anno$CaSCRD25_ave <- rowMeans(TE_cpm_anno[,c(53:54)])
TE_cpm_anno$A704D0_ave <- rowMeans(TE_cpm_anno[,c(55:56)])
TE_cpm_anno$A704D5_ave <- rowMeans(TE_cpm_anno[,c(57:58)])
TE_cpm_anno$A704D15_ave <- rowMeans(TE_cpm_anno[,c(59:60)])
TE_cpm_anno$A704D25_ave <- rowMeans(TE_cpm_anno[,c(61:62)])
TE_cpm_anno$P769D0_ave <- rowMeans(TE_cpm_anno[,c(63:64)])
TE_cpm_anno$P769D5_ave <- rowMeans(TE_cpm_anno[,c(65:66)])
TE_cpm_anno$P769D15_ave <- rowMeans(TE_cpm_anno[,c(67:68)])
TE_cpm_anno$P769D25_ave <- rowMeans(TE_cpm_anno[,c(69:70)])

TE_cpm_anno_ave_ACHN <- TE_cpm_anno[,c(1:10,71:74)]
TE_cpm_anno_ave_O786 <- TE_cpm_anno[,c(1:10,75:82)]
TE_cpm_anno_ave_CAKI <- TE_cpm_anno[,c(1:10,83:85)]
TE_cpm_anno_ave_A498 <- TE_cpm_anno[,c(1:10,86:88)]
TE_cpm_anno_ave_CAKI2 <- TE_cpm_anno[,c(1:10,89:92)]
TE_cpm_anno_ave_A704 <- TE_cpm_anno[,c(1:10,93:96)]
TE_cpm_anno_ave_P769 <- TE_cpm_anno[,c(1:10,97:100)]

filter_1cpm <- which(apply(TE_cpm_anno_ave_ACHN[,11:14], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_ACHN <- TE_cpm_anno_ave_ACHN[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_O786[,11:18], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_O786<- TE_cpm_anno_ave_O786[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_CAKI[,11:13], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_CAKI <- TE_cpm_anno_ave_CAKI[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_A498[,11:13], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_A498 <- TE_cpm_anno_ave_A498[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_CAKI2[,11:14], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_CAKI2<- TE_cpm_anno_ave_CAKI2[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_A704[,11:14], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_A704 <- TE_cpm_anno_ave_A704[filter_1cpm,]
filter_1cpm <- which(apply(TE_cpm_anno_ave_P769[,11:14], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave_P769 <- TE_cpm_anno_ave_P769[filter_1cpm,]

TE_cpm_anno_ave <- TE_cpm_anno[,c(1:10,63:88)]
filter_1cpm <- which(apply(TE_cpm_anno_ave[,11:36], 1, function (x) (max(x)>=1)))
TE_cpm_anno_ave <- TE_cpm_anno_ave[filter_1cpm,] 


##Calculate TE expression fold change
write.table(TE_cpm_anno_ave_ACHN[,11:14]+0.01, "SETD2_TE_averageCPM_V37anno_Days_ACHN.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_O786[,11:18]+0.01, "SETD2_TE_averageCPM_V37anno_Days_O786.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_CAKI[,11:13]+0.01, "SETD2_TE_averageCPM_V37anno_Days_CAKI.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_A498[,11:13]+0.01, "SETD2_TE_averageCPM_V37anno_Days_A498.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_CAKI2[,11:14]+0.01, "SETD2_TE_averageCPM_V37anno_Days_CAKI2.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_A704[,11:14]+0.01, "SETD2_TE_averageCPM_V37anno_Days_A704.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(TE_cpm_anno_ave_P769[,11:14]+0.01, "SETD2_TE_averageCPM_V37anno_Days_P769.txt", sep = "\t", col.names = F, row.names = F, quote = F)

#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1}' SETD2_TE_averageCPM_V37anno_Days_ACHN.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_ACHN.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1,$5/$5,$6/$5,$7/$5,$8/$5}' SETD2_TE_averageCPM_V37anno_Days_O786.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_O786.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1}' SETD2_TE_averageCPM_V37anno_Days_CAKI.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_CAKI.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1}' SETD2_TE_averageCPM_V37anno_Days_A498.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_A498.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1}' SETD2_TE_averageCPM_V37anno_Days_CAKI2.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_CAKI2.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1}' SETD2_TE_averageCPM_V37anno_Days_A704.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_A704.txt
#make fold change: awk -v OFS="\t" '{print $1/$1, $2/$1, $3/$1,$4/$1}' SETD2_TE_averageCPM_V37anno_Days_P769.txt > SETD2_TE_averageCPM_V37anno_Days_FoldChange_P769.txt

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


TE_cpm_anno_ave_FC_CAKI2 <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_CAKI2.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_CAKI2 <- cbind(TE_cpm_anno_ave_CAKI2[,c(1:10)],TE_cpm_anno_ave_FC_CAKI2)
colnames(TE_cpm_anno_ave_FC_CAKI2) <- colnames(TE_cpm_anno_ave_CAKI2)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_CAKI2[,11:14], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_CAKI2 <- TE_cpm_anno_ave_FC_CAKI2[filter_2FC,]

TE_cpm_anno_ave_FC2_both_CAKI2<-TE_cpm_anno_ave_FC2_CAKI2 #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_CAKI2<-TE_cpm_anno_ave_FC2_both_CAKI2[order(TE_cpm_anno_ave_FC2_both_CAKI2$TE_Class,TE_cpm_anno_ave_FC2_both_CAKI2$TE_Family,TE_cpm_anno_ave_FC2_both_CAKI2$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_CAKI2 <-TE_cpm_anno_ave_FC2_CAKI2[TE_cpm_anno_ave_FC2_CAKI2$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_CAKI2<-TE_cpm_anno_ave_FC2_intron_CAKI2[order(TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class,TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Family,TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_CAKI2 <-TE_cpm_anno_ave_FC2_CAKI2[TE_cpm_anno_ave_FC2_CAKI2$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_CAKI2<-TE_cpm_anno_ave_FC2_intergenic_CAKI2[order(TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class,TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Family,TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC_A704 <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_A704.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_A704 <- cbind(TE_cpm_anno_ave_A704[,c(1:10)],TE_cpm_anno_ave_FC_A704)
colnames(TE_cpm_anno_ave_FC_A704) <- colnames(TE_cpm_anno_ave_A704)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_A704[,11:14], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_A704 <- TE_cpm_anno_ave_FC_A704[filter_2FC,]

TE_cpm_anno_ave_FC2_both_A704<-TE_cpm_anno_ave_FC2_A704 #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_A704<-TE_cpm_anno_ave_FC2_both_A704[order(TE_cpm_anno_ave_FC2_both_A704$TE_Class,TE_cpm_anno_ave_FC2_both_A704$TE_Family,TE_cpm_anno_ave_FC2_both_A704$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_A704 <-TE_cpm_anno_ave_FC2_A704[TE_cpm_anno_ave_FC2_A704$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_A704<-TE_cpm_anno_ave_FC2_intron_A704[order(TE_cpm_anno_ave_FC2_intron_A704$TE_Class,TE_cpm_anno_ave_FC2_intron_A704$TE_Family,TE_cpm_anno_ave_FC2_intron_A704$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_A704 <-TE_cpm_anno_ave_FC2_A704[TE_cpm_anno_ave_FC2_A704$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_A704<-TE_cpm_anno_ave_FC2_intergenic_A704[order(TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class,TE_cpm_anno_ave_FC2_intergenic_A704$TE_Family,TE_cpm_anno_ave_FC2_intergenic_A704$TE_Subfamily),] #sort by TE name


TE_cpm_anno_ave_FC_P769 <- read.table("SETD2_TE_averageCPM_V37anno_Days_FoldChange_P769.txt",sep ="\t", header = F, stringsAsFactors = F)
TE_cpm_anno_ave_FC_P769 <- cbind(TE_cpm_anno_ave_P769[,c(1:10)],TE_cpm_anno_ave_FC_P769)
colnames(TE_cpm_anno_ave_FC_P769) <- colnames(TE_cpm_anno_ave_P769)

filter_2FC <- which(apply(TE_cpm_anno_ave_FC_P769[,11:14], 1, function (x) (max(x)>=2))) 
TE_cpm_anno_ave_FC2_P769 <- TE_cpm_anno_ave_FC_P769[filter_2FC,]

TE_cpm_anno_ave_FC2_both_P769<-TE_cpm_anno_ave_FC2_P769 #both Integenic and intronic TEs
TE_cpm_anno_ave_FC2_both_P769<-TE_cpm_anno_ave_FC2_both_P769[order(TE_cpm_anno_ave_FC2_both_P769$TE_Class,TE_cpm_anno_ave_FC2_both_P769$TE_Family,TE_cpm_anno_ave_FC2_both_P769$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intron_P769 <-TE_cpm_anno_ave_FC2_P769[TE_cpm_anno_ave_FC2_P769$Anno2 == "intron",]
TE_cpm_anno_ave_FC2_intron_P769<-TE_cpm_anno_ave_FC2_intron_P769[order(TE_cpm_anno_ave_FC2_intron_P769$TE_Class,TE_cpm_anno_ave_FC2_intron_P769$TE_Family,TE_cpm_anno_ave_FC2_intron_P769$TE_Subfamily),] #sort by TE name

TE_cpm_anno_ave_FC2_intergenic_P769 <-TE_cpm_anno_ave_FC2_P769[TE_cpm_anno_ave_FC2_P769$Anno2 == "Intergenic",]
TE_cpm_anno_ave_FC2_intergenic_P769<-TE_cpm_anno_ave_FC2_intergenic_P769[order(TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class,TE_cpm_anno_ave_FC2_intergenic_P769$TE_Family,TE_cpm_anno_ave_FC2_intergenic_P769$TE_Subfamily),] #sort by TE name



#####Bar plot of number of TEs that are >2FC (Fig. 3D)#######
bothTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a", "P769D5a","P769D15a","P769D25a"),
                       c(nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_O786[TE_cpm_anno_ave_FC2_both_O786$SK0D23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A498[TE_cpm_anno_ave_FC2_both_A498$A498MA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A498[TE_cpm_anno_ave_FC2_both_A498$A498MA15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI[TE_cpm_anno_ave_FC2_both_CAKI$CAKIWA12_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI[TE_cpm_anno_ave_FC2_both_CAKI$CAKIWA25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_ACHN[TE_cpm_anno_ave_FC2_both_ACHN$ACHNWA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_ACHN[TE_cpm_anno_ave_FC2_both_ACHN$ACHNWA14_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_ACHN[TE_cpm_anno_ave_FC2_both_ACHN$ACHNWA22_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI2[TE_cpm_anno_ave_FC2_both_CAKI2$CaSCRD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI2[TE_cpm_anno_ave_FC2_both_CAKI2$CaSCRD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_CAKI2[TE_cpm_anno_ave_FC2_both_CAKI2$CaSCRD25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A704[TE_cpm_anno_ave_FC2_both_A704$A704D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A704[TE_cpm_anno_ave_FC2_both_A704$A704D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_A704[TE_cpm_anno_ave_FC2_both_A704$A704D25_ave >=2,]), nrow(TE_cpm_anno_ave_FC2_both_P769[TE_cpm_anno_ave_FC2_both_P769$P769D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_P769[TE_cpm_anno_ave_FC2_both_P769$P769D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_both_P769[TE_cpm_anno_ave_FC2_both_P769$P769D25_ave >=2,])))
colnames(bothTE_FC2) <-c("Sample","Count")
bothTE_FC2<-merge(meta1[,c(1:6)],bothTE_FC2, by="Sample")
bothTE_FC2$Day <-factor(bothTE_FC2$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2$Cell <-factor(bothTE_FC2$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A704","A498"))
bothTE_FC2$CellType <- paste(bothTE_FC2$Cell,bothTE_FC2$Type,sep="_")
bothTE_FC2$CellType <-factor(bothTE_FC2$CellType, levels = c("ACHN_WT","O786_EV","O786_KO","CAKI_Mut","A498_Mut","CAKI2_WT","A704_Mut","P769_WT"))
bothTE_FC2<-bothTE_FC2[order(bothTE_FC2$CellType),]

p <- ggplot(bothTE_FC2, aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic+intron TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+ scale_y_continuous(limits = c(0, max(bothTE_FC2$Count)))+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(400,1700)) #save 5 x 14

####Make Barplot separating intron and intergenic TE >2FC
bothTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a", "P769D5a","P769D15a","P769D25a","OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a","P769D5a","P769D15a","P769D25a"),
                       c(rep("Intergenic",22),rep("Intronic",22)),
                       c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA15_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA12_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA25_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA14_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA22_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D25_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA15_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA12_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA25_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA14_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA22_ave >=2,]),
                         nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D25_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D25_ave >=2,])))
colnames(bothTE_FC2) <-c("Sample","anno","Count")
bothTE_FC2<-merge(meta1[,c(1:6)],bothTE_FC2, by="Sample")
bothTE_FC2$Day <-factor(bothTE_FC2$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2$Cell <-factor(bothTE_FC2$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A704","A498"))
bothTE_FC2$CellType <- paste(bothTE_FC2$Cell,bothTE_FC2$Type,sep="_")
bothTE_FC2$CellType <-factor(bothTE_FC2$CellType, levels = c("ACHN_WT","O786_EV","O786_KO","CAKI_Mut","A498_Mut","CAKI2_WT","A704_Mut","P769_WT"))
bothTE_FC2<-bothTE_FC2[order(bothTE_FC2$CellType),]


#Fig. S8A
p <- ggplot(bothTE_FC2[bothTE_FC2$anno == "Intergenic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("IntergenicTE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p <- p+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,1500))
p1 <- ggplot(bothTE_FC2[bothTE_FC2$anno == "Intronic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intron TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p1 <- p1+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,1500))
multiplot(p1,p,cols=1)

##########    Bar graph for TE class specific occurrence (Fig. S8B)   ################
bothTE_FC2_LTR<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a", "P769D5a","P769D15a","P769D25a","OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a","P769D5a","P769D15a","P769D25a"),
                           c(rep("Intergenic",22),rep("Intronic",22)),
                           c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "LTR",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "LTR",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "LTR",]) ))
colnames(bothTE_FC2_LTR) <-c("Sample","anno","Count")
bothTE_FC2_LTR<-merge(meta1[,c(1:6)],bothTE_FC2_LTR, by="Sample")
bothTE_FC2_LTR$Day <-factor(bothTE_FC2_LTR$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2_LTR$Cell <-factor(bothTE_FC2_LTR$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A704","A498"))
bothTE_FC2_LTR$CellType <- paste(bothTE_FC2_LTR$Cell,bothTE_FC2_LTR$Type,sep="_")
bothTE_FC2_LTR$CellType <-factor(bothTE_FC2_LTR$CellType, levels = c("ACHN_WT","O786_EV","O786_KO","CAKI_Mut","A498_Mut","CAKI2_WT","A704_Mut","P769_WT"))
bothTE_FC2_LTR<-bothTE_FC2_LTR[order(bothTE_FC2_LTR$CellType),]

p <- ggplot(bothTE_FC2_LTR[bothTE_FC2_LTR$anno == "Intergenic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic LTR >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p <-p+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")
p1 <- ggplot(bothTE_FC2_LTR[bothTE_FC2_LTR$anno == "Intronic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intron LTR >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p1 <-p1+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")

bothTE_FC2_LINE<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a", "P769D5a","P769D15a","P769D25a","OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a","P769D5a","P769D15a","P769D25a"),
                            c(rep("Intergenic",22),rep("Intronic",22)),
                            c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "LINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "LINE",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "LINE",]) ))
colnames(bothTE_FC2_LINE) <-c("Sample","anno","Count")
bothTE_FC2_LINE<-merge(meta1[,c(1:6)],bothTE_FC2_LINE, by="Sample")
bothTE_FC2_LINE$Day <-factor(bothTE_FC2_LINE$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2_LINE$Cell <-factor(bothTE_FC2_LINE$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A704","A498"))
bothTE_FC2_LINE$CellType <- paste(bothTE_FC2_LINE$Cell,bothTE_FC2_LINE$Type,sep="_")
bothTE_FC2_LINE$CellType <-factor(bothTE_FC2_LINE$CellType, levels = c("ACHN_WT","O786_EV","O786_KO","CAKI_Mut","A498_Mut","CAKI2_WT","A704_Mut","P769_WT"))
bothTE_FC2_LINE<-bothTE_FC2_LINE[order(bothTE_FC2_LINE$CellType),]

pa <- ggplot(bothTE_FC2_LINE[bothTE_FC2_LINE$anno == "Intergenic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic LINE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pa <-pa+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")
pa1 <- ggplot(bothTE_FC2_LINE[bothTE_FC2_LINE$anno == "Intronic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intron LINE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pa1 <-pa1+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")



bothTE_FC2_SINE<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a", "P769D5a","P769D15a","P769D25a","OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a","P769D5a","P769D15a","P769D25a"),
                            c(rep("Intergenic",22),rep("Intronic",22)),
                            c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "SINE",]),
                              nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "SINE",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "SINE",]) ))
colnames(bothTE_FC2_SINE) <-c("Sample","anno","Count")
bothTE_FC2_SINE<-merge(meta1[,c(1:6)],bothTE_FC2_SINE, by="Sample")
bothTE_FC2_SINE$Day <-factor(bothTE_FC2_SINE$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2_SINE$Cell <-factor(bothTE_FC2_SINE$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A704","A498"))
bothTE_FC2_SINE$CellType <- paste(bothTE_FC2_SINE$Cell,bothTE_FC2_SINE$Type,sep="_")
bothTE_FC2_SINE$CellType <-factor(bothTE_FC2_SINE$CellType, levels = c("ACHN_WT","O786_EV","O786_KO","CAKI_Mut","A498_Mut","CAKI2_WT","A704_Mut","P769_WT"))
bothTE_FC2_SINE<-bothTE_FC2_SINE[order(bothTE_FC2_SINE$CellType),]

pb <- ggplot(bothTE_FC2_SINE[bothTE_FC2_SINE$anno == "Intergenic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic SINE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pb <-pb+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")
pb1 <- ggplot(bothTE_FC2_SINE[bothTE_FC2_SINE$anno == "Intronic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intron SINE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pb1 <-pb1+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")



bothTE_FC2_DNA<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a", "P769D5a","P769D15a","P769D25a","OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23","A498MA5","A498MA15","CAKIWA12","CAKIWA25","ACHNWA5","ACHNWA14","ACHNWA22","CaSCRD5a","CaSCR15a","CaSCR25a","A704D5a","A704D15a","A704D25a","P769D5a","P769D15a","P769D25a"),
                           c(rep("Intergenic",22),rep("Intronic",22)),
                           c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_O786$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A498[TE_cpm_anno_ave_FC2_intergenic_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A498$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI[TE_cpm_anno_ave_FC2_intergenic_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_ACHN[TE_cpm_anno_ave_FC2_intergenic_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_ACHN$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_CAKI2[TE_cpm_anno_ave_FC2_intergenic_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_CAKI2$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_A704[TE_cpm_anno_ave_FC2_intergenic_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_A704$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intergenic_P769[TE_cpm_anno_ave_FC2_intergenic_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intergenic_P769$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2 & TE_cpm_anno_ave_FC2_intron_O786$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_A498[TE_cpm_anno_ave_FC2_intron_A498$A498MA15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A498$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA12_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI[TE_cpm_anno_ave_FC2_intron_CAKI$CAKIWA25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA5_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA14_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_ACHN[TE_cpm_anno_ave_FC2_intron_ACHN$ACHNWA22_ave >=2 & TE_cpm_anno_ave_FC2_intron_ACHN$TE_Class == "DNA",]),
                             nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD5_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD15_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_CAKI2[TE_cpm_anno_ave_FC2_intron_CAKI2$CaSCRD25_ave >=2 & TE_cpm_anno_ave_FC2_intron_CAKI2$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_A704[TE_cpm_anno_ave_FC2_intron_A704$A704D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_A704$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D5_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D15_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "DNA",]),nrow(TE_cpm_anno_ave_FC2_intron_P769[TE_cpm_anno_ave_FC2_intron_P769$P769D25_ave >=2 & TE_cpm_anno_ave_FC2_intron_P769$TE_Class == "DNA",]) ))
colnames(bothTE_FC2_DNA) <-c("Sample","anno","Count")
bothTE_FC2_DNA<-merge(meta1[,c(1:6)],bothTE_FC2_DNA, by="Sample")
bothTE_FC2_DNA$Day <-factor(bothTE_FC2_DNA$Day, levels = c("5","12","14","15","22","23","25"))
bothTE_FC2_DNA$Cell <-factor(bothTE_FC2_DNA$Cell, levels = c("ACHN","CAKI2","P769","O786","CAKI","A704","A498"))
bothTE_FC2_DNA$CellType <- paste(bothTE_FC2_DNA$Cell,bothTE_FC2_DNA$Type,sep="_")
bothTE_FC2_DNA$CellType <-factor(bothTE_FC2_DNA$CellType, levels = c("ACHN_WT","O786_EV","O786_KO","CAKI_Mut","A498_Mut","CAKI2_WT","A704_Mut","P769_WT"))
bothTE_FC2_DNA<-bothTE_FC2_DNA[order(bothTE_FC2_DNA$CellType),]

pc <- ggplot(bothTE_FC2_DNA[bothTE_FC2_DNA$anno == "Intergenic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic DNA >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pc <-pc+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")
pc1 <- ggplot(bothTE_FC2_DNA[bothTE_FC2_DNA$anno == "Intronic",], aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intron DNA >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+scale_fill_manual(values = c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+scale_color_manual(values= c(mypalette[c(5,9,10,7,1)],"#5cb85c","#c89f73","#999999"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pc1 <-pc1+ facet_wrap( ~ Cell, scales="free",ncol=7)+coord_cartesian(ylim = c(0,700))+ theme(legend.position="none")

multiplot(p,pb,pa,pc,p1,pb1,pa1,pc1,cols=2) ## save 10 x16




#Intergenic TE Bar plot of Number of TEs >2FC (Fig. 5A) FOR 786-O
intergenicTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23"),c(nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intergenic_O786[TE_cpm_anno_ave_FC2_intergenic_O786$SK0D23_ave >=2,])))
colnames(intergenicTE_FC2) <-c("Sample","Count")
intergenicTE_FC2<-merge(meta[,c(1:6)],intergenicTE_FC2, by="Sample")
intergenicTE_FC2$Day <-factor(intergenicTE_FC2$Day, levels = c("5","15","23"))
intergenicTE_FC2$Cell <-factor(intergenicTE_FC2$Cell, levels = c("O786"))
intergenicTE_FC2$CellType <- paste(intergenicTE_FC2$Cell,intergenicTE_FC2$Type,sep="_")

p <- ggplot(intergenicTE_FC2, aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Intergenic TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+ scale_y_continuous(limits = c(0, max(intergenicTE_FC2$Count)))+scale_fill_manual(values = mypalette[c(9,10)])+scale_color_manual(values= mypalette[c(9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))


#Intron TE Bar plot of Number of TEs >2FC (Fig. 5A)
intronTE_FC2<-data.frame(c("OEVD5","OEVD15","OEVD23","SK0D5","SK0D15","SK0D23"),c(nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$OEVD23_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D5_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D15_ave >=2,]),nrow(TE_cpm_anno_ave_FC2_intron_O786[TE_cpm_anno_ave_FC2_intron_O786$SK0D23_ave >=2,])))
colnames(intronTE_FC2) <-c("Sample","Count")
intronTE_FC2<-merge(meta[,c(1:6)],intronTE_FC2, by="Sample")
intronTE_FC2$Day <-factor(intronTE_FC2$Day, levels = c("5","15","23"))
intronTE_FC2$Cell <-factor(intronTE_FC2$Cell, levels = c("O786"))
intronTE_FC2$CellType <- paste(intronTE_FC2$Cell,intronTE_FC2$Type,sep="_")

p <- ggplot(intronTE_FC2, aes(x=Day, y=Count, fill=CellType)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("intron TE >2 FC count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Days", y = "TE count")+ scale_y_continuous(limits = c(0, max(intronTE_FC2$Count)))+scale_fill_manual(values = mypalette[c(9,10)])+scale_color_manual(values= mypalette[c(9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))




###    Make Z-score of TEs >2FC, 786-O only (Fig. 3E)   ###############
O786_aveCPM_with_FC <- merge(TE_cpm_anno_ave_O786,TE_cpm_anno_ave_FC2_O786, by= c(colnames(TE_cpm_anno_ave)[1:10]))

O786_aveCPM_with_FC_Zscore<- (log2(O786_aveCPM_with_FC[,c(11:18)]+0.01)-rowMeans(log2(O786_aveCPM_with_FC[,c(11:18)]+0.01)))/(rowSds(as.matrix(log2(O786_aveCPM_with_FC[,c(11:18)]+0.01))))[row(log2(O786_aveCPM_with_FC[,c(11:18)]+0.01))]
name <- rownames(O786_aveCPM_with_FC_Zscore)
df.OG2 <- data.matrix(O786_aveCPM_with_FC_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(O786_aveCPM_with_FC_Zscore)
htname_786O =  Heatmap(df.OG2, column_title = "Timepoint",name= "O786 TEs\n(Zscore)",col = colorRamp2(c(-1.5,-0.8,0, 0.8, 1.5), c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000")), 
                       cluster_rows = T, cluster_columns = FALSE,show_row_names = F)
htname_786O


df.OG2 <- data.matrix(O786_aveCPM_with_FC[row_order(htname_786O),c(11:18)])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(O786_aveCPM_with_FC[,c(11:18)])
Heatmap(df.OG2, column_title = "Timepoint",name= "O786 TEs\n(CPM)",col = colorRamp2(c(0,0.25,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5), c(rev(mypalette5))), 
        cluster_rows = F, cluster_columns = FALSE,show_row_names = F) ## (Fig. S4D (row dendrogram from Fig. 3E))


#######Identify Activated TEs (>1cpm) in O786 EV and KO COMPARISION###############
TE_cpm_anno_ave_O786 <- TE_cpm_anno_ave[,c(1:10,23:30)]
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
TE_cpm_anno_ave_O786_alldata_FC2_allKOup <- TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 > 2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >1,] #344 intronic & 207 intergenic


## Make plot comparing the CPM expression of induced INTRON TEs (>2FC) in EV vs KO (Fig. 5B) ####
TE_cpm_anno_ave_O786_alldata_FC2_intron_allTE <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$OEVD5_ave_FC > 2 | TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2,] 
TE_cpm_anno_ave_O786_alldata_FC2_intron_allKOup <- TE_cpm_anno_ave_O786_alldata_FC2_intron[TE_cpm_anno_ave_O786_alldata_FC2_intron$EVvKO_D5 > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave_FC > 2 & TE_cpm_anno_ave_O786_alldata_FC2_intron$SK0D5_ave >1,] #344 intronic 

p<-ggplot(TE_cpm_anno_ave_O786_alldata_FC2_intron_allTE, aes(x=log2(OEVD5_ave) , y=log2(SK0D5_ave)) ) +
  geom_point(color="grey") + 
  theme_bw()
p<-p + geom_point(data = TE_cpm_anno_ave_O786_alldata_FC2_intron_allKOup, aes(x=log2(OEVD5_ave) , y=log2(SK0D5_ave)), color = "#ffcc00")	
p+geom_abline(slope=1, intercept=0)+ggtitle("EV intron TE vs KO intron TE (KO D5 > 2 FC)")+coord_cartesian(xlim=c(-6,6),ylim=c(-6,6))



######### Enrichment of TEs across 786-O EV and KO cells (Fig. S8D) ##############
TE_rmsk <- read.table("GRCh38_GENCODE_rmsk_TE_GTF_parsed.txt",sep = "\t",header=T,stringsAsFactors=F)
TE_annot_full <- read.delim("GRCh38_GENCODE_rmsk_TE_GTF_parsed.V37annotate.bed", sep ="\t", header =T, stringsAsFactors = F)
TE_annot_full <- TE_annot_full[,c(2,3,4,8)]
TE_annot_full$Anno2 <- apply(TE_annot_full, 1, function(x) {unlist(strsplit(x[4], " \\("))[1]})
TE_annot_full$Start <- TE_annot_full$Start-1 #make into 0-based
TE_rmsk_anno<-merge(TE_rmsk,TE_annot_full[,c(1,2,3,5)],by=c("Chr","Start","End"),x.all=T)

#### INTRON+INTERGENIC TE ENRICHMENT, don't include SVA since low numbers
TotalTE <- nrow(TE_rmsk_anno[(TE_rmsk_anno$TEclass=="DNA" | TE_rmsk_anno$TEclass=="LINE" |TE_rmsk_anno$TEclass=="LTR" | TE_rmsk_anno$TEclass=="SINE") & (TE_rmsk_anno$Anno2=="Intergenic" | TE_rmsk_anno$Anno2=="intron"),])
totalDNA <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="DNA" & TE_rmsk_anno$TEclass !="Retroposon" & (TE_rmsk_anno$Anno2=="Intergenic" | TE_rmsk_anno$Anno2=="intron"),])
totalLINE <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="LINE"& TE_rmsk_anno$TEclass !="Retroposon" & (TE_rmsk_anno$Anno2=="Intergenic" | TE_rmsk_anno$Anno2=="intron"),])
totalSINE <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="SINE"& TE_rmsk_anno$TEclass !="Retroposon" & (TE_rmsk_anno$Anno2=="Intergenic" | TE_rmsk_anno$Anno2=="intron"),])
totalLTR <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="LTR"& TE_rmsk_anno$TEclass !="Retroposon" & (TE_rmsk_anno$Anno2=="Intergenic" | TE_rmsk_anno$Anno2=="intron"),])

#INTRON+INTERGENICTE Class Enrichment (# TEclass active/ #total TE active)/(# total TE class/# total TE)
Class_KO5 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[(TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786_FC2$Anno2 == "intron") & TE_cpm_anno_ave_O786_FC2$SK0D5_ave >=2,]$TE_Class))
Class_KO5$total <- sum(Class_KO5$Freq)
Class_KO5$sample <- "KO_D5"
Class_EV5 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[(TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786_FC2$Anno2 == "intron") & TE_cpm_anno_ave_O786_FC2$OEVD5_ave >=2,]$TE_Class))
Class_EV5$total <- sum(Class_EV5$Freq)
Class_EV5$sample <- "EV_D5"
Class_KO15 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[(TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786_FC2$Anno2 == "intron") & TE_cpm_anno_ave_O786_FC2$SK0D15_ave >=2,]$TE_Class))
Class_KO15$total <- sum(Class_KO15$Freq)
Class_KO15$sample <- "KO_D15"
Class_EV15 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[(TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786_FC2$Anno2 == "intron") & TE_cpm_anno_ave_O786_FC2$OEVD15_ave >=2,]$TE_Class))
Class_EV15$total <- sum(Class_EV15$Freq)
Class_EV15$sample <- "EV_D15"
Class_KO23 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[(TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786_FC2$Anno2 == "intron") & TE_cpm_anno_ave_O786_FC2$SK0D23_ave >=2,]$TE_Class))
Class_KO23$total <- sum(Class_KO23$Freq)
Class_KO23$sample <- "KO_D23"
Class_EV23 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[(TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" | TE_cpm_anno_ave_O786_FC2$Anno2 == "intron") & TE_cpm_anno_ave_O786_FC2$OEVD23_ave >=2,]$TE_Class))
Class_EV23$total <- sum(Class_EV23$Freq)
Class_EV23$sample <- "EV_D23"


Class <- rbind(Class_EV5,Class_KO5,Class_EV15,Class_KO15,Class_EV23,Class_KO23)
Class <- Class[Class$Var1 != "Retroposon",]
Class$Enrich <- apply(Class,1,function(x) {
  if (x[1]=="DNA") {((as.numeric(x[2])/as.numeric(x[3]))/(totalDNA/TotalTE))} 
  else if (x[1]=="LINE") {((as.numeric(x[2])/as.numeric(x[3]))/(totalLINE/TotalTE))}
  else if (x[1]=="SINE") {((as.numeric(x[2])/as.numeric(x[3]))/(totalSINE/TotalTE))}
  else if (x[1]=="LTR") {((as.numeric(x[2])/as.numeric(x[3]))/(totalLTR/TotalTE))}
})
Class$log2Enrich <- log2(Class$Enrich)
Class$log2Enrich <- apply(Class,1,function(x) {if (as.numeric(x[6]) < -1) {-1} else if (as.numeric(x[6])>1) {1} else x[6]})#nothing higher than 2 log2 enrichment
Class$sample <- factor(Class$sample, levels = c("EV_D5","EV_D15","EV_D23","KO_D5","KO_D15","KO_D23"))
ggplot(Class, aes(sample,Var1, fill=as.numeric(log2Enrich)))+geom_tile()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black')) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_distiller(limits=c(-1,1), palette = "RdBu")



#Intergenic TE Class Enrichment, don't include SVA since low numbers
TotalTE <- nrow(TE_rmsk_anno[(TE_rmsk_anno$TEclass=="DNA" | TE_rmsk_anno$TEclass=="LINE" |TE_rmsk_anno$TEclass=="LTR" | TE_rmsk_anno$TEclass=="SINE") & TE_rmsk_anno$Anno2=="Intergenic",])
totalDNA <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="DNA" & TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="Intergenic",])
totalLINE <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="LINE"& TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="Intergenic",])
totalSINE <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="SINE"& TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="Intergenic",])
totalLTR <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="LTR"& TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="Intergenic",])

#Intergenic TE Class Enrichment (# TEclass active/ #total TE active)/(# total TE class/# total TE)
Class_KO5 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_FC2$SK0D5_ave >=2,]$TE_Class))
Class_KO5$total <- sum(Class_KO5$Freq)
Class_KO5$sample <- "KO_D5"
Class_EV5 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_FC2$OEVD5_ave >=2,]$TE_Class))
Class_EV5$total <- sum(Class_EV5$Freq)
Class_EV5$sample <- "EV_D5"
Class_KO15 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_FC2$SK0D15_ave >=2,]$TE_Class))
Class_KO15$total <- sum(Class_KO15$Freq)
Class_KO15$sample <- "KO_D15"
Class_EV15 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_FC2$OEVD15_ave >=2,]$TE_Class))
Class_EV15$total <- sum(Class_EV15$Freq)
Class_EV15$sample <- "EV_D15"
Class_KO23 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_FC2$SK0D23_ave >=2,]$TE_Class))
Class_KO23$total <- sum(Class_KO23$Freq)
Class_KO23$sample <- "KO_D23"
Class_EV23 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_FC2$OEVD23_ave >=2,]$TE_Class))
Class_EV23$total <- sum(Class_EV23$Freq)
Class_EV23$sample <- "EV_D23"


Class <- rbind(Class_EV5,Class_KO5,Class_EV15,Class_KO15,Class_EV23,Class_KO23)
Class <- Class[Class$Var1 != "Retroposon",]
Class$Enrich <- apply(Class,1,function(x) {
  if (x[1]=="DNA") {((as.numeric(x[2])/as.numeric(x[3]))/(totalDNA/TotalTE))} 
  else if (x[1]=="LINE") {((as.numeric(x[2])/as.numeric(x[3]))/(totalLINE/TotalTE))}
  else if (x[1]=="SINE") {((as.numeric(x[2])/as.numeric(x[3]))/(totalSINE/TotalTE))}
  else if (x[1]=="LTR") {((as.numeric(x[2])/as.numeric(x[3]))/(totalLTR/TotalTE))}
})
Class$log2Enrich <- log2(Class$Enrich)
Class$log2Enrich <- apply(Class,1,function(x) {if (as.numeric(x[6]) < -1) {-1} else if (as.numeric(x[6])>1) {1} else x[6]})#nothing higher than 2 log2 enrichment
Class$sample <- factor(Class$sample, levels = c("EV_D5","EV_D15","EV_D23","KO_D5","KO_D15","KO_D23"))
ggplot(Class, aes(sample,Var1, fill=as.numeric(log2Enrich)))+geom_tile()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black')) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_distiller(limits=c(-1,1), palette = "RdBu")


#intron TE Class Enrichment (# TEclass active/ #total TE active)/(# total TE class/# total TE)
#don't include SVA since low numbers
TotalTE <- nrow(TE_rmsk_anno[(TE_rmsk_anno$TEclass=="DNA" | TE_rmsk_anno$TEclass=="LINE" |TE_rmsk_anno$TEclass=="LTR" | TE_rmsk_anno$TEclass=="SINE") & TE_rmsk_anno$Anno2=="intron",])
totalDNA <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="DNA" & TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="intron",])
totalLINE <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="LINE"& TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="intron",])
totalSINE <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="SINE"& TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="intron",])
totalLTR <- nrow(TE_rmsk_anno[TE_rmsk_anno$TEclass=="LTR"& TE_rmsk_anno$TEclass !="Retroposon" & TE_rmsk_anno$Anno2=="intron",])

#intron TE Class Enrichment (# TEclass active/ #total TE active)/(# total TE class/# total TE)
Class_KO5 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_FC2$SK0D5_ave >=2,]$TE_Class))
Class_KO5$total <- sum(Class_KO5$Freq)
Class_KO5$sample <- "KO_D5"
Class_EV5 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_FC2$OEVD5_ave >=2,]$TE_Class))
Class_EV5$total <- sum(Class_EV5$Freq)
Class_EV5$sample <- "EV_D5"
Class_KO15 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_FC2$SK0D15_ave >=2,]$TE_Class))
Class_KO15$total <- sum(Class_KO15$Freq)
Class_KO15$sample <- "KO_D15"
Class_EV15 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_FC2$OEVD15_ave >=2,]$TE_Class))
Class_EV15$total <- sum(Class_EV15$Freq)
Class_EV15$sample <- "EV_D15"
Class_KO23 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_FC2$SK0D23_ave >=2,]$TE_Class))
Class_KO23$total <- sum(Class_KO23$Freq)
Class_KO23$sample <- "KO_D23"
Class_EV23 <- data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_FC2$OEVD23_ave >=2,]$TE_Class))
Class_EV23$total <- sum(Class_EV23$Freq)
Class_EV23$sample <- "EV_D23"


Class <- rbind(Class_EV5,Class_KO5,Class_EV15,Class_KO15,Class_EV23,Class_KO23)
Class <- Class[Class$Var1 != "Retroposon",]
Class$Enrich <- apply(Class,1,function(x) {
  if (x[1]=="DNA") {((as.numeric(x[2])/as.numeric(x[3]))/(totalDNA/TotalTE))} 
  else if (x[1]=="LINE") {((as.numeric(x[2])/as.numeric(x[3]))/(totalLINE/TotalTE))}
  else if (x[1]=="SINE") {((as.numeric(x[2])/as.numeric(x[3]))/(totalSINE/TotalTE))}
  else if (x[1]=="LTR") {((as.numeric(x[2])/as.numeric(x[3]))/(totalLTR/TotalTE))}
})
Class$log2Enrich <- log2(Class$Enrich)
Class$log2Enrich <- apply(Class,1,function(x) {if (as.numeric(x[6]) < -1) {-1} else if (as.numeric(x[6])>1) {1} else x[6]})#nothing higher than 2 log2 enrichment
Class$sample <- factor(Class$sample, levels = c("EV_D5","EV_D15","EV_D23","KO_D5","KO_D15","KO_D23"))
ggplot(Class, aes(sample,Var1, fill=as.numeric(log2Enrich)))+geom_tile()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black')) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_distiller(limits=c(-1,1), palette = "RdBu")





##make bargraph showing distribution of TE class and Genomic annotation (Fig. S8C)
a<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$OEVD5_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]$TE_Class))
c<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]$TE_Class))
a1<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$OEVD5_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]$TE_Class))
c1<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]$TE_Class))


TE_class_freq <- rbind(a,c,a1,c1)
TE_class_freq$Day <- c(rep("Day5",16))
TE_class_freq$comp <- c("EV","EV","EV","EV","KO","KO","KO","KO","EV","EV","EV","EV","KO","KO","KO","KO")
TE_class_freq$comp <- factor(TE_class_freq$comp, levels = c("EV","KO"))
TE_class_freq$anno <- c(rep("intron",8),rep("Intergenic",8))
TE_class_freq$TE_anno <- apply(TE_class_freq, 1, function(x) {paste(as.character(x[1]),x[5],sep="_")})
TE_class_freq$TE_anno <-factor(TE_class_freq$TE_anno, levels = c("LINE_intron","LINE_Intergenic","LTR_intron","LTR_Intergenic","SINE_intron","SINE_Intergenic","DNA_intron","DNA_Intergenic"))

p<-ggplot(TE_class_freq, aes(x=comp, y=Freq, fill=TE_anno)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("TE annotation (FC >2): Day 0 vs Day 5")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Upregulated TEs", y = "Fraction of total")+scale_fill_manual(values = c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+scale_color_manual(values=c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p
pp<-ggplot(TE_class_freq, aes(x=comp, y=Freq, fill=TE_anno)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("TE annotation (FC >2): Day 0 vs Day 5")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Upregulated TEs", y = "Number of TEs")+scale_fill_manual(values = c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+scale_color_manual(values=c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp


a<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$OEVD15_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]$TE_Class))
c<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$SK0D15_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]$TE_Class))
a1<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$OEVD15_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]$TE_Class))
c1<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$SK0D15_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]$TE_Class))

TE_class_freq <- rbind(a,c,a1,c1)
TE_class_freq$Day <- c(rep("Day15",16))
TE_class_freq$comp <- c("EV","EV","EV","EV","KO","KO","KO","KO","EV","EV","EV","EV","KO","KO","KO","KO")
TE_class_freq$comp <- factor(TE_class_freq$comp, levels = c("EV","KO"))
TE_class_freq$anno <- c(rep("intron",8),rep("Intergenic",8))
TE_class_freq$TE_anno <- apply(TE_class_freq, 1, function(x) {paste(as.character(x[1]),x[5],sep="_")})
TE_class_freq$TE_anno <-factor(TE_class_freq$TE_anno, levels = c("LINE_intron","LINE_Intergenic","LTR_intron","LTR_Intergenic","SINE_intron","SINE_Intergenic","DNA_intron","DNA_Intergenic"))

p2<-ggplot(TE_class_freq, aes(x=comp, y=Freq, fill=TE_anno)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("TE annotation (FC >2): Day 0 vs Day 15")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Upregulated TEs", y = "Fraction of total")+scale_fill_manual(values = c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+scale_color_manual(values=c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p2
pp2<-ggplot(TE_class_freq, aes(x=comp, y=Freq, fill=TE_anno)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("TE annotation (FC >2): Day 0 vs Day 15")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Upregulated TEs", y = "Number of TEs")+scale_fill_manual(values = c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+scale_color_manual(values=c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp2



a<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$OEVD23_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]$TE_Class))
c<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$SK0D23_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "intron",]$TE_Class))
a1<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$OEVD23_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]$TE_Class))
c1<-data.frame(table(TE_cpm_anno_ave_O786_FC2[TE_cpm_anno_ave_O786_FC2$SK0D23_ave_FC >=2 & TE_cpm_anno_ave_O786_FC2$Anno2 == "Intergenic",]$TE_Class))

TE_class_freq <- rbind(a,c,a1,c1)
TE_class_freq$Day <- c(rep("Day23",16))
TE_class_freq$comp <- c("EV","EV","EV","EV","KO","KO","KO","KO","EV","EV","EV","EV","KO","KO","KO","KO")
TE_class_freq$comp <- factor(TE_class_freq$comp, levels = c("EV","KO"))
TE_class_freq$anno <- c(rep("intron",8),rep("Intergenic",8))
TE_class_freq$TE_anno <- apply(TE_class_freq, 1, function(x) {paste(as.character(x[1]),x[5],sep="_")})
TE_class_freq$TE_anno <-factor(TE_class_freq$TE_anno, levels = c("LINE_intron","LINE_Intergenic","LTR_intron","LTR_Intergenic","SINE_intron","SINE_Intergenic","DNA_intron","DNA_Intergenic"))

p1<-ggplot(TE_class_freq, aes(x=comp, y=Freq, fill=TE_anno)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("TE annotation (FC >2): Day 0 vs Day 23")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Upregulated TEs", y = "Fraction of total")+scale_fill_manual(values = c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+scale_color_manual(values=c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p1
pp1<-ggplot(TE_class_freq, aes(x=comp, y=Freq, fill=TE_anno)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("TE annotation (FC >2): Day 0 vs Day 23")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Upregulated TEs", y = "Number of TEs")+scale_fill_manual(values = c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+scale_color_manual(values=c("#ffa700","#ffeead","#008744","#88d8b0","#d62d20","#ff6f69","#428bca","#b3dcff","#428bca","#b3dcff"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp1


multiplot(pp,p,pp2,p2,pp1,p1,cols =3)





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

TE_denovo_allKO <- TE_denovo2[TE_denovo2$SK0D5_ave_FC >= 2 & TE_denovo2$EVvKO_D5 >= 2 & TE_denovo2$SK0D5_ave >=1,]  #all TEs induced and up-regulated in KO (>2FC than EV)


#iterate with different filters to get quantify what % of TEs overlap denovo or canonical exons transcript (Fig. S5A)
length(unique(TE_denovo_allKO[TE_denovo_allKO$Anno2 == "Intergenic",1])) #207
length(unique(TE_denovo_allKO[TE_denovo_allKO$Anno2 == "intron",1])) #344

TE_denovo_allKO_exon <- TE_denovo_allKO[(TE_denovo_allKO$Anno2 == "intron" | TE_denovo_allKO$Anno2 == "Intergenic") & TE_denovo_allKO$anno == "exon",]
TE_denovo_allKO_exon$qry_id <- apply(TE_denovo_allKO_exon ,1,function(x) {gsub(";","",unlist(strsplit(x[39], " "))[2])})
TE_denovo_allKO_exon$exon <- apply(TE_denovo_allKO_exon ,1,function(x) {gsub(";","",unlist(strsplit(x[39], " "))[6])})
TE_denovo_allKO_exon <- TE_denovo_allKO_exon[,c(1:34,37,41,42)]

TE_denovo_allKO_exon <- merge(TE_denovo_allKO_exon,tmap_denovo[,c("class_code","qry_gene_id","qry_id","num_exons")], by = "qry_id")
TE_denovo_allKO_exon <-TE_denovo_allKO_exon[!duplicated(TE_denovo_allKO_exon),]
TE_denovo_allKO_exon <- TE_denovo_allKO_exon[which(TE_denovo_allKO_exon$Strand == TE_denovo_allKO_exon$strand),] 
length(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "Intergenic",]$Geneid)) #143 or 69.1%
length(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron",]$Geneid)) #196 or 57.0%

write.table(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron",c(3:5)]),"Misspliced_intronTE_KOvEV_D5_2FC_location.bed", sep = "\t", col.names =F, row.names = F, quote=F)
write.table(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron" & TE_denovo_allKO_exon$TE_Class == "LINE",c(3:5)]),"Misspliced_intronLINE_KOvEV_D5_2FC_location.bed", sep = "\t", col.names =F, row.names = F, quote=F)
write.table(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron" & TE_denovo_allKO_exon$TE_Class == "LTR",c(3:5)]),"Misspliced_intronLTR_KOvEV_D5_2FC_location.bed", sep = "\t", col.names =F, row.names = F, quote=F)
write.table(unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron" & TE_denovo_allKO_exon$TE_Class == "SINE",c(3:5)]),"Misspliced_intronSINE_KOvEV_D5_2FC_location.bed", sep = "\t", col.names =F, row.names = F, quote=F)



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

write.table(unique(t_filter_longest_exon_with_info_only_exons[t_filter_longest_exon_with_info_only_exons$ExonName %in% t_filter_longest_exon_with_info_only_exons2[t_filter_longest_exon_with_info_only_exons2$class_code != "=",]$ExonName,c(1:3)]),"Misspliced_LongestExon_wTE_location.bed",sep= "\t", col.names = F, row.names = F,quote =F)


##SETD2-KO up-reg: Combine transcript class - Group 0 (=), Group 1 exon discrepancy (j,k,o), intron retention Group 2 (m,n), Group 3 novel transcript (i,x,y). Priority: Group1 > Group 2 >Group 3 > Group 0 
= =,j,k   =,k =,k,j     i i,k,j     j   j,i   j,k j,n,=     k k,=,j   k,j k,j,=     m     n     o     x     y 
14     2     1     1    29     1    29     1     3     1    18     1     3     1     5     7     2     7     1 

group 0: 14
group 1: 2+1+1+1+29+1+3+1+18+1+3+1+2 = 64
group 2: 7+5 = 12
group 3: 29+7+1 = 37

#make figure for frequency of different groups of exons (Fig 5C)
denovo <- data.frame(c(14,64,12,37),c(rep("allKO",4)),c(rep(c("g0","g1","g2","g3"),1)))
colnames(denovo) <- c("Count","Type","Group")
figurecolor <- c("#44bec7","#fa3c4c","#ffc300","#d696bb")

p1 <- ggplot(denovo, aes(x=Type, y=Count, fill=Group)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("De novo Exon type")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Type", y = "TE count")+scale_fill_manual(values = mypalette[c(1,5,7,3)])+scale_color_manual(values= figurecolor)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))


############################ HISTONE CHIP-seq ANALYSIS ############################
#R version 4.0.3 
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
mypalette7 <- brewer.pal(11,"Spectral")

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


##plot FRiP (Fig. S5B)
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



##plot PCA from Diffbind (Fig. S5C)
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

#Annotate all peaks using Homer (Fig. S5D)
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

#deeptools heatmap and linegraph of histone signal across gene body (Fig. 4B) . H3K9me3 example.
computeMatrix scale-regions --missingDataAsZero -bs 100 -p 12 -R /genomes/hg38/gencode.v37.AutosomalTranscriptsOnly.forDeepTools.gtf -S EV0K9_Combined_RPGC.bw EV5K9_Combined_RPGC.bw EV23K9_Combined_RPGC.bw SK0K9_Combined_RPGC.bw SK5K9_Combined_RPGC.bw SK23K9_Combined_RPGC.bw -m 5000 -b 3000 -a 3000 -out H3K9me3_GeneBody_O786_RPGC_combined.tab.gz 
plotHeatmap  --heatmapWidth 6 -m H3K9me3_GeneBody_O786_RPGC_combined.tab.gz --perGroup -out H3K9me3_GeneBody_O786_RPGC.pdf --missingDataColor "grey" --colorMap Greens --heatmapHeight 10




#### PERFORM TE centric analysis to quantify histone changes ###########
setwd("../TE_centric")
write.table(TE_cpm_anno_ave_O786_alldata_FC2[,c(2:4,1,5:30)],"SETD2_TE_averageTPM_Only786_FoldChangeGreaterThan2_V37.txt",sep = "\t", col.names = F, row.names =F, quote =F)
write.table(TE_cpm_anno_ave_O786_alldata_FC2[,c(2:4,1)],"TE_FC2_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F)
write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intergenic_TE_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F)#207
write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_alldata_FC2$OEVD5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$OEVD5_ave >=1 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave <1,c(2:4,1)],"Intronic_TE_EVonly_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #344


write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$TE_Class == "SINE" & TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intronic_SINE_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #98
write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$TE_Class == "LTR" & TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intronic_LTR_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #81
write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$TE_Class == "LINE" & TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "intron" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intronic_LINE_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #136

write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$TE_Class == "SINE" & TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intergenic_SINE_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #60
write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$TE_Class == "LTR" & TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intergenic_LTR_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #51
write.table(TE_cpm_anno_ave_O786_alldata_FC2[TE_cpm_anno_ave_O786_alldata_FC2$TE_Class == "LINE" & TE_cpm_anno_ave_O786_alldata_FC2$Anno2 == "Intergenic" & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave_FC >=2 & TE_cpm_anno_ave_O786_alldata_FC2$EVvKO_D5 >=2 & TE_cpm_anno_ave_O786_alldata_FC2$SK0D5_ave >=1,c(2:4,1)],"Intergenic_LINE_KOvEV_D5_FC2_1CPM_O786_location.bed",sep = "\t", col.names = F, row.names =F, quote =F) #83




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

###all time points comparison but order on D5 TE expression (Zscore)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2<- cbind(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(1:3)],log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(11:18)]+0.01))
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5 <- TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$SK0D5_ave

TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2[,c(4:11)] -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2[,c(4:11)]))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2[,c(4:11)])))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2[,c(4:11)])]

####### ALL TE expression Zscore heatmap (Fig. 4C)
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore) 
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),])
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
colnames(df.OG2) <- colnames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),])
d<-pheatmap(df.OG2,color = colorRampPalette(c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000"))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

#CPM heatmap (of Z-score, Fig. S6C)
Heatmap(data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_TE2$diffEV5_KO5),c(11:18)]), column_title = "Timepoint",name= "all TEs\n(CPM)",col = colorRamp2(c(0,0.25,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5), c(rev(mypalette7))), 
        cluster_rows = F, cluster_columns = FALSE,show_row_names = F)

#TE expression CPM boxplot (Fig. 4C)
TE_aveCPM_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(11:18)])
colnames(TE_aveCPM_m) <- c("Sample","CPM")
TE_aveCPM_m$Cell <- apply(TE_aveCPM_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
TE_aveCPM_m$Sample <- factor(TE_aveCPM_m$Sample, levels = c("OEVD0_ave","OEVD5_ave","OEVD15_ave","OEVD23_ave","SK0D0_ave","SK0D5_ave","SK0D15_ave","SK0D23_ave"))
TE_aveCPM_m$Cell <- factor(TE_aveCPM_m$Cell, levels = c("OEVD0","OEVD5","OEVD15","OEVD23","SK0D0","SK0D5","SK0D15","SK0D23"))
ppppp<-ggplot(TE_aveCPM_m, aes(x=Sample, y=CPM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("TE Expression Distribution (CPM)")+ coord_cartesian(ylim=c(0, 3))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "CPM")+scale_fill_manual(values = c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+scale_color_manual(values=c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot



### histone Zscore heatmaps by TE Zscore (Fig. 4C & S6C) ###
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
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K9)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K36)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K4)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_K27)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_Ac)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


##Boxplot
mypalette3 <- c("#2a4d69", "#4b86b4","#adcbe3","#9c0000","#fd0000","#ff9797")

K36me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(115:120)])
colnames(K36me3_rpgc_m) <- c("Sample","RPGC")
K36me3_rpgc_m$Cell <- apply(K36me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K36me3_rpgc_m$Sample <- factor(K36me3_rpgc_m$Sample, levels = c("EV0K36_ave","EV5K36_ave","EV23K36_ave","SK0K36_ave","SK5K36_ave","SK23K36_ave"))
K36me3_rpgc_m$Cell <- factor(K36me3_rpgc_m$Cell, levels = c("EV0K36","EV5K36","EV23K36","SK0K36","SK5K36","SK23K36"))
ppppp<-ggplot(K36me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K36me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 26))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

K9me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(109:114)])
colnames(K9me3_rpgc_m) <- c("Sample","RPGC")
K9me3_rpgc_m$Cell <- apply(K9me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K9me3_rpgc_m$Sample <- factor(K9me3_rpgc_m$Sample, levels = c("EV0K9_ave","EV5K9_ave","EV23K9_ave","SK0K9_ave","SK5K9_ave","SK23K9_ave"))
K9me3_rpgc_m$Cell <- factor(K9me3_rpgc_m$Cell, levels = c("EV0K9","EV5K9","EV23K9","SK0K9","SK5K9","SK23K9"))
pppp<-ggplot(K9me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K9me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pppp##boxplot


K27me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(103:108)])
colnames(K27me3_rpgc_m) <- c("Sample","RPGC")
K27me3_rpgc_m$Cell <- apply(K27me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K27me3_rpgc_m$Sample <- factor(K27me3_rpgc_m$Sample, levels = c("EV0K27_ave","EV5K27_ave","EV23K27_ave","SK0K27_ave","SK5K27_ave","SK23K27_ave"))
K27me3_rpgc_m$Cell <- factor(K27me3_rpgc_m$Cell, levels = c("EV0K27","EV5K27","EV23K27","SK0K27","SK5K27","SK23K27"))
ppp<-ggplot(K27me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K27me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 4))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppp##boxplot

K4me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(91:96)])
colnames(K4me3_rpgc_m) <- c("Sample","RPKM")
K4me3_rpgc_m$Cell <- apply(K4me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K4me3_rpgc_m$Sample <- factor(K4me3_rpgc_m$Sample, levels = c("EV0K4_ave","EV5K4_ave","EV23K4_ave","SK0K4_ave","SK5K4_ave","SK23K4_ave"))
K4me3_rpgc_m$Cell <- factor(K4me3_rpgc_m$Cell, levels = c("EV0K4","EV5K4","EV23K4","SK0K4","SK5K4","SK23K4"))
pp<-ggplot(K4me3_rpgc_m, aes(x=Sample, y=RPKM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K4me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 22))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPKM")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pp##boxplot


K27ac_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(97:102)])
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

####### WILCOXON RANK-SUM TEST for statistical comparison of TE expression and histone PTM distribution #########
m_TE<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(11:18)])
pairwise.wilcox.test(m_TE$value,m_TE$variable)

m_K4<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(91:96)])
pairwise.wilcox.test(m_K4$value,m_K4$variable)

m_Ac<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(97:102)])
pairwise.wilcox.test(m_Ac$value,m_Ac$variable)

m_K27<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(103:108)])
pairwise.wilcox.test(m_K27$value,m_K27$variable)

m_K9<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(109:114)])
pairwise.wilcox.test(m_K9$value,m_K9$variable)

m_K36<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[,c(115:120)])
pairwise.wilcox.test(m_K36$value,m_K36$variable)


#######Intergenic TE  comparison but order on D5 TE expression (Fig. S9B)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE <-TE_cpm_anno_aveCPM_O786_FC2_both_histone[TE_cpm_anno_aveCPM_O786_FC2_both_histone$Anno2 == "Intergenic",]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2<- cbind(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(1:3)],log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(11:18)]+0.01))
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5 <- TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$SK0D5_ave
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5_CPM <- TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$SK0D5_ave

TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2[,c(4:11)] -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2[,c(4:11)]))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2[,c(4:11)])))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2[,c(4:11)])]

# Intergenic TE expression Zscore heatmap 
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),])
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
colnames(df.OG2) <- colnames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),])
d<-pheatmap(df.OG2,color = colorRampPalette(c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000"))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

#CPM heatmap
Heatmap(data.matrix(2^TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),c(4:11)]), column_title = "Timepoint",name= "intergenic TEs\n(CPM)",col = colorRamp2(c(0,0.25,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5), c(rev(mypalette7))), 
        cluster_rows = F, cluster_columns = FALSE,show_row_names = F)

#intergenic TE expression CPM boxplot 
TE_aveCPM_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[TE_cpm_anno_aveCPM_O786_FC2_both_histone$Anno2 == "Intergenic",c(11:18)])
colnames(TE_aveCPM_m) <- c("Sample","CPM")
TE_aveCPM_m$Cell <- apply(TE_aveCPM_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
TE_aveCPM_m$Sample <- factor(TE_aveCPM_m$Sample, levels = c("OEVD0_ave","OEVD5_ave","OEVD15_ave","OEVD23_ave","SK0D0_ave","SK0D5_ave","SK0D15_ave","SK0D23_ave"))
TE_aveCPM_m$Cell <- factor(TE_aveCPM_m$Cell, levels = c("OEVD0","OEVD5","OEVD15","OEVD23","SK0D0","SK0D5","SK0D15","SK0D23"))
ppppp<-ggplot(TE_aveCPM_m, aes(x=Sample, y=CPM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("Intergenic TE Expression Distribution (CPM)")+ coord_cartesian(ylim=c(0, 3))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "CPM")+scale_fill_manual(values = c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+scale_color_manual(values=c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

### histone Zscore heatmaps by intergenic TE Zscore ###
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),c(91:96)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),c(97:102)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),c(103:108)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),c(109:114)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intergenicTE2$diffEV5_KO5),c(115:120)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36)]

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K9)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K36)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K4)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_K27)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE_log2_Ac)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


##Boxplot
mypalette3 <- c("#2a4d69", "#4b86b4","#adcbe3","#9c0000","#fd0000","#ff9797")

K36me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(115:120)])
colnames(K36me3_rpgc_m) <- c("Sample","RPGC")
K36me3_rpgc_m$Cell <- apply(K36me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K36me3_rpgc_m$Sample <- factor(K36me3_rpgc_m$Sample, levels = c("EV0K36_ave","EV5K36_ave","EV23K36_ave","SK0K36_ave","SK5K36_ave","SK23K36_ave"))
K36me3_rpgc_m$Cell <- factor(K36me3_rpgc_m$Cell, levels = c("EV0K36","EV5K36","EV23K36","SK0K36","SK5K36","SK23K36"))
ppppp<-ggplot(K36me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K36me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 26))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

K9me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(109:114)])
colnames(K9me3_rpgc_m) <- c("Sample","RPGC")
K9me3_rpgc_m$Cell <- apply(K9me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K9me3_rpgc_m$Sample <- factor(K9me3_rpgc_m$Sample, levels = c("EV0K9_ave","EV5K9_ave","EV23K9_ave","SK0K9_ave","SK5K9_ave","SK23K9_ave"))
K9me3_rpgc_m$Cell <- factor(K9me3_rpgc_m$Cell, levels = c("EV0K9","EV5K9","EV23K9","SK0K9","SK5K9","SK23K9"))
pppp<-ggplot(K9me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K9me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pppp##boxplot


K27me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(103:108)])
colnames(K27me3_rpgc_m) <- c("Sample","RPGC")
K27me3_rpgc_m$Cell <- apply(K27me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K27me3_rpgc_m$Sample <- factor(K27me3_rpgc_m$Sample, levels = c("EV0K27_ave","EV5K27_ave","EV23K27_ave","SK0K27_ave","SK5K27_ave","SK23K27_ave"))
K27me3_rpgc_m$Cell <- factor(K27me3_rpgc_m$Cell, levels = c("EV0K27","EV5K27","EV23K27","SK0K27","SK5K27","SK23K27"))
ppp<-ggplot(K27me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K27me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 4))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppp##boxplot

K4me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(91:96)])
colnames(K4me3_rpgc_m) <- c("Sample","RPKM")
K4me3_rpgc_m$Cell <- apply(K4me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K4me3_rpgc_m$Sample <- factor(K4me3_rpgc_m$Sample, levels = c("EV0K4_ave","EV5K4_ave","EV23K4_ave","SK0K4_ave","SK5K4_ave","SK23K4_ave"))
K4me3_rpgc_m$Cell <- factor(K4me3_rpgc_m$Cell, levels = c("EV0K4","EV5K4","EV23K4","SK0K4","SK5K4","SK23K4"))
pp<-ggplot(K4me3_rpgc_m, aes(x=Sample, y=RPKM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K4me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 22))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPKM")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pp##boxplot


K27ac_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(97:102)])
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

## Wilcoxon Rank Sum test for statistical comparison

m_TE<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(11:18)])
pairwise.wilcox.test(m_TE$value,m_TE$variable)


m_K4<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(91:96)])
pairwise.wilcox.test(m_K4$value,m_K4$variable)

m_Ac<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(97:102)])
pairwise.wilcox.test(m_Ac$value,m_Ac$variable)

m_K27<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(103:108)])
pairwise.wilcox.test(m_K27$value,m_K27$variable)

m_K9<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(109:114)])
pairwise.wilcox.test(m_K9$value,m_K9$variable)

m_K36<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intergenicTE[,c(115:120)])
pairwise.wilcox.test(m_K36$value,m_K36$variable)


#######Intron TE  comparison but order on D5 TE expression (Fig. S9C)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE <-TE_cpm_anno_aveCPM_O786_FC2_both_histone[TE_cpm_anno_aveCPM_O786_FC2_both_histone$Anno2 == "intron",]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2<- cbind(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(1:3)],log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(11:18)]+0.01))
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5 <- TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$SK0D5_ave
TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5_CPM <- TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$SK0D5_ave

TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2[,c(4:11)] -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2[,c(4:11)]))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2[,c(4:11)])))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2[,c(4:11)])]

#Intronic TE expression Zscore heatmap 
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE_Zscore) 
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),])
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
colnames(df.OG2) <- colnames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),])
d<-pheatmap(df.OG2,color = colorRampPalette(c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000"))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

#CPM heatmap
Heatmap(data.matrix(2^TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),c(4:11)]), column_title = "Timepoint",name= "intron TEs\n(CPM)",col = colorRamp2(c(0,0.25,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5), c(rev(mypalette7))), 
        cluster_rows = F, cluster_columns = FALSE,show_row_names = F)

#intron TE expression CPM boxplot (Fig 4C)
TE_aveCPM_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone[TE_cpm_anno_aveCPM_O786_FC2_both_histone$Anno2 == "intron",c(11:18)])
colnames(TE_aveCPM_m) <- c("Sample","CPM")
TE_aveCPM_m$Cell <- apply(TE_aveCPM_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
TE_aveCPM_m$Sample <- factor(TE_aveCPM_m$Sample, levels = c("OEVD0_ave","OEVD5_ave","OEVD15_ave","OEVD23_ave","SK0D0_ave","SK0D5_ave","SK0D15_ave","SK0D23_ave"))
TE_aveCPM_m$Cell <- factor(TE_aveCPM_m$Cell, levels = c("OEVD0","OEVD5","OEVD15","OEVD23","SK0D0","SK0D5","SK0D15","SK0D23"))
ppppp<-ggplot(TE_aveCPM_m, aes(x=Sample, y=CPM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("Intergenic TE Expression Distribution (CPM)")+ coord_cartesian(ylim=c(0, 3))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "CPM")+scale_fill_manual(values = c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+scale_color_manual(values=c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

### histone Zscore heatmaps by intron TE Zscore ###
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),c(91:96)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),c(97:102)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),c(103:108)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),c(109:114)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9)]
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36 <- log2(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[order(-TE_cpm_anno_aveCPM_O786_FC2_both_histone_log2_intronTE2$diffEV5_KO5),c(115:120)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36)))[row(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36)]

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K9)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =T, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K36)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K4)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_K27)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE_log2_Ac)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


##Boxplot
mypalette3 <- c("#2a4d69", "#4b86b4","#adcbe3","#9c0000","#fd0000","#ff9797")

K36me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(115:120)])
colnames(K36me3_rpgc_m) <- c("Sample","RPGC")
K36me3_rpgc_m$Cell <- apply(K36me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K36me3_rpgc_m$Sample <- factor(K36me3_rpgc_m$Sample, levels = c("EV0K36_ave","EV5K36_ave","EV23K36_ave","SK0K36_ave","SK5K36_ave","SK23K36_ave"))
K36me3_rpgc_m$Cell <- factor(K36me3_rpgc_m$Cell, levels = c("EV0K36","EV5K36","EV23K36","SK0K36","SK5K36","SK23K36"))
ppppp<-ggplot(K36me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K36me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 26))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

K9me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(109:114)])
colnames(K9me3_rpgc_m) <- c("Sample","RPGC")
K9me3_rpgc_m$Cell <- apply(K9me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K9me3_rpgc_m$Sample <- factor(K9me3_rpgc_m$Sample, levels = c("EV0K9_ave","EV5K9_ave","EV23K9_ave","SK0K9_ave","SK5K9_ave","SK23K9_ave"))
K9me3_rpgc_m$Cell <- factor(K9me3_rpgc_m$Cell, levels = c("EV0K9","EV5K9","EV23K9","SK0K9","SK5K9","SK23K9"))
pppp<-ggplot(K9me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K9me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pppp##boxplot


K27me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(103:108)])
colnames(K27me3_rpgc_m) <- c("Sample","RPGC")
K27me3_rpgc_m$Cell <- apply(K27me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K27me3_rpgc_m$Sample <- factor(K27me3_rpgc_m$Sample, levels = c("EV0K27_ave","EV5K27_ave","EV23K27_ave","SK0K27_ave","SK5K27_ave","SK23K27_ave"))
K27me3_rpgc_m$Cell <- factor(K27me3_rpgc_m$Cell, levels = c("EV0K27","EV5K27","EV23K27","SK0K27","SK5K27","SK23K27"))
ppp<-ggplot(K27me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K27me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 4))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppp##boxplot

K4me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(91:96)])
colnames(K4me3_rpgc_m) <- c("Sample","RPKM")
K4me3_rpgc_m$Cell <- apply(K4me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K4me3_rpgc_m$Sample <- factor(K4me3_rpgc_m$Sample, levels = c("EV0K4_ave","EV5K4_ave","EV23K4_ave","SK0K4_ave","SK5K4_ave","SK23K4_ave"))
K4me3_rpgc_m$Cell <- factor(K4me3_rpgc_m$Cell, levels = c("EV0K4","EV5K4","EV23K4","SK0K4","SK5K4","SK23K4"))
pp<-ggplot(K4me3_rpgc_m, aes(x=Sample, y=RPKM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K4me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 22))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPKM")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pp##boxplot


K27ac_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(97:102)])
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


### Wilcoxon Rank Sum Test for statistical comparision
m_TE<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(11:18)])
pairwise.wilcox.test(m_TE$value,m_TE$variable)

m_K4<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(91:96)])
pairwise.wilcox.test(m_K4$value,m_K4$variable)

m_Ac<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(97:102)])
pairwise.wilcox.test(m_Ac$value,m_Ac$variable)

m_K27<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(103:108)])
pairwise.wilcox.test(m_K27$value,m_K27$variable)

m_K9<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(109:114)])
pairwise.wilcox.test(m_K9$value,m_K9$variable)

m_K36<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_both_histone_intronTE[,c(115:120)])
pairwise.wilcox.test(m_K36$value,m_K36$variable)



############# Mis-spliced TE expression and histone analysis (N=196) (Fig. S11A) #######
#all time points comparison but order on D5 TE expression (Zscore)
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone <- TE_cpm_anno_aveCPM_O786_FC2_both_histone[TE_cpm_anno_aveCPM_O786_FC2_both_histone$"Geneid" %in% unique(TE_denovo_allKO_exon[TE_denovo_allKO_exon$Anno2 == "intron",c("Geneid")]),]

TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2<- cbind(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(1:3)],log2(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(11:18)]+0.01))
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5 <- TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$OEVD5_ave -TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$SK0D5_ave

TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2[,c(4:11)] -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2[,c(4:11)]))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2[,c(4:11)])))[row(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2[,c(4:11)])]

#Missplicd TE expression Zscore heatmap
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE_Zscore) #N = 196 TEs
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),])
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
colnames(df.OG2) <- colnames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE_Zscore[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),])
d<-pheatmap(df.OG2,color = colorRampPalette(c("#005073","#71c7ec", "#fffbea","#ffbaba","#ff0000"))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d

#CPM heatmap (of Z-score)
Heatmap(data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),c(11:18)]), column_title = "Timepoint",name= "all TEs\n(CPM)",col = colorRamp2(c(0,0.25,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5), c(rev(mypalette7))), 
        cluster_rows = F, cluster_columns = FALSE,show_row_names = F)

#TE expression CPM boxplot 
TE_aveCPM_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(11:18)])
colnames(TE_aveCPM_m) <- c("Sample","CPM")
TE_aveCPM_m$Cell <- apply(TE_aveCPM_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
TE_aveCPM_m$Sample <- factor(TE_aveCPM_m$Sample, levels = c("OEVD0_ave","OEVD5_ave","OEVD15_ave","OEVD23_ave","SK0D0_ave","SK0D5_ave","SK0D15_ave","SK0D23_ave"))
TE_aveCPM_m$Cell <- factor(TE_aveCPM_m$Cell, levels = c("OEVD0","OEVD5","OEVD15","OEVD23","SK0D0","SK0D5","SK0D15","SK0D23"))
ppppp<-ggplot(TE_aveCPM_m, aes(x=Sample, y=CPM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("TE Expression Distribution (CPM)")+ coord_cartesian(ylim=c(0, 8))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "CPM")+scale_fill_manual(values = c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+scale_color_manual(values=c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot


### histone Zscore heatmaps by TE Zscore###
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4 <- log2(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),c(91:96)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4)))[row(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4)]
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac <- log2(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),c(97:102)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac)))[row(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac)]
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27 <- log2(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),c(103:108)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27)))[row(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27)]
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9 <- log2(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),c(109:114)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9)))[row(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9)]
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36 <- log2(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[order(-TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_TE2$diffEV5_KO5),c(115:120)]+0.01)
TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36_Zscore<- (TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36 -rowMeans(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36))/(rowSds(as.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36)))[row(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36)]


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =T, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K9)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =T, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1

name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K36)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =T, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K4)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_K27)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac_Zscore)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac_Zscore)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d
name <- rownames(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac)
df.OG2 <- data.matrix(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone_log2_Ac)
row.names(df.OG2) <- name
df.OG2[is.na(df.OG2)] <-0
d1<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = F,cluster_rows =F, cluster_cols=F,border_color=NA,clustering_method ="ward.D2")
d1


##Boxplot
mypalette3 <- c("#2a4d69", "#4b86b4","#adcbe3","#9c0000","#fd0000","#ff9797")

K36me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(115:120)])
colnames(K36me3_rpgc_m) <- c("Sample","RPGC")
K36me3_rpgc_m$Cell <- apply(K36me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K36me3_rpgc_m$Sample <- factor(K36me3_rpgc_m$Sample, levels = c("EV0K36_ave","EV5K36_ave","EV23K36_ave","SK0K36_ave","SK5K36_ave","SK23K36_ave"))
K36me3_rpgc_m$Cell <- factor(K36me3_rpgc_m$Cell, levels = c("EV0K36","EV5K36","EV23K36","SK0K36","SK5K36","SK23K36"))
ppppp<-ggplot(K36me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K36me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 26))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppppp##boxplot

K9me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(109:114)])
colnames(K9me3_rpgc_m) <- c("Sample","RPGC")
K9me3_rpgc_m$Cell <- apply(K9me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K9me3_rpgc_m$Sample <- factor(K9me3_rpgc_m$Sample, levels = c("EV0K9_ave","EV5K9_ave","EV23K9_ave","SK0K9_ave","SK5K9_ave","SK23K9_ave"))
K9me3_rpgc_m$Cell <- factor(K9me3_rpgc_m$Cell, levels = c("EV0K9","EV5K9","EV23K9","SK0K9","SK5K9","SK23K9"))
pppp<-ggplot(K9me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K9me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pppp##boxplot


K27me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(103:108)])
colnames(K27me3_rpgc_m) <- c("Sample","RPGC")
K27me3_rpgc_m$Cell <- apply(K27me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K27me3_rpgc_m$Sample <- factor(K27me3_rpgc_m$Sample, levels = c("EV0K27_ave","EV5K27_ave","EV23K27_ave","SK0K27_ave","SK5K27_ave","SK23K27_ave"))
K27me3_rpgc_m$Cell <- factor(K27me3_rpgc_m$Cell, levels = c("EV0K27","EV5K27","EV23K27","SK0K27","SK5K27","SK23K27"))
ppp<-ggplot(K27me3_rpgc_m, aes(x=Sample, y=RPGC)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K27me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 4))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPGC")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ppp##boxplot

K4me3_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(91:96)])
colnames(K4me3_rpgc_m) <- c("Sample","RPKM")
K4me3_rpgc_m$Cell <- apply(K4me3_rpgc_m, 1, function(x) {unlist(strsplit(x[1],"_"))[1]})
K4me3_rpgc_m$Sample <- factor(K4me3_rpgc_m$Sample, levels = c("EV0K4_ave","EV5K4_ave","EV23K4_ave","SK0K4_ave","SK5K4_ave","SK23K4_ave"))
K4me3_rpgc_m$Cell <- factor(K4me3_rpgc_m$Cell, levels = c("EV0K4","EV5K4","EV23K4","SK0K4","SK5K4","SK23K4"))
pp<-ggplot(K4me3_rpgc_m, aes(x=Sample, y=RPKM)) +
  geom_boxplot(aes(fill=Cell),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA) +ggtitle("H3K4me3 Peak Size Distribution (RPGC)")+ coord_cartesian(ylim=c(0, 22))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "RPKM")+scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pp##boxplot


K27ac_rpgc_m<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(97:102)])
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

######Statistics 
m_TE<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(11:18)])
pairwise.wilcox.test(m_TE$value,m_TE$variable)

m_K4<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(91:96)])
pairwise.wilcox.test(m_K4$value,m_K4$variable)

m_Ac<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(97:102)])
pairwise.wilcox.test(m_Ac$value,m_Ac$variable)

m_K27<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(103:108)])
pairwise.wilcox.test(m_K27$value,m_K27$variable)

m_K9<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(109:114)])
pairwise.wilcox.test(m_K9$value,m_K9$variable)

m_K36<-reshape2::melt(TE_cpm_anno_aveCPM_O786_FC2_misspliced_histone[,c(115:120)])
pairwise.wilcox.test(m_K36$value,m_K36$variable)





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

#Use Deeptools to generate Fig. 4A.
computeMatrix scale-regions -bs 200 -p 12 -R /genomes/hg38/gencode.v37.AutosomalTranscriptsOnly.forDeepTools.gtf -S EV_D0_CpG_Methylation_Beta_values.sorted.bw EV_D5_CpG_Methylation_Beta_values.sorted.bw EV_D15_CpG_Methylation_Beta_values.sorted.bw EV_D39_CpG_Methylation_Beta_values.sorted.bw KO_D0_CpG_Methylation_Beta_values.sorted.bw KO_D5_CpG_Methylation_Beta_values.sorted.bw KO_D15_CpG_Methylation_Beta_values.sorted.bw KO_D39_CpG_Methylation_Beta_values.sorted.bw -m 5000 -b 3000 -a 3000 -out DNAmethylation_GeneBody_O786_BetaValues.tab.gz 
plotHeatmap  --heatmapWidth 6 -m DNAmethylation_GeneBody_O786_BetaValues.tab.gz --perGroup -out DNAmethylation_GeneBody_O786_BetaValues.pdf --missingDataColor "white" -min 0 -max 1 --colorList "blue,yellow,red" --heatmapHeight 10

## Distribution of probe methylation
meth_anno_CpG <- meth_anno[meth_anno$probeType == "cg" & meth_anno$MASK_general==F,c("EV_D0","EV_D5","EV_D15","EV_D39","KO_D0","KO_D5","KO_D15","KO_D39")] #filter masked probes
m_meth_anno_CpG <-reshape2::melt(meth_anno_CpG)
m_meth_anno_CpG$variable <- factor(m_meth_anno_CpG$variable, levels = c("EV_D0","EV_D5","EV_D15","EV_D39","KO_D0","KO_D5","KO_D15","KO_D39"))

mypalette6 <- c("#2a4d69", "#4b86b4","#73a7d0","#adcbe3","#9c0000","#fd0000","#fb6161","#ff9797")
#Fig. S5A
ggplot(m_meth_anno_CpG, aes(x=variable, y=value)) +
  geom_violin(aes(fill=variable),adjust =1)+geom_boxplot(fill="#e7e7e7",position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA,width = 0.1) +ggtitle("DNA methylation distribution (Beta)")+ coord_cartesian(ylim=c(0, 1))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "Beta")+scale_fill_manual(values = mypalette6)+scale_color_manual(values=mypalette6)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")


#####TE centric methylation analysis#####
bedtools intersect -wo -a TE_FC2_O786_location.bed -b Allsample_CpG_Probe_EPICmethylation_BetaValues.bed > TE_FC2_O786_allsample_CpG_Probe_EPICmethylation_BetaValues_V37.bed

bothTE_2FC_BS <-read.delim("TE_FC2_O786_allsample_CpG_Probe_EPICmethylation_BetaValues_V37.bed", header = F, sep = "\t",stringsAsFactors = F)

#Fig. S9A
pheatmap(data.matrix(bothTE_2FC_BS[bothTE_2FC_BS$V9 == "Intergenic",c(10:13,15:18)]),color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),breaks = c(seq(0,1,by=0.01)),show_rownames = F,cluster_rows =T, cluster_cols=F,clustering_method="ward.D2")
pheatmap(data.matrix(bothTE_2FC_BS[bothTE_2FC_BS$V9 == "intron",c(10:13,15:18)]),color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),breaks = c(seq(0,1,by=0.01)),show_rownames = F,cluster_rows =T, cluster_cols=F,clustering_method="ward.D2")








