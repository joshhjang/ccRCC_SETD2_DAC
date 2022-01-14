# hemantgujar@yahoo.om
# 2019_10_07
# KIRC hERV data analysis 

wdir="C:/Users/Hemant Gujar/Desktop/Kidney_cancer/KIRC/RNASeq_EData/analysis/hERV"
setwd(wdir)
hERV <- read.table(file = "C:/Users/Hemant Gujar/Desktop/pHGG/scripts/RcFun/GDC/DataSet/hERV/Smith2018/JCI121476.sdt12.txt", sep = "\t", header = T)

KIRC_metadata <- read.table("C:/Users/Hemant Gujar/Desktop/Kidney_cancer/KIRC/RNASeq_EData/metadata/KIRC_metadata_S.csv", sep = ",", header = T)
samples <- as.character(KIRC_metadata$TCGA_ID)
samples <- gsub(".{1}$", "", samples)                     
samples <- samples[!duplicated(samples)]

hERV_KIRC <- hERV[(hERV$Sample_ID %in% samples),]
rownames(hERV_KIRC) <- hERV_KIRC$Sample_ID
hERV_KIRC <- hERV_KIRC[,-1]
hERV_KIRC <- t(hERV_KIRC)
hERV_KIRC <- as.data.frame(hERV_KIRC)

lowSETD2 <- as.character(KIRC_metadata[(KIRC_metadata$SETD2_lvl == "low"),"TCGA_ID"])
lowSETD2 <- gsub(".{1}$", "", lowSETD2)            
lowSETD2 <- lowSETD2[!duplicated(lowSETD2)]
lowSETD2 <- hERV_KIRC[,lowSETD2]
lowSETD2 <- as.data.frame(lowSETD2)

highSETD2 <- as.character(KIRC_metadata[(KIRC_metadata$SETD2_lvl == "high"),"TCGA_ID"])
highSETD2 <- gsub(".{1}$", "", highSETD2)   
highSETD2 <- highSETD2[!duplicated(highSETD2)]
highSETD2 <- hERV_KIRC[,(names(hERV_KIRC) %in% highSETD2)]
hERV_KIRC <- cbind(lowSETD2, highSETD2)

############################################################################################
# Make PCA plot 

tmp <- t(hERV_KIRC)
tmp <- tmp[,hERV_KIRC_names]

metadata <- KIRC_metadata[grep("low|high", KIRC_metadata$SETD2_lvl),c(1,3,6,8)]
metadata$TCGA_ID <- gsub(".{1}$", "", metadata$TCGA_ID)

tmp <- merge(metadata, tmp, by.x = 1, by.y = "row.names", all=F)
tmp <- tmp[!duplicated(tmp$TCGA_ID),]

library(ggfortify)
autoplot(prcomp(tmp[,c(5:ncol(tmp))]), data = tmp, label =F, colour = "SETD2_lvl", shape=19,size=2) + 
  theme_bw() +  scale_color_manual(values = c("red", "green4"))

# Make heat map 
rownames(tmp) <- tmp$TCGA_ID

hERV_KIRC_sel <- t(tmp[,-c(1:4)])
hERV_KIRC_sel <- as.data.frame(hERV_KIRC_sel)
metadata <- tmp[,c(1:4)]
metadata <- metadata[order(metadata$SETD2_lvl),]
col_order <- as.character(metadata$TCGA_ID)
metadata <- metadata[,c(4),drop=F]

hERV_KIRC_sel <- hERV_KIRC_sel[,col_order]
hERV_KIRC_sel <- log2(hERV_KIRC_sel + 0.005)

hERV_KIRC_mean <- colMeans(hERV_KIRC_sel)
hERV_KIRC_sd <- apply(hERV_KIRC_sel,2,sd)

library(ComplexHeatmap)
library(circlize)
hc1 <- HeatmapAnnotation(df = metadata, which = "column", col = list(SETD2_lvl = c("high"="red", "low"="green4")),
                         AverageExp = anno_lines(hERV_KIRC_mean, gp = gpar(col = 1:4), add_points = F), height = unit(2, "cm"))

h_hERV_KIRC <- Heatmap(hERV_KIRC_sel, km=4, col = colorRamp2(c(-4.50, -2.00, 0.50),  c("green4", "white", "red")),
                       cluster_columns = F, clustering_distance_columns = "euclidean",  
                       clustering_method_columns = "complete", column_dend_side = c("top"), column_dend_height = unit(10, "mm"),
                       show_column_dend = TRUE, show_column_names = FALSE, cluster_rows = T, clustering_distance_rows = "euclidean",
                       clustering_method_rows = "complete", row_dend_side = c("left"), row_dend_width = unit(10, "mm"),
                       show_row_dend = TRUE, row_dend_reorder = TRUE, row_dend_gp = gpar(), 
                       column_title = "KIRC", column_title_side = "bottom", row_title = "hERV", 
                       name = "log2count", show_row_names = TRUE, top_annotation = hc1)
draw(h_hERV_KIRC)

##### make boxplot 
hERV_KIRC_sel <- t(hERV_KIRC_sel)
hERV_KIRC_sel <- as.data.frame(hERV_KIRC_sel)
hERV_KIRC_sel <- merge(metadata, hERV_KIRC_sel, by="row.names", all = F)
colnames(hERV_KIRC_sel)[1] <- "TCGA_ID"
rownames(hERV_KIRC_sel) <- hERV_KIRC_sel$Row.names

######## make box plot with average values 

tmp <- hERV_KIRC_sel[,c(1:2)]
tmp$mean <- rowMeans(hERV_KIRC_sel[,-c(1:2)])

library(ggplot2)
ggplot(tmp, aes(x=SETD2_lvl, y=mean, fill=SETD2_lvl)) + scale_fill_manual(values=c("red", "green4")) +
  geom_boxplot(width=0.5) + 
  geom_jitter(shape=16, position=position_jitter(0.1), size=0.01) + 
  theme_classic() + labs(title="KIRC hERV expression",x="SETD2 lvl", y = "Mean_log2NCount") + theme(plot.title = element_text(hjust = 0.5))

###########################################################################################

# Genes involved in RIG-I/MDA5 mediated induction of IFN-alpha/beta pathways
# from http://software.broadinstitute.org/gsea/msigdb/cards/REACTOME_RIG_I_MDA5_MEDIATED_INDUCTION_OF_IFN_ALPHA_BETA_PATHWAYS.html 

wdir="C:/Users/Hemant Gujar/Desktop/Kidney_cancer/KIRC/RNASeq_EData/analysis/hERV"

KIRC_cpm <- read.table("C:/Users/Hemant Gujar/Desktop/Kidney_cancer/KIRC/RNASeq_EData/analysis/tables/KIRC_cpm.tsv", sep="\t", header=T)

KIRC_metadata <- read.table("C:/Users/Hemant Gujar/Desktop/Kidney_cancer/KIRC/RNASeq_EData/metadata/KIRC_metadata_S.csv", sep = ",", header = T)
samples <- as.character(KIRC_metadata$TCGA_ID)
samples <- gsub(".{1}$", "", samples)                            # subsitute last character in the string
samples <- samples[!duplicated(samples)]

lowSETD2 <- as.character(KIRC_metadata[(KIRC_metadata$SETD2_lvl == "low"),"TCGA_ID"])
lowSETD2 <- lowSETD2[!duplicated(lowSETD2)]
length(lowSETD2)
lowSETD2 <- KIRC_cpm[,lowSETD2]
dim(lowSETD2)

highSETD2 <- as.character(KIRC_metadata[(KIRC_metadata$SETD2_lvl == "high"),"TCGA_ID"])
highSETD2 <- highSETD2[!duplicated(highSETD2)]
highSETD2 <- KIRC_cpm[,highSETD2]
anno <- KIRC_cpm[,c(1:5)]
head(anno) 

RIG <- c("AGER", "APP", "ATG12", "ATG5", "CASP10", "CASP8", "CHUK", "CREBBP", "CYLD", "DDX58", "DHX58", "EP300", "FADD", "HERC5", "HMGB1", "IFIH1", "IFNA1", "IFNA10", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNB1", "IKBKB", "IKBKE", "IKBKG", "IRF1", "IRF2", "IRF3", "IRF7", "ISG15", "MAP3K1", "MAVS", "NFKB2", "NFKBIA", "NFKBIB", "NLRC5", "NLRX1", "OTUD5", "PCBP2", "PIN1", "RELA", "RIPK1", "RNF125", "RNF135", "RPS27A", "RPS27AP11", "S100A12", "S100B", "SAA1", "SIKE1", "TANK", "TAX1BP1", "TBK1", "TKFC", "TNFAIP3", "TRAF2", "TRAF3", "TRAF6", "TRIM25", "UBA52", "UBA7", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2K", "UBE2L6")
RIG_anno <- anno[(anno$hgnc_symbol %in% RIG),]
RIG_ENSG <- as.character(RIG_anno$ensembl_gene_id)

KIRC_stat <- read.table("../tables/KIRC_SETD2_stats.csv", sep = ",", header = T)
KIRC_stat <- KIRC_stat[(KIRC_stat$ENSG %in% RIG_ENSG),]
KIRC_stat <- KIRC_stat[(KIRC_stat$p_adj < 0.05),]
KIRC_stat <- subset(KIRC_stat, KIRC_stat$fold_change > 1.5 | KIRC_stat$fold_change < 0.67)
RIG_ENSG <- as.character(KIRC_stat$ENSG)

KIRC_metadata <- read.table("C:/Users/Hemant Gujar/Desktop/Kidney_cancer/KIRC/RNASeq_EData/metadata/KIRC_metadata_S.csv", sep = ",", header = T)
KIRC_metadata$TCGA_ID <- gsub(".{1}$", "", KIRC_metadata$TCGA_ID)
KIRC_metadata <- KIRC_metadata[,c(1,3,8)]
KIRC_metadata <- KIRC_metadata[!(KIRC_metadata$SETD2_lvl == "mid"),]
KIRC_metadata <- KIRC_metadata[order(KIRC_metadata$TCGA_ID),]
KIRC_metadata <- KIRC_metadata[!duplicated(KIRC_metadata$TCGA_ID),]
KIRC_metadata <- KIRC_metadata[order(KIRC_metadata$SETD2_lvl, decreasing = F),]
rownames(KIRC_metadata) <- KIRC_metadata$TCGA_ID
KIRC_metadata <- KIRC_metadata[(KIRC_metadata$sample_type == "Primary"),]
KIRC_metadata <- KIRC_metadata[,-c(1,2),drop=F]
metadata <- KIRC_metadata

col_order <- rownames(KIRC_metadata)

# high SETD2 Vs low SETD2
dim(lowSETD2)
dim(highSETD2)
KIRC_cpm_sel <- cbind(lowSETD2, highSETD2)
KIRC_cpm_sel <- KIRC_cpm_sel[(rownames(KIRC_cpm_sel) %in% RIG_ENSG),]
KIRC_cpm_sel <- merge(RIG_anno, KIRC_cpm_sel, by.x = "ensembl_gene_id", by.y = "row.names", all = F)
rownames(KIRC_cpm_sel) <- KIRC_cpm_sel$hgnc_symbol
KIRC_cpm_sel <- KIRC_cpm_sel[,-c(1:5)]

names(KIRC_cpm_sel) <- gsub(".{1}$", "", names(KIRC_cpm_sel))
KIRC_cpm_sel <- KIRC_cpm_sel[,col_order]

library(ComplexHeatmap)
library(circlize)
hc1 <- HeatmapAnnotation(df = metadata, which = "column", col = list(SETD2_lvl = c("high"="red", "low"="green4")))
h_RIG_KIRC <- Heatmap(log2(KIRC_cpm_sel + 0.25), km=1, col = colorRamp2(c(3.50, 5.50, 8.00),  c("green4", "white", "red")),
                      cluster_columns = F, clustering_distance_columns = "euclidean",  
                      clustering_method_columns = "complete", column_dend_side = c("top"), column_dend_height = unit(10, "mm"),
                      show_column_dend = TRUE, show_column_names = FALSE, cluster_rows = T, clustering_distance_rows = "euclidean",
                      clustering_method_rows = "complete", row_dend_side = c("left"), row_dend_width = unit(10, "mm"),
                      show_row_dend = TRUE, row_dend_reorder = TRUE, row_dend_gp = gpar(), 
                      column_title = "KIRC", column_title_side = "bottom", row_title = "RIG genes", 
                      name = "log2cpm", show_row_names = TRUE, top_annotation = hc1)
draw(h_RIG_KIRC)

#############################################################################################

