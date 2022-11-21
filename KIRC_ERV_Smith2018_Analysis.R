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
# Genes involved in KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY were utilized to generate heatmap of RIG-I pathway

tum <- read.delim("KIRC_RNA_Final-Highlow.txt",head=TRUE,sep="\t")
head(tum)

row.names(tum) <- tum[,1]

data <-tum[,c(2:233)]
head (data)

data <- as.matrix (data)
#a <- data[,c(2,3)]
b <- scale(t(data))


b <-ifelse(b < -3, -3, b)
b <-ifelse(b > 3, 3, b)

e <-t(b) 

d1 <- dist(e, method = "euclidean")
fit <- hclust(d1, method="ward.D")
plot(fit) # display dendogram
fit_den <-as.dendrogram(fit)


png("TCGA KIRC immune SETD2 Fin.pdf")
heatmap.2(e, scale="none",col="bluered", density.info="none",trace="none", dendrogram="row", 
          Rowv=fit_den,  symbreaks=T,keysize = 1, key=T)

dev.off()
#############################################################################################

