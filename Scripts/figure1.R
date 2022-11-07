library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Signac)
library(ensembldb)

#FIGURE 1A - GD15.5
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-ATACwithRNAlables.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/scRNA-seqFiles/gd15.5_res0.8.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")

data <- RenameIdents(data, `0` = "Other cells", `1` = "Other cells",
                     `2` = "Other cells", `3` = "Other cells", `4` = "Other cells",
                     `5` = "Other cells", `6` = "Natural killer cells",
                     `7` = "Macrophage cells", `8` = "Macrophage cells",
                     `9` = "Other cells", `10` = "Endothelial cells",
                     `11` = "Macrophage cells", `12` = "Other cells",
                     `13` = "Endothelial cells", `14` = "Smooth muscle cells",
                     `15` = "Endothelial cells", `16` = "Other cells",
                     `17` = "Smooth muscle cells", `18` = "Other cells",
                     `19` = "Natural killer cells", `20` = "Other cells",
                     `21` = "Other cells", `22` = "Natural killer cells",
                     `23` = "Invasive trophoblast cells", `24` = "Other cells",
                     `25` = "Other cells", `26` = "Smooth muscle cells",
                     `27` = "Other cells", `28` = "Other cells")
Idents(data) <- gsub(pattern = " ", replacement = "_", x = Idents(data))
Idents(data) <- factor(Idents(data), levels = c("Other_cells", "Macrophage_cells", "Smooth_muscle_cells",
                                                "Endothelial_cells", "Natural_killer_cells", "Invasive_trophoblast_cells"))
rats.RNA <- data

Idents(rats) <- rats@meta.data$predicted.id
Idents(rats) <- gsub(pattern = " ", replacement = "_", x = Idents(rats))
Idents(rats) <- factor(Idents(rats), levels = c("Other_cells", "Macrophage_cells", "Smooth_muscle_cells",
                                                "Endothelial_cells", "Natural_killer_cells", "Invasive_trophoblast_cells"))

plot1 <- DimPlot(object = rats.RNA, cols = c("grey70", "orange3", "#718748", "#183D61", "#79704C", "red3")) + NoLegend() + ggtitle('GD15.5 scRNA-seq')
plot2 <- DimPlot(object=rats, cols = c("grey70", "orange3", "#718748", "#183D61", "#79704C", "red3")) + ggtitle('GD15.5 scATAC-seq')

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure1a-gd15.5-dimplot-ATAC-RNAlabels.pdf", width=18, height=8)
plot1 + plot2
dev.off()

#FIGURE 1A - GD19.5
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-ATACwithRNAlables.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/scRNA-seqFiles/gd19.5-4-5-6-7_res0.8.rda")
data <- data2

data <- RenameIdents(data, `0` = "Other cells", `1` = "Macrophage cells",
                     `2` = "Other cells", `3` = "Other cells", `4` = "Other cells",
                     `5` = "Other cells", `6` = "Invasive trophoblast cells",
                     `7` = "Macrophage cells", `8` = "Smooth muscle cells",
                     `9` = "Natural killer cells", `10` = "Endothelial cells",
                     `11` = "Macrophage cells", `12` = "Other cells",
                     `13` = "Other cells", `14` = "Macrophage cells",
                     `15` = "Other cells", `16` = "Other cells",
                     `17` = "Macrophage cells", `18` = "Macrophage cells",
                     `19` = "Other cells", `20` = "Endothelial cells",
                     `21` = "Other cells", `22` = "Other cells",
                     `23` = "Other cells", `24` = "Other cells",
                     `25` = "Other cells", `26` = "Smooth muscle cells")

Idents(data) <- gsub(pattern = " ", replacement = "_", x = Idents(data))
Idents(data) <- factor(Idents(data), levels = c("Other_cells", "Macrophage_cells", "Smooth_muscle_cells",
                                                "Endothelial_cells", "Natural_killer_cells", "Invasive_trophoblast_cells"))
rats.RNA <- data

Idents(rats) <- rats@meta.data$predicted.id
Idents(rats) <- gsub(pattern = " ", replacement = "_", x = Idents(rats))
Idents(rats) <- factor(Idents(rats), levels = c("Other_cells", "Macrophage_cells", "Smooth_muscle_cells",
                                                "Endothelial_cells", "Natural_killer_cells", "Invasive_trophoblast_cells"))

plot1 <- DimPlot(object = rats.RNA, cols = c("grey70", "orange3", "#718748", "#183D61", "#79704C", "red3")) + NoLegend() + ggtitle('GD19.5 scRNA-seq')
plot2 <- DimPlot(object=rats, cols = c("grey70", "orange3", "#718748", "#183D61", "#79704C", "red3")) + ggtitle('GD19.5 scATAC-seq')

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure1a-gd19.5-dimplot-ATAC-RNAlabels.pdf", width=18, height=8)
plot1 + plot2
dev.off()

#figure1b
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-ATACwithRNAlables.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")
Idents(rats) <- rats@meta.data$predicted.id
DefaultAssay(rats) <- "RNA"
features <- c("Prl5a1", "Prl7b1", "Krt7", "Krt8", "Krt18",
              "Nkg7","Prf1","Gzmbl2","Gzmm","Cdh5",
              "Plvap", "Adgrl4", "Egfl7", 
              "Acta2", "Tagln", "Myh11", "Myl9",
              "C1qa", "Lyz2", "Aif1", "Cybb")
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/Fig1b-gd15.5_dotPlotMarkers.pdf", width = 12, height = 8)
#DotPlot(rats, features = features, idents = c("Invasive trophoblast cells", "Natural killer cells", "Endothelial cells", "Smooth muscle cells", "Macrophage cells"), 
#        cols = c("red3", "#79704C", "#183D61", "#4E691A", "orange3"), split.by = "predicted.id") + 
#  ggtitle("GD15.5 Cell Type Markers") +
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 25),
#        panel.grid.major= element_blank(),panel.grid.minor = element_blank())
DotPlot(rats, features = features, cols = c("#c6dbef", "#08519c"),
        idents = c("Invasive trophoblast cells", "Natural killer cells", "Endothelial cells", "Smooth muscle cells", "Macrophage cells")) + 
  ggtitle("GD15.5 Cell Type Markers") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 25),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())

dev.off()

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-ATACwithRNAlables.rda")
Idents(rats) <- rats@meta.data$predicted.id
DefaultAssay(rats) <- "RNA"

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/Fig1b-gd19.5_dotPlotMarkers.pdf", width = 12, height = 8)
#DotPlot(rats, features = features, idents = c("Invasive trophoblast cells", "Natural killer cells", "Endothelial cells", "Smooth muscle cells", "Macrophage cells"), 
#        cols = c("red3", "#79704C", "#183D61", "#4E691A", "orange3"), split.by = "predicted.id") + 
#  ggtitle("GD19.5 Cell Type Markers") +
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 25),
#        panel.grid.major= element_blank(),panel.grid.minor = element_blank())
DotPlot(rats, features = features, cols = c("#c6dbef", "#08519c"),
        idents = c("Invasive trophoblast cells", "Natural killer cells", "Endothelial cells", "Smooth muscle cells", "Macrophage cells")) + 
  ggtitle("GD19.5 Cell Type Markers") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 25),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())
dev.off()


#figure1c
library(ggplot2)
df <- data.frame(cell=rep(c("Invasive trophoblast", "Natural killer", "Endothelium", "Macrophage", "Smooth muscle"),2),
                 no=c(3842, 1503, 670, 652, 267, 2594, 1172, 565, 534, 237),
                 type=c(rep("Peaks",5), rep("Genes",5)))
p1<-ggplot(data=df, aes(x=cell, y=no, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=no), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Number of cell specific peaks and closest genes - gd 15.5")

df <- data.frame(cell=rep(c("Invasive trophoblast", "Natural killer", "Endothelium", "Macrophage", "Smooth muscle"),2),
                 no=c(2033, 279, 686, 99, 345, 1535, 253, 574, 92, 311),
                 type=c(rep("Peaks",5), rep("Genes",5)))

p2<-ggplot(data=df, aes(x=cell, y=no, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=no), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Number of cell specific peaks and closest genes - gd 19.5")

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/Fig1c-numberPeaksGenes.pdf", width = 10, height = 5)
p1+p2
dev.off()


#figure 1d
library(ChIPseeker)
library(GenomicFeatures)
library(RMariaDB)
library(dplyr)

txdb <- makeTxDbFromEnsembl(organism= "Rattus norvegicus", release = 98)

#gd15.5
file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-nk.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pNK <- plotDistToTSS(peakAnno)
pNK <- as.data.frame(pNK$data)

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-sm.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pSM <- plotDistToTSS(peakAnno)
pSM <- as.data.frame(pSM$data)
df <- inner_join(pNK, pSM, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-endo.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pEndo <- plotDistToTSS(peakAnno)
pEndo <- as.data.frame(pEndo$data)
df <- inner_join(df, pEndo, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-macro.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pMacro <- plotDistToTSS(peakAnno)
pMacro <- as.data.frame(pMacro$data)
df <- inner_join(df, pMacro, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-tbc.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pTBC <- plotDistToTSS(peakAnno)
pTBC <- as.data.frame(pTBC$data)
df <- inner_join(df, pTBC, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

df$mean <- rowMeans(df[,c("freq.x", "freq.y", "freq.x.x", "freq.y.y", "freq")])
df <- df[,c(1,2,4,9)]
df$Feature <- factor(df$Feature, levels=levels(df$Feature))
df <- df[order(-df$.group),]
color <- c("#d73027", "#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4")
names(color) <- unique(df$Feature)
p <- ggplot(df, aes(x=1, fill=Feature))
p <- p + geom_bar(data=subset(df, sign==1), aes(y=mean, fill = Feature), stat="identity") +
  geom_bar(data=subset(df, sign==-1), aes(y=-mean, fill = Feature), stat="identity") +
  scale_fill_manual(values=color)
p <- p + geom_hline(yintercept = 0, colour = "black") +
  coord_flip() + theme_bw() + ggtitle("Average distance distribution at gd 15.5") +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) + 
  ylab("Binding sites (%) (5' -> 3')")
p


#gd19.5
file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-nk.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pNK <- plotDistToTSS(peakAnno)
pNK <- as.data.frame(pNK$data)

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-sm.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pSM <- plotDistToTSS(peakAnno)
pSM <- as.data.frame(pSM$data)
df <- inner_join(pNK, pSM, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-endo.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pEndo <- plotDistToTSS(peakAnno)
pEndo <- as.data.frame(pEndo$data)
df <- inner_join(df, pEndo, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-macro.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pMacro <- plotDistToTSS(peakAnno)
pMacro <- as.data.frame(pMacro$data)
df <- inner_join(df, pMacro, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

file <- "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-tbc.bed"
peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000), TxDb=txdb)
pTBC <- plotDistToTSS(peakAnno)
pTBC <- as.data.frame(pTBC$data)
df <- inner_join(df, pTBC, by = c("Feature" = "Feature", "sign" = "sign", ".group" = ".group"))

df$mean <- rowMeans(df[,c("freq.x", "freq.y", "freq.x.x", "freq.y.y", "freq")])
df <- df[,c(1,2,4,9)]
df$Feature <- factor(df$Feature, levels=levels(df$Feature))
df <- df[order(-df$.group),]
color <- c("#d73027", "#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4")
names(color) <- unique(df$Feature)
p1 <- ggplot(df, aes(x=1, fill=Feature))
p1 <- p1 + geom_bar(data=subset(df, sign==1), aes(y=mean, fill = Feature), stat="identity") +
  geom_bar(data=subset(df, sign==-1), aes(y=-mean, fill = Feature), stat="identity") +
  scale_fill_manual(values=color)
p1 <- p1 + geom_hline(yintercept = 0, colour = "black") +
  coord_flip() + theme_bw() + ggtitle("Average distance distribution at gd 19.5") +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) + 
  ylab("Binding sites (%) (5' -> 3')")
p1

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure1c.pdf", width = 15, height = 4)
p+p1
dev.off()

