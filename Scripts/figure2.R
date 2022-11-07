library(Seurat)
library(Signac)
library(Seurat)
library(S4Vectors)
library(GenomicRanges)
library(patchwork)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(ensembldb)
library(biovizBase)
library(EnsDb.Rnorvegicus.v98)
library(AnnotationHub)
library(AnnotationForge)
library(ggsignif)
set.seed(1234)



#figure 2A - peak number distribution
specificPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-tbc-toGenes.txt", sep = "\t")
t <- as.data.frame(table(specificPeaks$gene_name))
p1<-ggplot(t, aes(x = Freq)) + geom_histogram(bins = 9, fill = "#c6dbef", color="grey") +
  scale_x_continuous(breaks=c(0:9)) +
  stat_bin(aes(y=..count.., label=ifelse(..count..==0,"",..count..)), geom="text", vjust=-.5) +
  theme_bw() +
  ggtitle("Distribution of peak number per gene - gd 15.5")


specificPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-tbc-toGenes.txt", sep = "\t")
t <- as.data.frame(table(specificPeaks$gene_name))
p2<-ggplot(t, aes(x = Freq)) + geom_histogram(bins = 7, fill = "#08519c", color="grey") +
  scale_x_continuous(breaks=c(0:7)) +
  stat_bin(aes(y=..count.., label=ifelse(..count..==0,"",..count..)), geom="text", vjust=-.5) +
  theme_bw() +
  ggtitle("Distribution of peak number per gene - gd 19.5")

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2A-peakDis.pdf", width = 8, height = 6)
p1+p2
dev.off()

#figure 2B - high number of peaks, high expression
#gd15.5
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/scRNA-seqFiles/gd15.5_res0.8.rda")
exp <- AverageExpression(data)
exp <- as.data.frame(exp$RNA)
avgRNA <- exp
specificPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-tbc-toGenes.txt", sep = "\t")
t <- table(specificPeaks$gene_name)
genes <- names(t[t >= 3])
exp <- exp[order(rownames(exp)),]
exp$gene <- rownames(exp)
df <- exp[exp$gene %in% names(t), c("gene", "23")]
colnames(df) <- c("gene", "value")
df2 <- data.frame(gene = setdiff(names(t), df$gene), value = rep(0, length(setdiff(names(t), df$gene))))
names(df2) <- colnames(df)
df <- rbind(df, df2)
df$category <- ifelse(df$gene %in% genes, ">= 3 peaks", "< 3 peaks")
p1<-ggplot(df, aes(x=category, y=log10(value+10^(-5)))) +
  geom_boxplot(fill = c("#c6dbef", "#08519c")) +
  ggtitle("Cutoff = 3, TBC peaks all genes GD15.5") +
  geom_signif(stat="identity",
              data=data.frame(x=c(1), xend=c(2),
                              y=c(2.3), annotation=c("*")),
              aes(x=x, xend=xend, y=y, yend=y, annotation=annotation),
              tip_length = 2) +theme_bw()

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/scRNA-seqFiles/gd19.5-4-5-6-7_res0.8.rda")
exp <- AverageExpression(data2)
exp <- as.data.frame(exp$RNA)
avgRNA <- exp
data <- data2
specificPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-tbc-toGenes.txt", sep = "\t")
t <- table(specificPeaks$gene_name)
genes <- names(t[t >= 3])
exp <- exp[order(rownames(exp)),]
exp$gene <- rownames(exp)
df <- exp[exp$gene %in% names(t), c("gene", "6")]
colnames(df) <- c("gene", "value")
df2 <- data.frame(gene = setdiff(names(t), df$gene), value = rep(0, length(setdiff(names(t), df$gene))))
names(df2) <- colnames(df)
df <- rbind(df, df2)
df$category <- ifelse(df$gene %in% genes, ">= 3 peaks", "< 3 peaks")
p2<-ggplot(df, aes(x=category, y=log10(value+10^(-5)))) +
  geom_boxplot(fill = c("#c6dbef", "#08519c")) +
  ggtitle("Cutoff = 3, TBC peaks all genes GD19.5")  +
  geom_signif(stat="identity",
              data=data.frame(x=c(1), xend=c(2),
                              y=c(2.8), annotation=c("*")),
              aes(x=x, xend=xend, y=y, yend=y, annotation=annotation),
              tip_length = 2) +theme_bw()

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2B-expression.pdf", width = 6, height = 6)
p1+p2
dev.off()


#figure 2C - example of genes with many peaks
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-ATACwithRNAlables.rda")

Idents(rats) <- rats@meta.data$predicted.id # change from old identities to those predicted by scRNA-seq
DefaultAssay(rats) <- 'MACSpeaks'

annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v98))
annotations <- subset(annotations, gene_biotype=="protein_coding")
seqlevelsStyle(annotations) <- "Ensembl"
Annotation(rats) <- annotations

cols <- c("orange3", "#718748", "#183D61", "#79704C", "red3")

cov_plot <- CoveragePlot( 
  object = rats,
  region = "1-12823363-12825786",
  feature = "Cited2",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 200000,
  extend.downstream = 10000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Cited2-15.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 1 position")) + ggtitle("GD15.5 Cited2") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()

cov_plot <- CoveragePlot( 
  object = rats,
  region = "11-36075709-36092495",
  feature = "Ets2",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 50000,
  extend.downstream = 120000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Ets2-15.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 11 position")) + ggtitle("GD15.5 Ets2") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()

cov_plot <- CoveragePlot( 
  object = rats,
  region = "3-170550314-170558194",
  feature = "Tfap2c",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 50000,
  extend.downstream = 120000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Tfap2c-15.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 3 position")) + ggtitle("GD15.5 Tfap2c") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()

#Prl5a1 17-38323674-38330944
cov_plot <- CoveragePlot( 
  object = rats,
  region = "17-38323674-38330944",
  feature = "Prl5a1",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 300000,
  extend.downstream = 10000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Prl5a1-15.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 17 position")) + ggtitle("GD15.5 Prl5a1") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()


#gd19.5
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-ATACwithRNAlables.rda")
Idents(rats) <- rats@meta.data$predicted.id # change from old identities to those predicted by scRNA-seq
DefaultAssay(rats) <- 'MACSpeaks'
Annotation(rats) <- annotations
Idents(rats) <- factor(Idents(rats), levels = c("Other cells", "Natural killer cells", "Endothelial cells", "Macrophage cells",         
                                                "Smooth muscle cells", "Invasive trophoblast cells"))
cols <- c("orange3", "#718748", "#183D61", "#79704C", "red3")

cov_plot <- CoveragePlot( 
  object = rats,
  region = "1-12823363-12825786",
  feature = "Cited2",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 215000,
  extend.downstream = 10000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Cited2-19.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 1 position")) + ggtitle("GD19.5 Cited2") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()

cov_plot <- CoveragePlot( 
  object = rats,
  region = "11-36075709-36092495",
  feature = "Ets2",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 70000,
  extend.downstream = 120000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Ets2-19.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 11 position")) + ggtitle("GD19.5 Ets2") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()

cov_plot <- CoveragePlot( 
  object = rats,
  region = "3-170550314-170558194",
  feature = "Tfap2c",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 40000,
  extend.downstream = 120000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Tfap2c-19.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 3 position")) + ggtitle("GD19.5 Tfap2c") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()

#Prl5a1 17-38323674-38330944
cov_plot <- CoveragePlot( 
  object = rats,
  region = "17-38323674-38330944",
  feature = "Prl5a1",
  annotation = T,
  peaks = TRUE,
  heights = c(6, 2, 2),
  extend.upstream = 300000,
  extend.downstream = 10000,
  idents = c("Natural killer cells", "Endothelial cells", "Macrophage cells",         
             "Smooth muscle cells", "Invasive trophoblast cells")
)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2C-Prl5a1-19.5.pdf", width = 8, height = 5)
wrap_elements(cov_plot & scale_fill_manual(values = cols) & xlab("Chromosome 17 position")) + ggtitle("GD19.5 Prl5a1") + theme(plot.title = element_text(size=18, face = "bold")) 
dev.off()



#library(stringr)
#test <- rownames(atac3)
#test <- test[grep("^17-", test, value = F)]
#test <- as.data.frame(test)
#test[c('chr', 'start', 'stop')] <- str_split_fixed(test$test, '-', 3)
#39153133-39153572
#test2 <- test[test$start > 39150000 & test$stop < 39270000,]

#figure 2D -  number of EVT peaks
#gd15.5
specificPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-tbc-toGenes.txt", sep = "\t")
t <- table(specificPeaks$gene_name)
library(stringr)
geneToPeak <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/EVT_GREAT/EVT_genesToPeaks.txt", header = F, sep = "\t", fill = T, stringsAsFactors = F, quote = "")
geneToPeak$noPeaks <- str_count(geneToPeak$V2, "EVT_peaks")
quantile(geneToPeak$noPeaks)
closest_genes <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd15.5-tbc-toGenes.txt", sep = "\t")
closest_genes$gene_name <- toupper(closest_genes$gene_name)

load(file = "/work/LAS/geetu-lab/hhvu/combine-test-expression1.Rdata")
ratHumanOrthologs <- dataset$GRCH38$ratHumanOrthologs
closest_genes2 <- inner_join(closest_genes, ratHumanOrthologs[,c(2,4)], by = c("gene_name" = "Gene.name"))
t <- table(closest_genes2$Human.gene.name)
t <- as.data.frame(t)
t$Var1 <- toupper(t$Var1)

geneToPeak3 <- geneToPeak[geneToPeak$V1 %in% closest_genes2$Human.gene.name,]
geneToPeak3 <- inner_join(geneToPeak3, t, by = c("V1" = "Var1"))
geneToPeak3$category <- ifelse(geneToPeak3$Freq >= 3, ">= 3 rats", "< 3 rats")

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2D-evtPeaks.pdf", width = 6, height = 6)
ggplot(geneToPeak3, aes(x=category, y=noPeaks)) +
  geom_boxplot(fill = c("#c6dbef", "#08519c")) +
  ggtitle("Cutoff = 3, invasive trophoblast GD15.5") +
  geom_signif(stat="identity",
              data=data.frame(x=c(1), xend=c(2),
                              y=c(120), annotation=c("*")),
              aes(x=x, xend=xend, y=y, yend=y, annotation=annotation),
              tip_length = 2) +theme_bw()
dev.off()

#19.5
specificPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-tbc-toGenes.txt", sep = "\t")
t <- table(specificPeaks$gene_name)
library(stringr)
geneToPeak <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/EVT_GREAT/EVT_genesToPeaks.txt", header = F, sep = "\t", fill = T, stringsAsFactors = F, quote = "")
geneToPeak$noPeaks <- str_count(geneToPeak$V2, "EVT_peaks")

closest_genes <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-tbc-toGenes.txt", sep = "\t")
closest_genes$gene_name <- toupper(closest_genes$gene_name)

load(file = "/work/LAS/geetu-lab/hhvu/combine-test-expression1.Rdata")
ratHumanOrthologs <- dataset$GRCH38$ratHumanOrthologs
closest_genes2 <- inner_join(closest_genes, ratHumanOrthologs[,c(2,4)], by = c("gene_name" = "Gene.name"))
t <- table(closest_genes2$Human.gene.name)
t <- as.data.frame(t)
t$Var1 <- toupper(t$Var1)

geneToPeak3 <- geneToPeak[geneToPeak$V1 %in% closest_genes2$Human.gene.name,]
geneToPeak3 <- inner_join(geneToPeak3, t, by = c("V1" = "Var1"))
geneToPeak3$category <- ifelse(geneToPeak3$Freq >= 3, ">= 3 rats", "< 3 rats")

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/figure2D-evtPeaks-19.5.pdf", width = 6, height = 6)
ggplot(geneToPeak3, aes(x=category, y=noPeaks)) +
  geom_boxplot(fill = c("#c6dbef", "#08519c")) +
  ggtitle("Cutoff = 3, invasive trophoblast GD19.5") +
  geom_signif(stat="identity",
              data=data.frame(x=c(1), xend=c(2),
                              y=c(120), annotation=c("*")),
              aes(x=x, xend=xend, y=y, yend=y, annotation=annotation),
              tip_length = 2) +theme_bw()
dev.off()
