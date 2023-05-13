library(ggplot2)

#figure 3A after revision
r <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes-Webgestal/enrichment_results_wg_result1683039055.txt", header = T, sep = "\t")
r <- r[r$enrichmentRatio >= 1.5 & r$FDR <= 0.05 & r$overlap >= 5,]
#write.table(r, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes-Webgestal/goFiltered.txt", quote = F, sep = "\t", row.names = F)

r <- r[r$description %in% c("cell-cell adhesion",
                            "positive regulation of cell migration",
                            "female pregnancy",
                            "canonical Wnt signaling pathway",
                            "positive regulation of JAK-STAT cascade"),]
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/commPeaks-GO.pdf", width = 15, height = 5)
ggplot(data=r, aes(x=-log10(FDR), y=description,
                   color=log2(`enrichmentRatio`), size=`overlap`)) + 
  geom_point(show.legend = TRUE) + ylab("") + xlab("-Log10(FDR)") + scale_color_gradient(low = "#fb6a4a", high = "#a50f15") +
  scale_size_continuous(range = c(5, 20)) +
  theme(legend.text = element_text(size=10), legend.title = element_text(size=10),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)) + ggtitle("Interesting GO terms enriched with genes associated to common peaks")
dev.off()


#TBC.markers <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/daPeaks/TBC-unfilteredDApeaks.txt", header = T)
#names(TBC.markers)[names(TBC.markers) == "p_val_adj"] <- 'padj'
#names(TBC.markers)[names(TBC.markers) == "avg_log2FC"] <- 'log2FoldChange'
#TBC.markers$log2FoldChange <- -TBC.markers$log2FoldChange

#TBC.markers$threshold <- ifelse(TBC.markers$padj <= 0.05
#                                & TBC.markers$log2FoldChange >= log2(1.5), "#7b3294",
#                                ifelse(TBC.markers$padj <= 0.05
#                                       & TBC.markers$log2FoldChange <= -log2(1.5), "#008837", "grey70"))

#input <- cbind(TBC.markers, gene=rownames(TBC.markers))
#ggplot(TBC.markers) +
#  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
#  ggtitle("gd 15.5 vs gd 19.5 invasive trophoblast cell cluster") +
#  xlab("log2(Fold Change)") + 
#  ylab("-log10(adjusted p-value)") +
#  scale_colour_identity("Legend", labels = c("gd 15.5", "Insignificant", "gd 19.5"),
#                        breaks = c("#008837", "grey70", "#7b3294"), guide = "legend") +
#  #scale_y_continuous(limits = c(0,50)) +
#  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
#  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "black", size = 1) +
#  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "black", size = 1) +
#  theme(#legend.position = "right", 
#    plot.title = element_text(size = 20, face = "bold"),
#    axis.title = element_text(size = 20),
#    axis.text = element_text(size = 20),
#    legend.title = element_text(size=20),
#    legend.text = element_text(size=20))

#figure 3a - motifs
motifPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaksGD19.5withFinalMotifs.txt", header = T)
motifs19.5 <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor.txt", sep = "\t")
motifPeaks <- inner_join(motifPeaks, motifs19.5[,c("motif", "motif.name")], by = c("motif" = "motif"))
jaspar <- read.table("/work/LAS/geetu-lab/hhvu/JASPAR-library.txt", header = T, sep = "\t")
motifs19.5 <- inner_join(motifs19.5, jaspar, by = c("motif" = "ID"))

subColon <- function(x) {
  ifelse(grepl("::", x, ignore.case = TRUE), strsplit(x, "::"), x)
}
keep <- motifs19.5$Family
keep <- subColon(keep)
for (i in 1:length(keep)) {
  keep[i] <- list(unique(keep[i][[1]]))
  for (j in setdiff(1:length(keep), which(keep %in% keep[i]))) {
    if (sum(keep[j][[1]] %in% keep[i][[1]]) > 0) {
      keep[j] <- keep[i]
    }
  }
}
keep2 <- unlist(lapply(keep, toString))
motifs19.5$Family2 <- keep2

r <- data.frame(matrix(ncol=ncol(motifs19.5), nrow = 0))
for (i in unique(keep2)) {
  sub <- motifs19.5[motifs19.5$Family2 == i,]
  t <- sub[sub$fold.enrichment == max(sub$fold.enrichment),]
  r <- rbind(r, t)
}
r <- r[-3,]

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/Fig3a-motifs.pdf", width = 8, height = 10)
ggplot(data=r, aes(x=reorder(motif.name, log2(fold.enrichment)), y=log2(fold.enrichment))) +
  geom_bar(stat="identity", fill="#D40000", width = 0.7) + 
  coord_flip() +
  labs(title = "Motifs of enriched TF families",
       y = "log2(fold change)",
       x = "") +
  theme(plot.title = element_text(size = 25, face = "bold"), legend.text=element_text(size=25),
        axis.text=element_text(size=25), axis.title.x = element_text(size = 25))
dev.off()

#figure 3b - heatmap
#see 5_TFcombo.R


