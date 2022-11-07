library(circlize)
library(dplyr)

motifPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaksGD19.5withFinalMotifs.txt", header = T)
motifs19.5 <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor.txt", sep = "\t")
motifPeaks <- inner_join(motifPeaks, motifs19.5[,c("motif", "motif.name")], by = c("motif" = "motif"))
jaspar <- read.table("/work/LAS/geetu-lab/hhvu/JASPAR-library.txt", header = T, sep = "\t")
motifs19.5 <- inner_join(motifs19.5, jaspar, by = c("motif" = "ID"))

dat <- data.frame(matrix(ncol = 24, nrow = 0))
for (i in unique(motifPeaks$motif.name)) {
  sub <- motifPeaks[motifPeaks$motif.name == i,]
  sub2 <- motifPeaks[motifPeaks$peak %in% unique(sub$peak),]
  t <- rep(0, length(unique(motifPeaks$motif.name)))
  names(t) <- unique(motifPeaks$motif.name)
  for (j in unique(motifPeaks$motif.name)) {
    sub3 <- sub2[sub2$motif.name == j,]
    t[j] <- length(unique(sub3$peak))
  }
  
  dat <- rbind(dat, t)
}
rownames(dat) <- unique(motifPeaks$motif.name)
colnames(dat) <- unique(motifPeaks$motif.name)

dat <- as.matrix(dat)

col_fun = colorRamp2(range(dat), c("#ffffbf", "#d7191c"), transparency = 0)
grid.col = c(`TFAP2C` = "red", `TFAP2C(var.3)` = "red", `TFAP2C(var.2)` = "red",
             `CREB3` = "blue", `Creb5` = "blue",
             `JUNB(var.2)` = "grey", `JUN` = "grey",
             `SNAI1` = "pink",
             `FOSB::JUNB(var.2)` = "black", `FOSL2::JUN(var.2)` = "black",
             `FOS::JUN(var.2)` = "black", `FOSL2::JUNB(var.2)` = "black",
             `FOSL1::JUN(var.2)` = "black", `FOSB::JUN` = "black",
             `JUN::JUNB(var.2)` = "brown",
             `TCF4` = "green",
             `NR2F6(var.3)` = "yellow", `Nr2f6` = "yellow",
             `Pparg::Rxra` = "orange",
             `ARNT::HIF1A` = "darkolivegreen", `HIF1A` = "darkgreen",
             `STAT3` = "purple",
             `ETS2` = "grey70",
             `CDX2` = "lightblue"
             )
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/20221024_circosPlotTF.pdf", height = 15, width = 15)
chordDiagram(dat, grid.col = grid.col, col = col_fun, transparency = 0, symmetric = T, preAllocateTracks = 1, annotationTrack = "grid")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

df = data.frame(from = rep(rownames(dat), times = ncol(dat)),
                to = rep(colnames(dat), each = nrow(dat)),
                value = as.vector(dat),
                stringsAsFactors = FALSE)
df2 <- df[df$from != df$to,]
df3 <- df2[!duplicated(apply(df2, 1, function(x) paste0(sort(x), collapse = ""))),]
hist(df3$value, breaks = seq(0, max(df3$value), 1))
df4 <- df3[df3$value >= 50,]


####plot by family

motifPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaksGD19.5withFinalMotifs.txt", header = T)
motifs19.5 <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor.txt", sep = "\t")
motifPeaks <- inner_join(motifPeaks, motifs19.5[,c("motif", "motif.name")], by = c("motif" = "motif"))
jaspar <- read.table("/work/LAS/geetu-lab/hhvu/JASPAR-library.txt", header = T, sep = "\t")
motifPeaks <- inner_join(motifPeaks, jaspar, by = c("motif" = "ID"))

dat <- data.frame(matrix(ncol = length(unique(motifPeaks$Family)), nrow = 0))

for (i in unique(motifPeaks$Family)) {
  sub <- motifPeaks[motifPeaks$Family == i,]
  sub2 <- motifPeaks[motifPeaks$peak %in% unique(sub$peak),]
  t <- rep(0, length(unique(motifPeaks$Family)))
  names(t) <- unique(motifPeaks$Family)
  for (j in unique(motifPeaks$Family)) {
    sub3 <- sub2[sub2$Family == j,]
    t[j] <- length(unique(sub3$peak))
  }
  
  dat <- rbind(dat, t)
}
rownames(dat) <- unique(motifPeaks$Family)
colnames(dat) <- unique(motifPeaks$Family)

dat <- as.matrix(dat)

df = data.frame(from = rep(rownames(dat), times = ncol(dat)),
                to = rep(colnames(dat), each = nrow(dat)),
                value = as.vector(dat),
                stringsAsFactors = FALSE)
df2 <- df[df$from != df$to,]
df3 <- df2[!duplicated(apply(df2, 1, function(x) paste0(sort(x), collapse = ""))),]
hist(df3$value, breaks = seq(0, max(df3$value), 1))
df4 <- df3[df3$value >= 100,]

col_fun = colorRamp2(range(dat), c("#ffffbf", "#d7191c"), transparency = 0.3)

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/20221024_circosPlotTFfam.pdf", height = 15, width = 15)
chordDiagram(dat, col = col_fun, transparency = 0, symmetric = T, link.sort = TRUE, link.decreasing = TRUE, preAllocateTracks = 1, annotationTrack = "grid")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()


###merge family

motifPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaksGD19.5withFinalMotifs.txt", header = T)
motifs19.5 <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor.txt", sep = "\t")
motifPeaks <- inner_join(motifPeaks, motifs19.5[,c("motif", "motif.name")], by = c("motif" = "motif"))
jaspar <- read.table("/work/LAS/geetu-lab/hhvu/JASPAR-library.txt", header = T, sep = "\t")
motifPeaks <- inner_join(motifPeaks, jaspar, by = c("motif" = "ID"))


subColon <- function(x) {
  ifelse(grepl("::", x, ignore.case = TRUE), strsplit(x, "::"), x)
}
keep <- motifPeaks$Family
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
motifPeaks$Family2 <- keep2

#library(Signac)
#library(Seurat)
#library(JASPAR2020)
#library(TFBSTools)
#library(BSgenome.Rnorvegicus.Ensembl.rn6)
#library(patchwork)
#library(dplyr)
#library(org.Rn.eg.db)
#library(clusterProfiler)
#set.seed(1234)
#load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")
#load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-ATACwithRNAlables.rda")
#Idents(rats) <- rats@meta.data$predicted.id # change from old identities to those predicted by scRNA-seq
#DefaultAssay(rats) <- 'MACSpeaks'
#Annotation(rats) <- rnor6.ranges

#peaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5.bed", header=F)
#closest_genes <- ClosestFeature(rats, regions = peaks$V4)
#write.table(closest_genes, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes.txt", sep = "\t", quote = F, row.names = F)

toGene <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes.txt", header = T)
motifPeaks <- inner_join(motifPeaks, toGene, by = c("peak" = "query_region"))
#write.table(motifPeaks, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor-toGenes.txt", sep = "\t", quote = F, row.names = F)

dat <- data.frame(matrix(ncol = length(unique(motifPeaks$Family2)), nrow = 0))

for (i in unique(motifPeaks$Family2)) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  sub2 <- motifPeaks[motifPeaks$peak %in% unique(sub$peak),]
  t <- rep(0, length(unique(motifPeaks$Family2)))
  names(t) <- unique(motifPeaks$Family2)
  for (j in unique(motifPeaks$Family2)) {
    sub3 <- sub2[sub2$Family2 == j,]
    t[j] <- length(unique(sub3$peak))
  }
  
  dat <- rbind(dat, t)
}
rownames(dat) <- unique(motifPeaks$Family2)
colnames(dat) <- unique(motifPeaks$Family2)

dat <- as.matrix(dat)

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/20221024_circosPlotTFfamMerged.pdf", height = 20, width = 20)
chordDiagram(dat, col = col_fun, transparency = 0, symmetric = T, link.sort = TRUE, link.decreasing = TRUE, preAllocateTracks = 1, annotationTrack = "grid")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

df = data.frame(from = rep(rownames(dat), times = ncol(dat)),
                to = rep(colnames(dat), each = nrow(dat)),
                value = as.vector(dat),
                stringsAsFactors = FALSE)
df2 <- df[df$from != df$to,]
df3 <- df2[!duplicated(apply(df2, 1, function(x) paste0(sort(x), collapse = ""))),]
hist(df3$value, breaks = seq(0, max(df3$value), 10))
df4 <- df3[df3$value >= 50,]

#test significance
#number of samples
k <- c()
for (i in df2$from) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  k <- c(k, length(unique(sub$peak)))
}
df2$k <- k

#number of successes
m <- c()
for (i in df2$to) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  m <- c(m, length(unique(sub$peak)))
}
df2$m <- m

#number of population
length(unique(motifPeaks$peak)) #number of peaks with any enriched + expressed motifs

#number of failures
df2$n <- length(unique(motifPeaks$peak)) - df2$m

sig <- phyper(df2$value, df2$m, df2$n, df2$k, lower.tail = FALSE, log.p = FALSE)
df2$pval <- sig
df2$upper.adj.pval <- p.adjust(df2$pval, "BH")
df2$lower.adj.pval <- 1-df2$upper.adj.pval

sub <- df[df$from == df$to, c("from", "to")]
sub$upper.adj.pval <- 1
sub$lower.adj.pval <- 1
r<- df2[,c("from", "to", "upper.adj.pval", "lower.adj.pval")]
r <- rbind(r, sub)
test <- data.frame(matrix(ncol=10, nrow=10))
rownames(test) <- sort(unique(r$from))
colnames(test) <- sort(unique(r$to))
for (i in rownames(test)) {
  for (j in colnames(test)) {
    test[i, j] <- r[r$from == i & r$to == j, "upper.adj.pval"]
  }
}
test <- apply(test, 1, as.numeric)
rownames(test) <- sort(unique(r$from))
test[lower.tri(test)] <- 1-test[lower.tri(test)]
heatmap(t(test), scale="none", Rowv=NA, Colv=NA)

r<- df2[,c("from", "to", "adj.pval")]
r$adj.pval <- -log10(r$adj.pval)
r <- r[!duplicated(apply(r, 1, function(x) paste0(sort(x), collapse = ""))),]
col_fun = colorRamp2(range(r$adj.pval), c("#ffffbf", "#d7191c"), transparency = 0.3)

pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/20221024_circosPlotTFfamMerged_significance.pdf", height = 20, width = 20)
circos.par(gap.after = c(rep(6, 10)))
chordDiagram(r, col = col_fun, transparency = 0, link.sort = TRUE, link.decreasing = TRUE, preAllocateTracks = 1, annotationTrack = "grid")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()