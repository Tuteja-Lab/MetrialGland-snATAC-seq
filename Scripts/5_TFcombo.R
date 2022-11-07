library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


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
#set.seed(1234)
#load("/work/LAS/geetu-lab/hhvu/project3_scATAC/rnor6.ranges.rda")
#load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-ATACwithRNAlables.rda")
#Idents(rats) <- rats@meta.data$predicted.id # change from old identities to those predicted by scRNA-seq
#DefaultAssay(rats) <- 'MACSpeaks'
#Annotation(rats) <- rnor6.ranges

#peaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5.bed", header=F)
#closest_genes <- ClosestFeature(rats, regions = peaks$V4)
#write.table(closest_genes, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes.txt", sep = "\t", quote = F, row.names = F)

#toGene <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes.txt", header = T)
#motifPeaks <- inner_join(motifPeaks, toGene, by = c("peak" = "query_region"))
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

df = data.frame(from = rep(rownames(dat), times = ncol(dat)),
                to = rep(colnames(dat), each = nrow(dat)),
                value = as.vector(dat),
                stringsAsFactors = FALSE)
df2 <- df[df$from != df$to,]
#df3 <- df2[!duplicated(apply(df2, 1, function(x) paste0(sort(x), collapse = ""))),]
#hist(df3$value, breaks = seq(0, max(df3$value), 10))
#df4 <- df3[df3$value >= 50,]

#test significance
#number of samples
k <- c()
for (i in df2$from) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  k <- c(k, length(unique(sub$peak)))
}
df2$k <- k #number of peaks with binding sites of family A

#number of successes
m <- c()
for (i in df2$to) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  m <- c(m, length(unique(sub$peak)))
}
df2$m <- m #number of peaks with binding sites of family B

#number of failures
df2$n <- 1242 - df2$m #number of peaks without binding sites of family B

#df2$value = number of peaks with binding sites of both families

sig <- phyper(df2$value, df2$m, df2$n, df2$k, lower.tail = FALSE, log.p = FALSE)
df2$pval <- sig
df2$upper.adj.pval <- p.adjust(df2$pval, "BH")
sig <- phyper(df2$value, df2$m, df2$n, df2$k, lower.tail = TRUE, log.p = FALSE)
df2$lower.adj.pval <- p.adjust(sig, "BH")

test <- data.frame(matrix(ncol=10, nrow=10))
rownames(test) <- sort(unique(df2$from))
colnames(test) <- sort(unique(df2$to))
for (i in rownames(test)) {
  for (j in colnames(test)) {
    if (i != j) {
      test[i, j] <- df2[df2$from == i & df2$to == j, "upper.adj.pval"] 
    } else {
      test[i, j] <- NA
    }
  }
}
test <- apply(test, 1, as.numeric)
rownames(test) <- sort(unique(df2$from))


test2 <- data.frame(matrix(ncol=10, nrow=10))
rownames(test2) <- sort(unique(df2$from))
colnames(test2) <- sort(unique(df2$to))
for (i in rownames(test2)) {
  for (j in colnames(test2)) {
    if (i != j) {
      test2[i, j] <- df2[df2$from == i & df2$to == j, "lower.adj.pval"] 
    } else {
      test2[i, j] <- NA
    }
  }
}
test2 <- apply(test2, 1, as.numeric)
rownames(test2) <- sort(unique(df2$from))
test[lower.tri(test)] <- test2[lower.tri(test2)]
#heatmap(test, scale="none", Rowv=NA, Colv=NA)

test2 <- reshape2::melt(test)

col1 = circlize::colorRamp2(c(1, 0), c("#fcae91", "#a50f15"))
col2 = circlize::colorRamp2(c(1, 0), c("#eff3ff", "#084594"))

# here reordering the symmetric matrix is necessary
ht = Heatmap(test, rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE,
             cluster_rows = FALSE, cluster_columns = FALSE,
             layer_fun = function(j, i, x, y, w, h, fill) {
               l = i > j
               grid.rect(x[l], y[l], w[l], h[l], 
                         gp = gpar(fill = col2(pindex(test, i[l], j[l])), col = NA))
               
               l = i < j
               grid.rect(x[l], y[l], w[l], h[l], 
                         gp = gpar(fill = col1(pindex(test, i[l], j[l])), col = NA))
               
               l = i == j
               grid.rect(x[l], y[l], w[l], h[l], 
                         gp = gpar(fill = "grey", col = NA))
               
               v = pindex(test, i, j)
               grid.text(sprintf("%.3f", v), x, y, gp = gpar(fontsize = 15))
             })

pdf(file="/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/fig3b-heatmap2.pdf", height = 10, width = 11)
draw(ht, heatmap_legend_list = list(
  Legend(title = "Over-representation", col_fun = col1),
  Legend(title = "Under-representation", col_fun = col2)
))
dev.off()

motifPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor-toGenes.txt", header = T, sep="\t")
fos <- motifPeaks[motifPeaks$Family2 == "Fos-related factors, Jun-related factors",]
ap2 <- motifPeaks[motifPeaks$Family2 == "AP-2",]
length(unique(intersect(fos$gene_name, ap2$gene_name)))
length(unique(fos$gene_name))
length(unique(ap2$gene_name))
phyper(q=length(unique(intersect(fos$gene_name, ap2$gene_name))),
       m=length(unique(fos$gene_name)),
       n=length(unique(motifPeaks$gene_name)) - length(unique(fos$gene_name)),
       k=length(unique(ap2$gene_name)), lower.tail = F)
phyper(q=length(unique(intersect(fos$gene_name, ap2$gene_name))),
       m=length(unique(ap2$gene_name)),
       n=length(unique(motifPeaks$gene_name)) - length(unique(ap2$gene_name)),
       k=length(unique(fos$gene_name)), lower.tail = F)
