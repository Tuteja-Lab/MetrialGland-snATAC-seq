library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


###merge family
toGenes <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-toGenes.txt", header = T, sep = "\t")
motifPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/finalMotifs-commonPeaksGD19.5coor-toGenes.txt", header = T, sep = "\t")

df = data.frame(from = rep(unique(motifPeaks$Family2), times = 10),
                to = rep(unique(motifPeaks$Family2), each = 10),
                q = NA,
                stringsAsFactors = FALSE)

#test significance
#number of samples
k <- c()
for (i in df$from) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  k <- c(k, length(unique(sub$gene_name)))
}
df$k <- k #number of genes regulated by family A

#number of successes
m <- c()
for (i in df$to) {
  sub <- motifPeaks[motifPeaks$Family2 == i,]
  m <- c(m, length(unique(sub$gene_name)))
}
df$m <- m #number of genes regulated by family B

#number of failures
df$n <- length(unique(toGenes$gene_name)) - df$m #number of peaks without binding sites of family B

#number of genes regulated by both families
q <- c()
for (i in unique(df$to)) {
  for (j in unique(df$from)) {
    sub1 <- motifPeaks[motifPeaks$Family2 == i,]
    sub2 <- motifPeaks[motifPeaks$Family2 == j,]
    q <- c(q, length(intersect(sub1$gene_name, sub2$gene_name)))
  }
}
df$q <- q #number of genes regulated by both families

df <- df[df$from != df$to,]
sig <- phyper(df$q, df$m, df$n, df$k, lower.tail = FALSE, log.p = FALSE)
df$upper.adj.pval <- p.adjust(sig, "BH")
sig <- phyper(df$q, df$m, df$n, df$k, lower.tail = TRUE, log.p = FALSE)
df$lower.adj.pval <- p.adjust(sig, "BH")

test <- data.frame(matrix(ncol=10, nrow=10))
rownames(test) <- sort(unique(df$from))
colnames(test) <- sort(unique(df$to))
for (i in rownames(test)) {
  for (j in colnames(test)) {
    if (i != j) {
      test[i, j] <- df[df$from == i & df$to == j, "upper.adj.pval"] 
    } else {
      test[i, j] <- NA
    }
  }
}
test <- apply(test, 1, as.numeric)
rownames(test) <- sort(unique(df$from))


test2 <- data.frame(matrix(ncol=10, nrow=10))
rownames(test2) <- sort(unique(df$from))
colnames(test2) <- sort(unique(df$to))
for (i in rownames(test2)) {
  for (j in colnames(test2)) {
    if (i != j) {
      test2[i, j] <- df[df$from == i & df$to == j, "lower.adj.pval"] 
    } else {
      test2[i, j] <- NA
    }
  }
}
test2 <- apply(test2, 1, as.numeric)
rownames(test2) <- sort(unique(df$from))
test[lower.tri(test)] <- test2[lower.tri(test2)]

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

pdf(file="/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/FIGURES/20220829_draftV3/fig3b-heatmap2-genes.pdf", height = 10, width = 11)
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
