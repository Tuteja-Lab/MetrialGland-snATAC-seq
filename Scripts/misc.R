## Conserved common peak - conserved UNION of gene networks:
Build list of genes that are TBC markers at gd15.5 OR gd19.5 that are conserved in human EVT using scRNA-seq data, and these markers are also near to the common peaks.
```{r}
marker15.5 <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/scRNA-seqFiles/gd15.5-wilcoxin-cluster23.txt", header = T, sep = "\t")
marker19.5 <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/scRNA-seqFiles/gd19.5-wilcoxin-cluster6.txt", header = T, sep = "\t")
markers <- intersect(union(marker15.5$gene, marker19.5$gene), closest_genes$gene_name)
conservedRats <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/conservedRat.txt", header = T, sep = "\t")

markers <- markers[which(toupper(markers) %in% conservedRats$ratGenes == T)]

length(markers)
```

Convert common peaks (gd19.5 coordinates) to human coordinates using LiftOver, then find the regions that have ATAC-seq peaks in EVT. EVT peaks were obtained from https://www.medrxiv.org/content/10.1101/2022.05.25.22275520v1.full-text
```{Bash, eval = F}
cd /work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks
module load bedtools2
bedtools intersect -a commonPeaks-gd19.5-humanLiftOver.bed -b EVTpeaks.bed -wa | sort -u > commonPeaks-gd19.5-humanLiftOver-inEVT.bed
```

Find peaks with enriched motifs and are also conserved
```{r}
conservedPeaks <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/commonPeaks-gd19.5-humanLiftOver-inEVT.bed", header = F)
motifPeaks <- motifPeaks[motifPeaks$peak %in% conservedPeaks$V4,]
```

Associate peaks/motifs to genes
```{r}
network <- inner_join(closest_genes, motifPeaks, by = c("query_region" = "peak"))
network <- inner_join(network, motifs19.5[,c("motif", "motif.name")], by = c("motif" = "motif"))
network <- network[network$gene_name %in% markers,]
dim(network)
```

Grouping motifs of the same gene families together
```{r}
jaspar <- read.table("/work/LAS/geetu-lab/hhvu/JASPAR-library.txt", header = T, sep = "\t")
jaspar <- jaspar[jaspar$ID %in% network$motif,]
network <- inner_join(network, jaspar[,c("ID", "Family")], by = c("motif" = "ID"))
keep <- network$Family
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
network$Family2 <- keep2

all <- data.frame(tfGene=character(), regGene=character())
for (i in unique(network$Family2)) {
  tfGene <- unique(network[which(network$Family2 == i), "motif.name"])
  tf <- toString(tfGene)
  regGene <- unique(network[which(network$Family2 == i), "gene_name"])
  tbl <- data.frame(matrix(nrow = length(tfGene)*length(regGene), ncol = 2))
  for (j in 1:length(tfGene)) {
    tbl[((j-1)*length(regGene)+1):(j*length(regGene)),1] <- tf
    tbl[((j-1)*length(regGene)+1):(j*length(regGene)),2] <- as.character(regGene)
  }
  all <- rbind(all, tbl)
}
all$X2 <- toupper(all$X2)
all <- dplyr::distinct(all)

g <- unique(all$X2)
for (k in 1:length(g)) {
  if (length(grep(g[k], all$X1)) != 0) {
    all[which(g[k] == all$X2), "X2"] <- all[grep(g[k], all$X1)[1], "X1"]
  }
}
all <- dplyr::distinct(all)

#write.table(all, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/conservedCommonPeaksUNIONGenes-gd9.5.txt", sep = "\t", row.names = F, quote = F)
```

Build network in Cytoscape.
```{r}
knitr::include_graphics("Files/conservedPeaksUNIONGenes-gd19.5.PNG")
```

Genes with most associated motifs:
  ```{r}
sort(table(all$X2))[(length(table(all$X2))-5):length(table(all$X2))]
```

```{r}
cov_plot <- CoveragePlot(
  object = unintegrated,
  region = "9-98231011-98232391",
  feature = "Lrrfip1",
  annotation = TRUE,
  peaks = TRUE,
  heights = c(6, 1, 1),
  extend.upstream = 7000,
  extend.downstream = 7000,
  idents = c("GD15.5_Trophoblast_cells", "GD19.5_Trophoblast_cells")
)

cov_plot
```



```{r}
cov_plot <- CoveragePlot(
  object = unintegrated,
  region = "4-147273825-147275519",
  feature = "Pparg",
  annotation = TRUE,
  peaks = TRUE,
  heights = c(6, 1, 1),
  extend.upstream = 7000,
  extend.downstream = 7000,
  idents = c("GD15.5_Trophoblast_cells", "GD19.5_Trophoblast_cells")
)

cov_plot
```

```{r}
t2 <- closest_genes[closest_genes$gene_name %in% df[df$category == ">= 3 peaks", "gene"],]
#go(unique(t2$ensembl_gene_id)) #no terms enriched

gd15.5tbc <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/res0.8/gd15.5-wilcoxin-cluster23.txt", header = T)
length(intersect(df[df$category == ">= 3 peaks", "gene"], gd15.5tbc$gene))
length(unique(df[df$category == ">= 3 peaks", "gene"]))
intersect(df[df$category == ">= 3 peaks", "gene"], gd15.5tbc$gene) #57 genes out of 255 genes are TBC genes
genes <- names(t[t >= 4])
length(intersect(genes, gd15.5tbc$gene)) #25 out of 57 genes have >= 4 peaks
genes <- names(t[t >= 5])
length(genes)
length(intersect(genes, gd15.5tbc$gene)) #16 out of 57 genes have >= 5 peaks
```


library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Rnorvegicus.Ensembl.rn6)
library(patchwork)
library(dplyr)
library(org.Rn.eg.db)
library(clusterProfiler)
library(chromVAR)
set.seed(1234)

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/6_networks/ccans-gd19.5.rda")
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/6_networks/conns-gd19.5.rda")

t <- conns[conns$Peak1 == "1-12808761-12809434" | conns$Peak2 == "1-12808761-12809434",]
t <- t[t$coaccess > 0,]

enriched.motifs <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/GD19.5_motifAnalysis/common19.5-ratsMotifs.txt", header = T, sep = "\t")

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/GD19.5_motifAnalysis/gd19.5_ATACwithMousePFM.rda")
rats <- RunChromVAR(rats, genome = BSgenome.Rnorvegicus.Ensembl.rn6)
rats2 <- rats

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/4_RNAintegration/2_MACS2Peaks/annotation_scRNAseq/gd19.5-ATACwithRNAlables.rda")

chromvar.data <- GetAssayData(object =  rats2, slot = 'data')
rats[["chromvar"]] <- CreateAssayObject(data = chromvar.data)
DefaultAssay(rats) <- 'chromvar'


FeaturePlot(
  object = rats,
  features = "MA0524.2",
  pt.size = 0.1,
  label = T,
  repel = T)
FeaturePlot(
  object = rats,
  features = "MA0814.2",
  pt.size = 0.1)
FeaturePlot(
  object = rats,
  features = "MA0815.1",
  pt.size = 0.1)


aver_chromvar <- AverageExpression(rats, assays = "chromvar", features = enriched.motifs$motif[1:10]) %>%
  as.data.frame()

# sort the heatmap prior to plotting find column index for maxium value for each TF activity
aver_chromvar <- aver_chromvar[do.call(order, c(aver_chromvar, list(decreasing=TRUE))),]
aver_chromvar$max <- max.col(aver_chromvar)
aver_chromvar <- aver_chromvar[order(aver_chromvar$max), ]
aver_chromvar <- dplyr::select(aver_chromvar, -max)
colnames(aver_chromvar) <- levels(Idents(rats))
pheatmap::pheatmap(aver_chromvar,scale = "row",
                   cluster_cols=F,cluster_rows = F,
                   show_rownames=T) #450x640

# gather the footprinting information for sets of motifs
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/5_integrating2timepoints/annotation_scRNAseq/commonPeaks/GD19.5_motifAnalysis/gd19.5_ATACwithHumanPFM.rda")
peaks.1based <- GRangesToString(
  grange = StringToGRanges(
    regions = rownames(rats),
    sep = c(":","-"),
    starts.in.df.are.0based = TRUE
  ),
  sep = c(":", "-")
)
rownames(rats@assays[["MACSpeaks"]]@CountS) <- ATAC.1based
rownames(rats@assays[["MACSpeaks"]]@DaTa) <- ATAC.1based
rownames(rats@assays[["MACSpeaks"]]@meta.features) <- ATAC.1based

gorb <- StringToGRanges(ATAC.1based, sep = c("-", "-"))
aggr@assays[["peaks"]]@ranges <- gorb
rats <- Footprint(
  object = rats,
  motif.name = c("TFAP2C"),
  in.peaks = TRUE,
  genome = BSgenome.Rnorvegicus.Ensembl.rn6
)

# plot the footprint data for each group of cells
PlotFootprint(rats, features = c("TFAP2C"))