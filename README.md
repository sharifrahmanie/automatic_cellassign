# automatic_cellassign
Automatic cell type assingning for Single-cell RNA-Seq dataset by provideing count matrix and marker list
```r{}
require(scran)
require(igraph)
require(SingleCellExperiment)
require(AUCell)
require(reshape2)

automatic_cellassign <- function(counts, 
                                  markerlist) {
  scexp <- SingleCellExperiment(assays = list(counts = counts))
  scexp$cell <- colnames(scexp)
  scexp <- computeSumFactors(scexp, min.mean = 0.1)
  scexp <- logNormCounts(scexp)
  logcounts(scexp) <- as.matrix(logcounts(scexp))
  dec <- modelGeneVar(scexp)
  scexp <- denoisePCA(scexp, technical = dec, BSPARAM = IrlbaParam())
  scexp <- runTSNE(scexp, dimred = "PCA", perplexity = 30)
  cell_rankings <- AUCell_buildRankings(logcounts(scexp))
  cell_AUC <- AUCell_calcAUC(markerlist, cell_rankings)
  cell_assignment <- AUCell::AUCell_exploreThresholds(
    cell_AUC, plotHist = FALSE, assignCells = TRUE)
  cellsAssigned <- lapply(cell_assignment, function(x) x$assignment)
  new_cellsAssigned <- list()
  for (i in seq_along(cellsAssigned)) {
    if(length(cellsAssigned[[i]]) > 0) {
      new_cellsAssigned[[names(cellsAssigned)[i]]] <- cellsAssigned[[i]]
    }
  }
  assignmentTable <- melt(new_cellsAssigned, value.name = "cell")
  colnames(assignmentTable)[2] <- "geneSet"
  assignments <- assignmentTable %>% dplyr::group_by(cell) %>% 
    summarize(AUCellType = paste(geneSet, collapse = "/")) %>%
    group_by(AUCellType) %>%
    mutate(n = length(unique(cell))) %>%
    filter(n > 5)
  scexp$CellTypes <- assignments$AUCellType[match(scexp$cell, assignments$cell)]
  plotTSNE(scexp, colour_by = "CellTypes")
}

load("markerlist.RData")
load("pbmc3k_filtered.RData")
automatic_cellassign(counts = pbmc3k_filtered,
                      markerlist = markerlist)

```
