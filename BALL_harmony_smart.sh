#!/bin/bash
echo "🚀 Starting Unified Harmony Pipeline with Smart Checkpoints..."

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.3"

cat > BALL_harmony_smart.R << 'EOF'
library(Seurat)
library(harmony)
library(future)

plan("multicore", workers = 64)
options(future.globals.maxSize = 1024 * 1024^3)

dir.create("FastHarmony_Results", showWarnings = FALSE)
dir.create("FastHarmony_Results/checkpoints", showWarnings = FALSE)

log <- function(msg) {
  message(sprintf("[%s] %s", format(Sys.time(), "%F %T"), msg))
}

samples <- list.dirs("/data/BALL_processed/BALL_scRNA_matrix", full.names = FALSE, recursive = FALSE)

# Step 1: Load or resume checkpointed objects
load_or_preprocess <- function(id) {
  path <- file.path("FastHarmony_Results/checkpoints", paste0(id, "_pca.rds"))
  if (file.exists(path)) {
    log(paste("🔁 Resuming:", id))
    return(readRDS(path))
  }
  tryCatch({
    dir <- file.path("/data/BALL_processed/BALL_scRNA_matrix", id)
    genes <- read.table(file.path(dir, "genes.tsv"), sep = "\t", stringsAsFactors = FALSE)
    colnames(genes) <- c("ENSEMBL", "SYMBOL")
    genes$SYMBOL <- ifelse(duplicated(genes$SYMBOL) | duplicated(genes$SYMBOL, fromLast = TRUE),
                           paste0(genes$SYMBOL, "_", genes$ENSEMBL),
                           genes$SYMBOL)
    counts <- Read10X(dir, gene.column = 2)
    rownames(counts) <- genes$SYMBOL
    obj <- CreateSeuratObject(counts, project = id, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)

    log(paste("⚙️ Preprocessing:", id))
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, npcs = 50)
    saveRDS(obj, path)
    obj
  }, error = function(e) {
    log(paste("⚠️ Error loading", id, ":", conditionMessage(e)))
    NULL
  })
}

objs <- future.apply::future_lapply(samples, load_or_preprocess)
objs <- objs[!sapply(objs, is.null)]
log(sprintf("✅ Loaded %d total samples", length(objs)))

# Step 2: Merge and re-run PCA
log("🔗 Merging all preprocessed objects...")
combined <- merge(x = objs[[1]], y = objs[-1])
DefaultAssay(combined) <- "RNA"
combined <- RunPCA(combined, npcs = 50)
saveRDS(combined, "FastHarmony_Results/checkpoints/combined_pca.rds")

# Step 3: Harmony batch correction
log("🧠 Running Harmony...")
combined$orig.ident <- as.factor(combined$orig.ident)
combined <- RunHarmony(combined, group.by.vars = "orig.ident", max_iter = 20)
saveRDS(combined, "FastHarmony_Results/checkpoints/harmony_corrected.rds")

# Step 4: Clustering and visualization
log("📊 UMAP and clustering...")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = c(0.3, 0.6, 0.9))

saveRDS(combined, "FastHarmony_Results/harmony_clustered_data.rds")
write.csv(combined@meta.data, "FastHarmony_Results/metadata.csv")

pdf("FastHarmony_Results/clustering_plots.pdf", width = 14)
DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + ggtitle("Harmony Clustering")
VlnPlot(combined, features = c("nFeature_RNA", "percent.mt"), pt.size = 0)
dev.off()

log("✅ Pipeline complete. Final outputs saved to FastHarmony_Results/")
EOF

nohup Rscript BALL_harmony_smart.R > harmony_smart.log 2>&1 &
echo "🧬 Monitor with: tail -f harmony_smart.log"

