#!/bin/bash
echo "🚀 Launching Final RPCA Integration Pipeline..."

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.3"

cat > BALL_seurat_final.R << 'EOF'
library(Seurat)
library(harmony)
library(ggplot2)
library(future)

plan("multicore", workers = 64)
options(future.globals.maxSize = 1024 * 1024^3)

log <- function(msg) {
  message(sprintf("[%s] %s", format(Sys.time(), "%F %T"), msg))
}

dir.create("Core_Results", showWarnings = FALSE)

# Step 1: Load and filter samples
load_sample <- function(id) {
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
    subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
  }, error = function(e) {
    log(paste("⚠️ Failed loading", id, ":", conditionMessage(e)))
    NULL
  })
}

samples <- list.dirs("/data/BALL_processed/BALL_scRNA_matrix", full.names = FALSE, recursive = FALSE)
ball.list <- future.apply::future_lapply(samples, load_sample, future.seed = TRUE)
ball.list <- ball.list[!sapply(ball.list, is.null)]
log(sprintf("✅ Loaded %d samples", length(ball.list)))

# Step 2: Preprocess individually
for (i in seq_along(ball.list)) {
  ball.list[[i]] <- NormalizeData(ball.list[[i]])
  ball.list[[i]] <- FindVariableFeatures(ball.list[[i]], selection.method = "vst", nfeatures = 2000)
  ball.list[[i]] <- ScaleData(ball.list[[i]])
  ball.list[[i]] <- RunPCA(ball.list[[i]], npcs = 50)
}

# Step 3: Integration with RPCA
log("🔍 Finding RPCA anchors...")
features <- SelectIntegrationFeatures(ball.list, nfeatures = 2500)
anchors <- FindIntegrationAnchors(ball.list, anchor.features = features, reduction = "rpca", dims = 1:40)

log("🧬 Integrating samples...")
integrated <- IntegrateData(anchors, dims = 1:40)
DefaultAssay(integrated) <- "integrated"

# Step 4: Clustering and Harmony
log("⚙️ Clustering integrated object...")
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 50)
integrated$orig.ident <- as.factor(integrated$orig.ident)
integrated <- RunHarmony(integrated, group.by.vars = "orig.ident", max_iter = 20)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindClusters(integrated, resolution = c(0.3, 0.6, 0.9))

# Step 5: Save results
log("💾 Saving Seurat object and metadata...")
saveRDS(integrated, "Core_Results/integrated_data.rds")
write.csv(integrated@meta.data, "Core_Results/integration_metadata.csv")

# Step 6: Generate plots
log("📊 Generating UMAP and QC plots...")
pdf("Core_Results/core_clustering_plots.pdf", width = 14)
DimPlot(integrated, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clustering Overview")
VlnPlot(integrated, features = c("nFeature_RNA", "percent.mt"), pt.size = 0)
dev.off()

log("✅ Final integration pipeline complete! Results are in Core_Results/")
EOF

nohup Rscript BALL_seurat_final.R > seurat_final.log 2>&1 &
echo "🔬 Monitor with: tail -f seurat_final.log"
