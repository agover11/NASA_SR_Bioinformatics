# OSD-404 scRNA-seq Analysis Pipeline
# RRRM-2 Mission: PBMC from LAR flight and ground control mice
# Conditions: FLT vs GC, YNG (12wk) vs OLD (29wk), all female C57BL/6NTac
install.packages("harmony")
install.packages("clustree")
BiocManager::install("SingleR")
BiocManager::install("celldex")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(harmony)

setwd("/Users/lexigover/Desktop/OSD-404/GLDS404_STARsolo")

#sample metadata: 8 flight (FLT) and 8 ground control (GC) mice split by age
sample_info <- data.frame(
  sample    = c("FLT_LAR_YNG_FY9",  "FLT_LAR_YNG_FY11", "FLT_LAR_YNG_FY14", "FLT_LAR_YNG_FY15",
                "FLT_LAR_OLD_FO4",  "FLT_LAR_OLD_FO16", "FLT_LAR_OLD_FO19", "FLT_LAR_OLD_FO20",
                "GC_LAR_YNG_GY4",   "GC_LAR_YNG_GY6",   "GC_LAR_YNG_GY9",   "GC_LAR_YNG_GY10",
                "GC_LAR_OLD_GO14",  "GC_LAR_OLD_GO16",  "GC_LAR_OLD_GO19",  "GC_LAR_OLD_GO20"),
  condition = c(rep("FLT", 8), rep("GC", 8)),
  age       = c(rep("YNG", 4), rep("OLD", 4), rep("YNG", 4), rep("OLD", 4))
)

#load STARsolo GeneFull filtered matrices for each sample
#processed by NASA GeneLab pipeline GL-DPPD-7111, reference GRCm39 Ensembl 107
#using ReadMtx instead of Read10X bc the files are uncompressed tsv
seurat_list <- list()

for (i in seq_len(nrow(sample_info))) {
  cat("Loading", sample_info$sample[i], "\n")
  
  path <- file.path(paste0(sample_info$sample[i], "_GeneFull"),
                    paste0("RRRM2_PBMC_", sample_info$sample[i], "_GeneFull"),
                    "filtered")
  
  counts <- ReadMtx(mtx      = file.path(path, "matrix.mtx"),
                    cells    = file.path(path, "barcodes.tsv"),
                    features = file.path(path, "features.tsv"))
  
  #require gene detection in at least 3 cells and at least 200 genes per cell to filter out low quality cells and unreliable genes
  obj <- CreateSeuratObject(counts       = counts,
                            project      = sample_info$sample[i],
                            min.cells    = 3,
                            min.features = 200)
  
  obj$sample    <- sample_info$sample[i]
  obj$condition <- sample_info$condition[i]
  obj$age       <- sample_info$age[i]
  
  #calculate mitochondrial and ribosomal percentages per cell
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  
  seurat_list[[sample_info$sample[i]]] <- obj
}

#merge all 16 samples into one object
combined <- merge(seurat_list[[1]],
                  y            = seurat_list[2:length(seurat_list)],
                  add.cell.ids = sample_info$sample)

#filter low quality cells based on QC metrics
#percent.mt threshold set to 20 instead of 10 bc mouse PBMCs run higher than human (interesting!)
#46,341 cells retained after filtering
combined <- subset(combined,
                   subset = nFeature_RNA > 300  & nFeature_RNA < 5000 &
                     nCount_RNA   > 500  & nCount_RNA   < 30000 &
                     percent.mt   < 20   &
                     percent.ribo < 50)
print(combined)

#check correlation between counts and features after filtering - 0.86 correlation
FeatureScatter(combined, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "sample")

#normalize, find variable features, scale, and run PCA
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 3000)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

#elbow plot to determine number of PCs to use
#gradual slope with no clear elbow — 30 PCs selected as standard upper bound
ElbowPlot(combined, ndims = 50)

#run UMAP on uncorrected PCA to visualize initial clustering before integration
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")

#initial UMAP revealed samples clustering by mouse rather than cell type - indicates batch effect!
DimPlot(combined, reduction = "umap.unintegrated", group.by = "sample")

#use harmony for batch correction across 16 samples
#corrects technical variation while preserving biological signal
combined <- RunHarmony(combined,
                       group.by.vars  = "sample",
                       reduction      = "pca",
                       reduction.save = "harmony")

#cluster and embed using Harmony reduction
#30 PCs selected based on ElbowPlot
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)

#testing resolution before committing to one
library(clustree)
combined <- FindClusters(combined, resolution = c(0.3, 0.5, 0.8, 1.0))
clustree(combined)

#0.5 selected as most stable resolution from clustree comparison - 26 clusters found
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

#check that samples mix well after integration - they do!
DimPlot(combined, group.by = "sample")
DimPlot(combined, group.by = "seurat_clusters", label = TRUE)

#JoinLayers required before FindAllMarkers in Seurat v5
combined <- JoinLayers(combined)

#find top positive markers for each cluster
markers <- FindAllMarkers(combined,
                          only.pos        = TRUE,
                          min.pct         = 0.25,
                          logfc.threshold = 0.25)

markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% print(n = 130)


#download mouse immune reference
ref <- celldex::ImmGenData()

#get normalized counts
counts <- GetAssayData(combined, assay = "RNA", layer = "data")

#run SingleR
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

#add labels to Seurat object
combined$singler_labels <- pred$labels
DimPlot(combined, group.by = "singler_labels", label = TRUE, repel = TRUE)

table(combined$singler_labels)

