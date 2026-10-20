install.packages("Seurat")
install.packages("hdf5r")
install.packages("tidyverse")

library(Seurat)
library(hdf5r)
library(tidyverse)


setwd("/Users/lexigover/Desktop/OSD918_matrices")

#sample metadata
sample_info <- data.frame(
  sample = c("HC_LAR_19", "HC_LAR_20",
             "FL_LAR_19", "FL_LAR_20"),
  group  = c("Ground", "Ground",
             "Flight", "Flight"),
  age    = c("8mo", "8mo",
             "8mo", "8mo")
)


seurat_list <- list()


for (i in seq_len(nrow(sample_info))) {
  
  #print which sample we're on to can track progress
  cat("Loading", sample_info$sample[i], "\n")
  
  #read the .h5 file for sample i
  mat <- Read10X_h5(paste0(sample_info$sample[i], "_filtered_feature_bc_matrix.h5"))
  
  #turn it into a Seurat object
  #require gene detection in at least 3 cells and at least 200 genes per cell to filter out low quality cells and unreliable genes
  obj <- CreateSeuratObject(counts = mat, 
                            project = sample_info$sample[i], 
                            min.cells = 3, 
                            min.features = 200)
  
  #tag every cell with its sample, group, and age
  obj$sample <- sample_info$sample[i]
  obj$group  <- sample_info$group[i]
  obj$age    <- sample_info$age[i]
  
  #calculate mitochondrial and ribosomal percentages per cell
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  
  #save this sample's object into the list
  seurat_list[[sample_info$sample[i]]] <- obj
}

#merge all 8 into one object
combined <- merge(seurat_list[[1]], 
                  y = seurat_list[-1], 
                  add.cell.ids = sample_info$sample)

#quick check
print(combined)


VlnPlot(combined, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "sample",
        pt.size = 0.1,
        ncol = 2)

#Lenient filtering 
combined <- subset(combined,
                   subset = nFeature_RNA > 300 &
                     nFeature_RNA < 5000 &
                     nCount_RNA > 500 &
                     nCount_RNA < 30000 &
                     percent.mt < 10 &
                     percent.ribo < 40)

print(combined)

FeatureScatter(combined, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "sample")

table(combined$sample)

#stricter filtering
combined <- subset(combined,
                   subset = nFeature_RNA > 500 &
                     nFeature_RNA < 5000 &
                     nCount_RNA > 500 &
                     nCount_RNA < 30000 &
                     percent.mt < 10 &
                     percent.ribo < 40)

print(combined)

FeatureScatter(combined, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "sample")

table(combined$sample)

