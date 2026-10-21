#OSD-404 scRNA-seq Analysis Pipeline
#RRRM-2 Mission: PBMC from LAR flight and ground control mice
#Conditions: FLT vs GC, OLD (29wk), all female C57BL/6NTac

# install.packages("harmony")
# install.packages("clustree")
# install.packages("BiocManager")
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("scDblFinder")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(harmony)
library(clustree)
library(readxl)
library(scDblFinder)
library(SingleCellExperiment)
library(ggrepel)
library(DESeq2)
library(edgeR)
library(limma)


set.seed(42)
setwd("/Users/lexigover/Desktop/OSD-404/GLDS404_STARsolo")

#sample metadata: 4 flight (FLT) and 4 ground control (GC) mice, OLD age only
sample_info <- data.frame(
  sample    = c("FLT_LAR_OLD_FO4",  "FLT_LAR_OLD_FO16", "FLT_LAR_OLD_FO19", "FLT_LAR_OLD_FO20",
                "GC_LAR_OLD_GO14",  "GC_LAR_OLD_GO16",  "GC_LAR_OLD_GO19",  "GC_LAR_OLD_GO20"),
  condition = c(rep("FLT", 4), rep("GC", 4)),
  age       = rep("OLD", 8)
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


#merge all samples into one object
combined <- merge(seurat_list[[1]],
                  y            = seurat_list[2:length(seurat_list)],
                  add.cell.ids = sample_info$sample)


#filter low quality cells based on QC metrics - should have filtered for each sample individually 
#percent.mt threshold set to 20 instead of 10 bc mouse PBMCs run higher than human (interesting!)
#25,352 cells retained after filtering
combined <- subset(combined,
                   subset = nFeature_RNA > 300  & nFeature_RNA < 5000 &
                     nCount_RNA   > 500  & nCount_RNA   < 30000 &
                     percent.mt   < 20   &
                     percent.ribo < 50)
print(combined)

#check correlation between counts and features after filtering - 0.82 correlation!
FeatureScatter(combined,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               group.by = "sample")

#on seurat_list, before merging/filtering, if you still have it, or rerun QC metrics per sample
combined@meta.data %>%
  group_by(sample) %>%
  summarise(
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt),
    n = n()
  )

#Remove hemoglobin
grep("^Hb[ab]-", rownames(combined), value = TRUE)

is_hb <- grepl("^Hb[ab]-", rownames(combined))
cat("Hb genes removed:", sum(is_hb), "-", rownames(combined)[is_hb], "\n")
combined <- combined[!is_hb, ]

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
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30,
                    reduction.name = "umap.unintegrated")

#initial UMAP revealed samples clustering by mouse rather than cell type - indicates batch effect!
DimPlot(combined, reduction = "umap.unintegrated", group.by = "sample")


#use harmony for integration/batch correction across 8 samples
#corrects technical variation while preserving biological signal
combined <- RunHarmony(combined,
                       group.by.vars  = "sample",
                       reduction      = "pca",
                       reduction.save = "harmony")

#cluster and embed using Harmony reduction
#30 PCs selected based on ElbowPlot
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)

#testing multiple resolution before committing to one
combined <- FindClusters(combined, resolution = c(0.3, 0.5, 0.8, 1.0))
clustree(combined)

#0.5 selected as most stable resolution from clustree comparison - 23 clusters found
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

#check that samples mix well after integration - they do!
DimPlot(combined, group.by = "sample")
DimPlot(combined, group.by = "seurat_clusters", label = TRUE)



#MARKER GENE DETECTION
#JoinLayers required before FindAllMarkers in Seurat v5
combined <- JoinLayers(combined)

#find top positive markers for each cluster
markers <- FindAllMarkers(combined,
                          only.pos        = TRUE,
                          min.pct         = 0.25,
                          logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% print(n = 250)


combined@meta.data %>%
  filter(seurat_clusters == "22") %>%
  summarise(mean_nCount = mean(nCount_RNA), mean_nFeature = mean(nFeature_RNA), mean_pct_mt = mean(percent.mt), n = n())

#check how many cells in each cluster
table(combined$seurat_clusters)

#Remove cluster 22 (low-quality, ribosomal-only, 58 cells)
combined <- subset(combined, seurat_clusters != "22")
table(combined$seurat_clusters)  # should show 0-21, 22 gone

combined$seurat_clusters <- droplevels(combined$seurat_clusters)
table(combined$seurat_clusters)


#SINGLER AUTOMATED ANNOTATION
#download mouse immune reference
ref    <- celldex::ImmGenData()
#get normalized counts
counts <- GetAssayData(combined, assay = "RNA", layer = "data")

#run SingleR using ImmGen mouse immune reference
#assigns a cell type label to each individual cell based on reference gene expression profiles
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

#add SingleR labels to Seurat object as metadata
combined$singler_labels <- pred$labels


#view SingleR labels on Harmony-corrected UMAP
DimPlot(combined, reduction = "umap", group.by = "singler_labels",
        label = TRUE, repel = TRUE)

#check how many cells were assigned to each label
table(combined$singler_labels)


#view numbered clusters on same UMAP for comparison with SingleR labels
DimPlot(combined, reduction = "umap", group.by = "seurat_clusters",
        label = TRUE, repel = TRUE)

#remove label categories with fewer than 20 cells to reduce visual noise
small_labels <- names(table(combined$singler_labels))[table(combined$singler_labels) < 20]
combined$singler_labels_clean <- ifelse(combined$singler_labels %in% small_labels,
                                        NA, combined$singler_labels)

#view cleaned SingleR labels
DimPlot(combined, reduction = "umap", group.by = "singler_labels_clean",
        label = TRUE, repel = TRUE)

#compare SingleR labels against cluster numbers
table(combined$singler_labels, combined$seurat_clusters)


#FEATUREPLOTS FOR AMBIGUOUS CLUSTERS
#cluster 1 — monocyte vs proliferating population
FeaturePlot(combined,
            features = c("F13a1", "Csf1r", "Mki67", "Top2a"),
            reduction = "umap", ncol = 2)

#cluster 2 — progenitor identity (HSC vs GMP vs other)
FeaturePlot(combined,
            features = c("Kit", "Cd34", "Elane", "Mpo"),
            reduction = "umap", ncol = 2)

#cluster 8 — mature vs immature neutrophil signature
FeaturePlot(combined,
            features = c("S100a8", "Camp", "Elane", "Ms4a3"),
            reduction = "umap", ncol = 2)

#check cluster 5 to see about subclustering
cluster5 <- subset(combined, seurat_clusters == "5")
cluster5 <- FindNeighbors(cluster5, reduction = "harmony", dims = 1:30)
cluster5 <- FindClusters(cluster5, resolution = 0.3)
cluster5 <- RunUMAP(cluster5, reduction = "harmony", dims = 1:30, reduction.name = "umap")
DimPlot(cluster5, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
FeaturePlot(cluster5, features = c("Cd3e", "Ncr1", "Trdc", "Klrb1c"), reduction = "umap", ncol = 2)

cluster5_markers <- FindAllMarkers(JoinLayers(cluster5), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster5_markers %>% filter(cluster == "3") %>% arrange(desc(avg_log2FC)) %>% head(15)

#cell type annotation based on marker genes, SingleR, and FeaturePlot validation
cell_types <- c(
  "0"  = "B cells",
  "1"  = "Monocytes",
  "2"  = "HSCs/Progenitors",
  "3"  = "Monocytes (inflammatory)",
  "4"  = "B cells (naive)",
  "5"  = "T cells",
  "6"  = "Monocytes (patrolling)",
  "7"  = "B cells (atypical/memory)",
  "8"  = "Neutrophils (immature)",
  "9"  = "T cells (CD8+)",
  "10" = "Dendritic cells",
  "11" = "Plasmacytoid DC",
  "12" = "B cells (progenitor)",
  "13" = "Innate B cells",
  "14" = "Plasma cells",
  "15" = "Neutrophils",
  "16" = "Macrophages",
  "17" = "NK cells",
  "18" = "Platelets",
  "19" = "Neutrophils (inflammatory)",
  "20" = "Proerythroblasts",
  "21" = "Monocyte-derived Macrophages"
)

#apply cell type labels to each cell based on cluster assignment
combined$cell_type <- as.character(cell_types[as.character(combined$seurat_clusters)])

#confirm final cell type count
length(unique(combined$cell_type))

#ANNOTATED UMAP
#view final annotated UMAP
DimPlot(combined, reduction = "umap", group.by = "cell_type",
        label = TRUE, repel = TRUE)

#UMAP formatted for poster
p <- DimPlot(combined, reduction = "umap", group.by = "cell_type",
             label = TRUE, label.size = 5, repel = TRUE,
             pt.size = 0.5) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.35, "cm"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
print(p)

table(combined$seurat_clusters, combined$cell_type)

ncol(combined)
table(combined$condition)
grep("^Hb[ab]-", rownames(combined), value = TRUE)
length(unique(combined$cell_type))


#CANONICAL MARKER DOTPLOT
#canonical mouse immune marker panel, grouped by lineage
#check against data first so DotPlot doesn't choke on missing genes
canonical_markers <- c(
  # basophil / mast (relevant to the cluster 3 ambiguity)
  "Mcpt8", "Ms4a2", "Cpa3", "Prss34",
  # platelet / megakaryocyte
  "Itga2b", "Pf4", "Ppbp", "Gp9",
  # monocyte / monoblast
  "Csf1r", "Cd14", "Ly6c2", "Ccr2", "F13a1", "Nr4a1", "Cx3cr1", "Fcgr1",
  # neutrophil
  "S100a8", "S100a9", "Retnlg", "Camp", "Ltf", "Cxcr2",
  # dendritic cell
  "Flt3", "Zbtb46", "Itgax", "Xcr1", "Clec9a",
  # plasmacytoid DC
  "Siglech", "Bst2", "Irf8",
  # macrophage (relevant to the PBMC-macrophage question)
  "Adgre1", "Mrc1", "Mertk", "Spic", "C1qa",
  # T cell
  "Cd3e", "Cd3d", "Cd4", "Cd8a", "Il7r", "Ccr7",
  # cytotoxic / NK
  "Gzmb", "Nkg7", "Ncr1", "Klrb1c",
  # B cell (naive / activated / progenitor / plasma)
  "Ms4a1", "Cd19", "Cd79a", "Ighd", "Cd83", "Vpreb3", "Dntt", "Rag1", "Sdc1", "Jchain", "Prdm1"
)

missing_markers <- setdiff(canonical_markers, rownames(combined))
if (length(missing_markers) > 0) {
  cat("Not found in dataset, dropping:", paste(missing_markers, collapse = ", "), "\n")
}
canonical_markers <- intersect(canonical_markers, rownames(combined))

#raw cluster numbers - unbiased check before labels get involved
DotPlot(combined, features = canonical_markers, group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7))

#assigned cell_type labels, a "clean" annotation should show each row lighting up narrowly on its expected markers; a messy/bleeding
#pattern across unrelated marker groups flags a label worth rechecking
DotPlot(combined, features = canonical_markers, group.by = "cell_type") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7))

#POSTER DOTPLOT
cell_type_order <- c(
  "Plasma cells",
  "B cells",
  "B cells (naive)",
  "B cells (progenitor)",
  "B cells (atypical/memory)",
  "Innate B cells",
  "NK cells",
  "T cells (CD8+)",
  "T cells",
  "Macrophages",
  "Monocyte-derived Macrophages",
  "Plasmacytoid DC",
  "Dendritic cells",
  "Neutrophils (inflammatory)",
  "Neutrophils (immature)",
  "Neutrophils",
  "Monocytes (patrolling)",
  "Monocytes (inflammatory)",
  "Monocytes",
  "HSCs/Progenitors",
  "Proerythroblasts",
  "Platelets"
)

#one representative marker per cell type (roughly), ordered for diagonal expression pattern
poster_markers <- c(
  "Itga2b",   # platelets
  "Car2",     # proerythroblasts
  "Cd34",     # HSCs/Progenitors
  "F13a1",    # monocytes
  "Ccr2",     # monocytes (inflammatory)
  "Nr4a1",    # monocytes (patrolling)
  "Ly6g",     # neutrophils
  "Elane",    # neutrophils (immature)
  "Acod1",    # neutrophils (inflammatory)
  "Flt3",   # dendritic cells
  "Siglech",  # plasmacytoid DC
  "Mertk",    # monocyte-derived macrophages
  "Gata6",    # macrophages
  "Lef1",     # T cells
  "Cd8b1",    # T cells (CD8+)
  "Nkg7",     # NK cells
  "Sox5",    # B cells
  "Il5ra",    # innate B cells
  "Zbtb32",   # B cells (atypical/memory)
  "Vpreb3",   # B cells (progenitor)
  "Ighd",     # B cells (naive)
  "Jchain"    # plasma cells
)

p_dotplot <- DotPlot(combined,
                     features = poster_markers,
                     group.by = "cell_type") +
  scale_y_discrete(limits = cell_type_order) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 14))

p_dotplot

FetchData(combined, vars = c("Zbtb46", "Flt3", "cell_type")) %>%
  filter(cell_type %in% c("Dendritic cells", "Plasmacytoid DC", "Monocytes",
                          "Monocytes (inflammatory)", "Macrophages")) %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), list(mean = mean, pct = ~mean(. > 0))))

FetchData(combined, vars = c("Retnlg", "Acod1", "cell_type")) %>%
  filter(cell_type %in% c("Neutrophils", "Neutrophils (inflammatory)")) %>%
  group_by(cell_type) %>%
  summarise(mean_Retnlg = mean(Retnlg), pct_Retnlg = mean(Retnlg > 0),
            mean_Acod1 = mean(Acod1), pct_Acod1 = mean(Acod1 > 0))

FetchData(combined, vars = c("Mmp8", "Mmp9", "Ly6g", "Ltf", "cell_type")) %>%
  filter(cell_type %in% c("Neutrophils", "Neutrophils (immature)", "Neutrophils (inflammatory)")) %>%
  group_by(cell_type) %>%
  summarise(across(everything(), list(mean = mean, pct = ~mean(. > 0))))

FetchData(combined, vars = c("Sox5", "Ighd", "Vpreb3", "Zbtb32", "Il5ra", "cell_type")) %>%
  filter(cell_type %in% c("B cells", "B cells (naive)", "B cells (progenitor)",
                          "B cells (atypical/memory)", "Innate B cells")) %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), list(mean = mean, pct = ~mean(. > 0))))

FetchData(combined, vars = c("Ebf1", "Bank1", "Ralgps2", "Fchsd2", "cell_type")) %>%
  filter(cell_type %in% c("B cells", "B cells (naive)", "B cells (progenitor)",
                          "B cells (atypical/memory)", "Innate B cells")) %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), list(mean = mean, pct = ~mean(. > 0))))

FetchData(combined, vars = c("Cd4", "cell_type")) %>%
  filter(cell_type %in% c("T cells", "T cells (CD8+)")) %>%
  group_by(cell_type) %>%
  summarise(mean = mean(Cd4), pct = mean(Cd4 > 0))

markers %>% filter(cluster == "0") %>% arrange(desc(avg_log2FC)) %>% head(15)

FetchData(combined, vars = c("Pax5", "Btla", "cell_type")) %>%
  filter(cell_type %in% c("B cells", "B cells (naive)", "B cells (progenitor)",
                          "B cells (atypical/memory)", "Innate B cells")) %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), list(mean = mean, pct = ~mean(. > 0))))

#QC AND DOUBLET CHECKS
#basic QC metrics per cell type -- doublets often run higher on both
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        group.by = "cell_type", pt.size = 0)

#check whether "macrophage" cells also carry strong platelet signal
#(Hba-a1 removed - no longer in object after Hb gene filtering)
mac_cells <- subset(combined, cell_type == "Monocyte-derived Macrophages")
VlnPlot(mac_cells, features = c("Adgre1", "Pf4", "Itga2b", "Ppbp"), pt.size = 0)

#doublet check across final annotation
combined <- JoinLayers(combined)
sce <- as.SingleCellExperiment(combined)
sce <- scDblFinder(sce, samples = "sample")
combined$scDblFinder_class <- sce$scDblFinder.class

table(combined$cell_type, combined$scDblFinder_class)

#CELL COMPOSITION ANALYSIS
#check cluster sizes
table(combined$seurat_clusters)
#check overall cell counts per condition
table(combined$condition)
#check per-sample cell counts to assess within-group balance
table(combined$sample, combined$condition)
table(combined$cell_type, combined$condition)
#what fraction of total cells is this?
round(prop.table(table(combined$cell_type)) * 100, 2)

#normalize by total cells per condition to get proportions
prop_table <- combined@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = round(n / sum(n) * 100, 2)) %>%
  arrange(cell_type)
print(prop_table, n = 50)

#stacked bar plot of cell type proportions per sample
combined@meta.data %>%
  group_by(sample, condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ggplot(aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(title = "Cell type composition per sample",
       x = NULL, y = "Proportion (%)", fill = "Cell type") +
  scale_x_discrete(labels = c(
    "FLT_LAR_OLD_FO4"  = "FO4",  "FLT_LAR_OLD_FO16" = "FO16",
    "FLT_LAR_OLD_FO19" = "FO19", "FLT_LAR_OLD_FO20" = "FO20",
    "GC_LAR_OLD_GO14"  = "GO14", "GC_LAR_OLD_GO16"  = "GO16",
    "GC_LAR_OLD_GO19"  = "GO19", "GC_LAR_OLD_GO20"  = "GO20"
  ))

#flight vs control (no individual samples)
combined@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ggplot(aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(title = "Cell type composition: Flight vs Ground Control",
       x = NULL, y = "Proportion (%)", fill = "Cell type")

#poster cell proportion plot
p_composition <- combined@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ggplot(aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold")) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 5))) +
  labs(title = "Cell type composition: Flight vs Ground Control",
       x = NULL, y = "Proportion (%)", fill = "Cell type")
print(p_composition)

#cell type proportions
combined@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ggplot(aes(x = cell_type, y = proportion, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("FLT" = "#2166AC", "GC" = "#D6604D")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Cell type proportions: Flight vs Ground Control",
       x = NULL, y = "Proportion (%)", fill = "Condition")

#cell type counts
combined@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = cell_type, y = n, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("FLT" = "#2166AC", "GC" = "#D6604D")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Cell type counts: Flight vs Ground Control",
       x = NULL, y = "Cell count", fill = "Condition")

#visualize condition distribution across UMAP
DimPlot(combined,
        reduction = "umap",
        group.by = "condition",
        cols = c("FLT" = "#2166AC", "GC" = "#D6604D"),
        pt.size = 0.3) +
  ggtitle("Flight vs Ground Control")

#split by condition to compare distributions side by side
DimPlot(combined,
        reduction = "umap",
        group.by = "condition",
        cols = c("FLT" = "#2166AC", "GC" = "#D6604D"),
        split.by = "condition",
        pt.size = 0.3)

#MITOCARTA MODULE SCORING
mitocarta <- read_excel("/Users/lexigover/Downloads/Mouse.MitoCarta3.0.xls",
                        sheet = "A Mouse MitoCarta3.0")
mito_genes <- mitocarta$Symbol
mito_genes <- mito_genes[mito_genes %in% rownames(combined)]
cat("MitoCarta genes found in dataset:", length(mito_genes), "\n")

combined <- AddModuleScore(
  combined,
  features = list(mito_genes),
  name = "MitoCarta_score"
)
VlnPlot(combined, features = "MitoCarta_score1", group.by = "cell_type", pt.size = 0)
FeaturePlot(combined, features = "MitoCarta_score1", reduction = "umap")
VlnPlot(combined, features = "MitoCarta_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)

#test MitoCarta score difference between FLT and GC within each cell type (Wilcoxon)
score_test <- combined@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    FLT_mean = mean(MitoCarta_score1[condition == "FLT"]),
    GC_mean  = mean(MitoCarta_score1[condition == "GC"]),
    p_value  = wilcox.test(MitoCarta_score1[condition == "FLT"],
                           MitoCarta_score1[condition == "GC"])$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
print(score_test, n = 24)

#MITODAMP PATHWAY GENE SETS
read_gmt <- function(file) {
  lines <- readLines(file)
  genes <- strsplit(lines, "\t")[[1]]
  genes <- genes[-(1:2)]
  return(genes)
}
mitophagy_genes  <- read_gmt("~/Desktop/GOBP_MITOPHAGY.v2026.1.Mm.gmt")
cgas_sting_genes <- read_gmt("~/Desktop/GOBP_CGAS_STING_SIGNALING_PATHWAY.v2026.1.Mm.gmt")
nlrp3_genes      <- read_gmt("~/Desktop/GOBP_NLRP3_INFLAMMASOME_COMPLEX_ASSEMBLY.v2026.1.Mm.gmt")
aim2_genes       <- read_gmt("~/Desktop/GOBP_AIM2_INFLAMMASOME_COMPLEX_ASSEMBLY.v2026.1.Mm.gmt")
ifn_genes        <- read_gmt("~/Desktop/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2026.1.Mm.gmt")

cat("Mitophagy genes:", length(mitophagy_genes), "\n")
cat("cGAS-STING genes:", length(cgas_sting_genes), "\n")
cat("NLRP3 genes:", length(nlrp3_genes), "\n")
cat("AIM2 genes:", length(aim2_genes), "\n")
cat("IFN genes:", length(ifn_genes), "\n")

mitophagy_genes  <- mitophagy_genes[mitophagy_genes %in% rownames(combined)]
cgas_sting_genes <- cgas_sting_genes[cgas_sting_genes %in% rownames(combined)]
nlrp3_genes      <- nlrp3_genes[nlrp3_genes %in% rownames(combined)]
aim2_genes       <- aim2_genes[aim2_genes %in% rownames(combined)]
ifn_genes        <- ifn_genes[ifn_genes %in% rownames(combined)]

combined <- AddModuleScore(combined, features = list(mitophagy_genes),  name = "Mitophagy_score")
combined <- AddModuleScore(combined, features = list(cgas_sting_genes), name = "cGAS_STING_score")
combined <- AddModuleScore(combined, features = list(nlrp3_genes),      name = "NLRP3_score")
combined <- AddModuleScore(combined, features = list(aim2_genes),       name = "AIM2_score")
combined <- AddModuleScore(combined, features = list(ifn_genes),        name = "IFN_score")

VlnPlot(combined, features = "Mitophagy_score1", group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "cGAS_STING_score1", group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "NLRP3_score1", group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "AIM2_score1", group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "IFN_score1", group.by = "cell_type", split.by = "condition", pt.size = 0)

#Wilcoxon test: each pathway score FLT vs GC within each cell type
score_names <- c("Mitophagy_score1", "cGAS_STING_score1",
                 "NLRP3_score1", "AIM2_score1", "IFN_score1")
for (score in score_names) {
  cat("\n===", score, "===\n")
  result <- combined@meta.data %>%
    group_by(cell_type) %>%
    summarise(
      FLT_mean = mean(.data[[score]][condition == "FLT"]),
      GC_mean  = mean(.data[[score]][condition == "GC"]),
      p_value  = wilcox.test(.data[[score]][condition == "FLT"],
                             .data[[score]][condition == "GC"])$p.value
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
  print(result, n = 22)
}

#subset Macrophages for pathway and gene expression plots
mac_ifn <- subset(combined, cell_type == "Macrophages")
VlnPlot(mac_ifn, features = "IFN_score1", group.by = "condition",
        cols = c("FLT" = "#2166AC", "GC" = "#D6604D"), pt.size = 0.5) +
  labs(title = "Interferon Alpha Response Score", x = NULL, y = "IFN Module Score") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16), plot.title = element_text(size = 18, face = "bold"))

VlnPlot(mac_ifn, features = "Snca", group.by = "condition",
        cols = c("FLT" = "#2166AC", "GC" = "#D6604D"), pt.size = 0.5) +
  labs(title = "Snca Expression", x = NULL, y = "Expression Level") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16), plot.title = element_text(size = 18, face = "bold"))

#CELL-LEVEL DEG ANALYSIS
#overall DEG analysis: all FLT vs all GC cells regardless of cell type
Idents(combined) <- "condition"
all_degs <- FindMarkers(combined, ident.1 = "FLT", ident.2 = "GC",
                        min.pct = 0.25, logfc.threshold = 0.25)
all_degs$gene <- rownames(all_degs)
all_degs <- all_degs %>%
  mutate(direction = case_when(
    avg_log2FC > 0.25  & p_val_adj < 0.05 ~ "Up in FLT",
    avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Up in GC",
    TRUE ~ "NS"
  ))
top10_up   <- all_degs %>% filter(direction == "Up in FLT") %>% arrange(desc(avg_log2FC)) %>% head(10)
top10_down <- all_degs %>% filter(direction == "Up in GC") %>% arrange(avg_log2FC) %>% head(10)
top10 <- bind_rows(top10_up, top10_down)
ggplot(all_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(data = top10, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Up in FLT" = "#2166AC", "Up in GC" = "#D6604D", "NS" = "grey70")) +
  theme_classic() +
  labs(title = "Overall DEGs: Flight vs Ground Control",
       x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

#DEG analysis within each cell type (updated list to match current 21 cell types)
Idents(combined) <- "cell_type"
cell_types_to_test <- c(
  "Monocytes", "Monocytes (inflammatory)", "Monocytes (patrolling)",
  "Neutrophils", "Neutrophils (immature)", "Neutrophils (inflammatory)",
  "NK cells", "T cells", "T cells (CD8+)",
  "B cells", "B cells (naive)", "B cells (progenitor)",
  "B cells (atypical/memory)", "Innate B cells",
  "Macrophages", "Monocyte-derived Macrophages",
  "Dendritic cells", "Plasmacytoid DC"
)
deg_results <- list()
for (ct in cell_types_to_test) {
  cat("Running DEG for", ct, "\n")
  tryCatch({
    degs <- FindMarkers(combined, ident.1 = "FLT", ident.2 = "GC",
                        group.by = "condition", subset.ident = ct,
                        min.pct = 0.25, logfc.threshold = 0.25)
    degs$gene      <- rownames(degs)
    degs$cell_type <- ct
    degs <- degs %>%
      mutate(direction = case_when(
        avg_log2FC > 0.25  & p_val_adj < 0.05 ~ "Up in FLT",
        avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Up in GC",
        TRUE ~ "NS"
      ))
    deg_results[[ct]] <- degs

    top10_up   <- degs %>% filter(direction == "Up in FLT") %>% arrange(desc(avg_log2FC)) %>% head(10)
    top10_down <- degs %>% filter(direction == "Up in GC") %>% arrange(avg_log2FC) %>% head(10)
    top10_ct <- bind_rows(top10_up, top10_down)

    p <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_text_repel(data = top10_ct, aes(label = gene), size = 3, max.overlaps = 20) +
      scale_color_manual(values = c("Up in FLT" = "#2166AC", "Up in GC" = "#D6604D", "NS" = "grey70")) +
      theme_classic() +
      labs(title = paste("DEGs:", ct, "— Flight vs Ground Control"),
           x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")
    print(p)

  }, error = function(e) {
    cat("Skipping", ct, "- insufficient cells\n")
  })
}
all_cell_type_degs <- bind_rows(deg_results)
deg_summary <- bind_rows(deg_results) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cell_type, direction) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  arrange(cell_type)
print(deg_summary, n = 50)
deg_summary %>%
  ggplot(aes(x = reorder(cell_type, n_degs), y = n_degs, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Up in FLT" = "#2166AC", "Up in GC" = "#D6604D")) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title   = element_text(size = 18, face = "bold"),
        legend.text  = element_text(size = 16),
        legend.key.size = unit(0.5, "cm")) +
  labs(title = "Significant DEGs per cell type: Flight vs Ground Control",
       x = NULL, y = "Number of significant DEGs", fill = NULL)


#PSEUDOBULK DESeq2
pseudo_bulk <- AggregateExpression(
  object        = combined,
  group.by      = c("cell_type", "orig.ident", "condition"),
  assays        = "RNA",
  return.seurat = FALSE
)
cts <- pseudo_bulk$RNA
deg_results_deseq2 <- list()
for (ct in cell_types_to_test) {
  cat("\nDESeq2:", ct, "\n")
  ct_escaped <- gsub("([\\+\\(\\)])", "\\\\\\1", ct)
  ct_cols    <- colnames(cts)[grep(paste0("^", ct_escaped, "[-_]"), colnames(cts))]

  if (length(ct_cols) < 6) {
    cat("skipping", ct, "- fewer than 6 pseudobulk samples\n"); next
  }

  cluster_counts <- cts[, ct_cols]
  conditions     <- sapply(strsplit(colnames(cluster_counts), "[-_]"), function(x) tail(x, 1))
  meta           <- data.frame(row.names = colnames(cluster_counts),
                               condition = factor(conditions, levels = c("GC", "FLT")))

  if (length(unique(meta$condition)) < 2) {
    cat("skipping", ct, "- only one condition present\n"); next
  }

  tryCatch({
    dds  <- DESeqDataSetFromMatrix(countData = cluster_counts, colData = meta, design = ~ condition)
    #gene filter tightened for small n (4 vs 4): counts >=10 in at least 4 samples
    #(aligned with limma-voom filter below so both methods test the same gene set)
    keep <- rowSums(counts(dds) >= 10) >= 4
    dds  <- dds[keep, ]
    dds  <- DESeq(dds, quiet = TRUE)
    res  <- results(dds, contrast = c("condition", "FLT", "GC"))

    degs <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      rename(p_val_adj = padj, avg_log2FC = log2FoldChange, p_val = pvalue) %>%
      mutate(cell_type = ct) %>%
      filter(!is.na(p_val)) %>%
      mutate(direction = case_when(
        avg_log2FC > 0.25  & p_val_adj < 0.05 ~ "Up in FLT",
        avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Up in GC",
        TRUE ~ "NS"
      ))
    deg_results_deseq2[[ct]] <- degs

    top10_up   <- degs %>% filter(direction == "Up in FLT") %>% arrange(p_val_adj) %>% head(10)
    top10_down <- degs %>% filter(direction == "Up in GC")  %>% arrange(p_val_adj) %>% head(10)
    top10_ct   <- bind_rows(top10_up, top10_down)

    p <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_text_repel(data = top10_ct, aes(label = gene), size = 3, max.overlaps = 20) +
      scale_color_manual(values = c("Up in FLT" = "#2166AC", "Up in GC" = "#D6604D", "NS" = "grey70")) +
      theme_classic() +
      labs(title = paste("Pseudobulk DESeq2:", ct, "- Flight vs Ground Control"),
           x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")
    print(p)

  }, error = function(e) { cat("error for", ct, ":", e$message, "\n") })
}
all_deseq2_degs <- bind_rows(deg_results_deseq2)
deg_summary_deseq2 <- all_deseq2_degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cell_type, direction) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  arrange(cell_type)
print(deg_summary_deseq2, n = 50)

#p-value distribution sanity check
hist(deg_results_deseq2[["Monocytes"]]$p_val, breaks = 50, main = "Monocytes raw p-values")
hist(deg_results_deseq2[["T cells"]]$p_val, breaks = 50, main = "T cells raw p-values")

#LIMMA-VOOM
deg_results_limma <- list()
for (ct in cell_types_to_test) {
  cat("\nlimma:", ct, "\n")
  ct_escaped <- gsub("([\\+\\(\\)])", "\\\\\\1", ct)
  ct_cols    <- colnames(cts)[grep(paste0("^", ct_escaped, "[-_]"), colnames(cts))]

  if (length(ct_cols) < 6) {
    cat("skipping", ct, "- fewer than 6 pseudobulk samples\n"); next
  }

  cluster_counts <- cts[, ct_cols]
  conditions     <- factor(
    sapply(strsplit(colnames(cluster_counts), "[-_]"), function(x) tail(x, 1)),
    levels = c("GC", "FLT")
  )

  if (length(unique(conditions)) < 2) {
    cat("skipping", ct, "- only one condition present\n"); next
  }

  tryCatch({
    #gene filter aligned with DESeq2 above: counts >=10 in at least 4 samples
    keep           <- rowSums(cluster_counts >= 10) >= 4
    cluster_counts <- cluster_counts[keep, ]
    design <- model.matrix(~ conditions)
    v      <- voom(cluster_counts, design, plot = FALSE)
    fit    <- lmFit(v, design)
    fit    <- eBayes(fit)

    res <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
      rownames_to_column("gene") %>%
      mutate(cell_type = ct,
             direction = case_when(
               logFC > 0.25  & adj.P.Val < 0.05 ~ "Up in FLT",
               logFC < -0.25 & adj.P.Val < 0.05 ~ "Up in GC",
               TRUE ~ "NS"
             ))
    deg_results_limma[[ct]] <- res

    top10_up   <- res %>% filter(direction == "Up in FLT") %>% arrange(adj.P.Val) %>% head(10)
    top10_down <- res %>% filter(direction == "Up in GC")  %>% arrange(adj.P.Val) %>% head(10)
    top10_ct   <- bind_rows(top10_up, top10_down)

    p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_text_repel(data = top10_ct, aes(label = gene), size = 3, max.overlaps = 20) +
      scale_color_manual(values = c("Up in FLT" = "#2166AC", "Up in GC" = "#D6604D", "NS" = "grey70")) +
      theme_classic() +
      labs(title = paste("Pseudobulk limma:", ct, "- Flight vs Ground Control"),
           x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")
    print(p)

  }, error = function(e) { cat("error for", ct, ":", e$message, "\n") })
}
all_limma_degs <- bind_rows(deg_results_limma)
deg_summary_limma <- all_limma_degs %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(cell_type, direction) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  arrange(cell_type)
print(deg_summary_limma, n = 50)

#compare DESeq2 and limma hits
deseq2_sig <- all_deseq2_degs %>%
  filter(p_val_adj < 0.05) %>%
  dplyr::select(gene, cell_type, avg_log2FC, p_val_adj, direction)
limma_sig <- all_limma_degs %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::select(gene, cell_type, logFC, adj.P.Val, direction)
cat("\nDESeq2 significant genes:\n"); print(deseq2_sig)
cat("\nlimma significant genes:\n");  print(limma_sig)


#check for previously-flagged B cell composition-driven genes -
#confirms whether they're real within-celltype effects or composition artifacts
#(see earlier cell-level "Overall DEGs" discussion)
deseq2_sig %>% filter(gene %in% c("Ebf1", "Pax5", "Bank1", "Aff3", "Ikzf3", "Bach2", "Foxp1"))
limma_sig %>% filter(gene %in% c("Ebf1", "Pax5", "Bank1", "Aff3", "Ikzf3", "Bach2", "Foxp1"))

#check hemoglobin-expressing cell impact
#reload raw matrices (same as original loading step)
hb_check_list <- list()
for (i in seq_len(nrow(sample_info))) {
  path <- file.path(paste0(sample_info$sample[i], "_GeneFull"),
                    paste0("RRRM2_PBMC_", sample_info$sample[i], "_GeneFull"),
                    "filtered")
  counts <- ReadMtx(mtx = file.path(path, "matrix.mtx"),
                    cells = file.path(path, "barcodes.tsv"),
                    features = file.path(path, "features.tsv"))
  obj <- CreateSeuratObject(counts = counts, project = sample_info$sample[i],
                            min.cells = 3, min.features = 200)
  obj$sample    <- sample_info$sample[i]
  obj$condition <- sample_info$condition[i]
  hb_check_list[[sample_info$sample[i]]] <- obj
}
hb_check <- merge(hb_check_list[[1]], y = hb_check_list[2:length(hb_check_list)],
                  add.cell.ids = sample_info$sample)
#apply the SAME cell-level QC threshold used originally, so this reflects
#the same population you'd actually be filtering within
hb_check[["percent.mt"]] <- PercentageFeatureSet(hb_check, pattern = "^mt-")
hb_check <- subset(hb_check,
                   subset = nFeature_RNA > 300 & nFeature_RNA < 5000 &
                     nCount_RNA > 500 & nCount_RNA < 30000 &
                     percent.mt < 20)

#compute percent.hb
is_hb <- grepl("^Hb[ab]-", rownames(hb_check))
hb_check[["percent.hb"]] <- PercentageFeatureSet(hb_check, features = rownames(hb_check)[is_hb])
#look at the distribution before picking any threshold
summary(hb_check$percent.hb)
VlnPlot(hb_check, features = "percent.hb", group.by = "sample", pt.size = 0)
#how many cells would be affected at a few candidate thresholds
for (thresh in c(1, 2, 5, 10)) {
  n_affected <- sum(hb_check$percent.hb > thresh)
  pct_affected <- round(100 * n_affected / ncol(hb_check), 2)
  cat("Threshold >", thresh, "%: ", n_affected, "cells (", pct_affected, "% of total)\n")
}
#per-sample breakdown at a reasonable candidate threshold
hb_check@meta.data %>%
  group_by(sample) %>%
  summarise(n_total = n(),
            n_above_2pct = sum(percent.hb > 2),
            pct_above_2pct = round(100 * n_above_2pct / n_total, 2))

deg_results[["Macrophages"]] %>% filter(gene == "Snca")
all_cell_type_degs %>% filter(gene == "Snca") %>% arrange(p_val_adj)

#CONTAMINATION VALIDATION: does Snca / IFN score hold up after Hb correction?
#if hb_check still exists from the exploratory check above:
combined$percent.hb <- hb_check$percent.hb[match(colnames(combined), colnames(hb_check))]
#confirm the transfer worked (no NAs for cells present in both objects)
sum(is.na(combined$percent.hb))

mac_ifn <- subset(combined, cell_type == "Macrophages")
#does percent.hb differ between FLT and GC in Macrophages?
wilcox.test(percent.hb ~ condition, data = mac_ifn@meta.data)
#does Snca correlate with contamination level?
snca_expr <- GetAssayData(mac_ifn, layer = "data")["Snca", ]
cor.test(snca_expr, mac_ifn$percent.hb, method = "spearman")
#does Snca still differ by condition AFTER controlling for percent.hb?
mac_ifn$snca_expr <- snca_expr
summary(glm(snca_expr ~ condition + percent.hb, data = mac_ifn@meta.data))
#same for IFN score
summary(glm(IFN_score1 ~ condition + percent.hb, data = mac_ifn@meta.data))

#repeat across all cell types where Snca was significant in the pooled/per-celltype DEG results
cell_types_snca_all8 <- c("B cells", "T cells", "B cells (atypical/memory)",
                          "Macrophages", "B cells (naive)", "Innate B cells",
                          "Monocytes (patrolling)", "NK cells")
snca_contamination_check_all8 <- lapply(cell_types_snca_all8, function(ct) {
  sub_obj <- subset(combined, cell_type == ct)

  hb_test <- wilcox.test(percent.hb ~ condition, data = sub_obj@meta.data)

  snca_expr <- GetAssayData(sub_obj, layer = "data")["Snca", ]
  sub_obj$snca_expr <- snca_expr

  cor_test <- cor.test(snca_expr, sub_obj$percent.hb, method = "spearman")

  glm_fit <- glm(snca_expr ~ condition + percent.hb, data = sub_obj@meta.data)
  glm_summary <- summary(glm_fit)$coefficients

  data.frame(
    cell_type = ct,
    n_cells = ncol(sub_obj),
    hb_by_condition_p = hb_test$p.value,
    snca_hb_correlation = cor_test$estimate,
    condition_p_after_adjustment = glm_summary["conditionGC", "Pr(>|t|)"],
    condition_estimate = glm_summary["conditionGC", "Estimate"]
  )
})
snca_check_all8_df <- bind_rows(snca_contamination_check_all8)
print(snca_check_all8_df)

#How many significant DEGs does Macrophages actually have in the corrected object?
nrow(deg_results[["Macrophages"]] %>% filter(p_val_adj < 0.05))

#compare to other cell types to see if Macrophages is really the strongest
deg_summary %>%
  group_by(cell_type) %>%
  summarise(total_degs = sum(n_degs)) %>%
  arrange(desc(total_degs))

#Are interferon-stimulated genes actually enriched among Macrophage DEGs?
mac_degs_up <- deg_results[["Macrophages"]] %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  pull(gene)

#how many of the significantly upregulated FLT genes are IFN-response genes?
sum(mac_degs_up %in% ifn_genes)
length(mac_degs_up)
mac_degs_up[mac_degs_up %in% ifn_genes]  # actually list which ones

#repeated below with a fixed seed for a reproducible final value
wilcox.test(IFN_score1 ~ condition, data = mac_ifn@meta.data)

grep("^Hb[ab]-", rownames(mac_ifn), value = TRUE)  #should be empty if corrected object
"Erythrocytes" %in% unique(combined$cell_type)      #should be FALSE
nrow(mac_ifn)

#recompute IFN score with fixed seed for a clean test
set.seed(42)
combined <- AddModuleScore(combined, features = list(ifn_genes), name = "IFN_score")
mac_ifn <- subset(combined, cell_type == "Macrophages")
wilcox.test(IFN_score1 ~ condition, data = mac_ifn@meta.data)

#recompute the full per-cell-type IFN test with correction, using this freshly-seeded IFN_score1
ifn_score_test <- combined@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    FLT_mean = mean(IFN_score1[condition == "FLT"]),
    GC_mean  = mean(IFN_score1[condition == "GC"]),
    p_value  = wilcox.test(IFN_score1[condition == "FLT"],
                           IFN_score1[condition == "GC"])$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
print(ifn_score_test, n = 22)

#does the IFN module score itself support this claim?
ifn_score_test %>% filter(cell_type == "Macrophages")

b_cells_check <- subset(combined, cell_type == "B cells")

#does contamination differ between FLT/GC in B cells specifically?
wilcox.test(percent.hb ~ condition, data = b_cells_check@meta.data)

#does Snca correlate with contamination in B cells?
snca_expr_b <- GetAssayData(b_cells_check, layer = "data")["Snca", ]
cor.test(snca_expr_b, b_cells_check$percent.hb, method = "spearman")

#does condition survive after controlling for percent.hb?
b_cells_check$snca_expr <- snca_expr_b
summary(glm(snca_expr ~ condition + percent.hb, data = b_cells_check@meta.data))
summary(glm(IFN_score1 ~ condition + percent.hb, data = b_cells_check@meta.data))
ifn_expr_b <- b_cells_check$IFN_score1
cor.test(ifn_expr_b, b_cells_check$percent.hb, method = "spearman")

#Innate B cells
innate_b <- subset(combined, cell_type == "Innate B cells")
VlnPlot(innate_b, features = "IFN_score1", group.by = "condition",
        cols = c("FLT" = "#2166AC", "GC" = "#D6604D"), pt.size = 0.5) +
  labs(title = "Interferon Alpha Response Score", x = NULL, y = "IFN Module Score") +
  theme(legend.position = "none")

#confirm direction matches the claim (FLT > GC)
innate_b@meta.data %>%
  group_by(condition) %>%
  summarise(mean_ifn = mean(IFN_score1))

#confirm it survives the contamination check
summary(glm(IFN_score1 ~ condition + percent.hb, data = innate_b@meta.data))

#final poster-ready plots
VlnPlot(innate_b,
        features  = "IFN_score1",
        group.by  = "condition",
        cols      = c("FLT" = "#2166AC", "GC" = "#D6604D"),
        pt.size   = 0.5) +
  labs(title = "Interferon Alpha Response Score",
       x = NULL, y = "IFN Module Score") +
  theme(legend.position  = "none",
        axis.text.x      = element_text(size = 16),
        axis.text.y      = element_text(size = 14),
        axis.title.y     = element_text(size = 16),
        plot.title       = element_text(size = 18, face = "bold"))
mac_cells_final <- subset(combined, cell_type == "Macrophages")

VlnPlot(mac_cells_final,
        features  = "Snca",
        group.by  = "condition",
        cols      = c("FLT" = "#2166AC", "GC" = "#D6604D"),
        pt.size   = 0.5) +
  labs(title = "Snca Expression",
       x = NULL, y = "Expression Level") +
  theme(legend.position  = "none",
        axis.text.x      = element_text(size = 16),
        axis.text.y      = element_text(size = 14),
        axis.title.y     = element_text(size = 16),
        plot.title       = element_text(size = 18, face = "bold"))

#Do the top 5 upregulated and downregulated genes survive contamination?
up_genes <- c("Gpx1", "Ubb", "Prdx2", "Snca", "Slc25a37") 
gene_hb <- bind_rows(lapply(up_genes, function(g) {
  bind_rows(lapply(cell_types_snca_all8, function(ct) {
    s <- subset(combined, cell_type == ct)
    s$expr <- GetAssayData(s, layer = "data")[g, ]
    fit <- summary(glm(expr ~ condition + percent.hb + nCount_RNA, data = s@meta.data))$coefficients
    data.frame(gene = g, cell_type = ct,
               estimate = fit["conditionGC", "Estimate"], p = fit["conditionGC", "Pr(>|t|)"])
  }))
}))
gene_hb$p_adj <- p.adjust(gene_hb$p, method = "BH")
gene_hb %>% arrange(gene, p_adj)

down_genes <- c("Ebf1", "Bach2", "Foxp1", "Bank1", "Aff3")
down_hb <- bind_rows(lapply(down_genes, function(g) {
  bind_rows(lapply(cell_types_snca_all8, function(ct) {
    s <- subset(combined, cell_type == ct)
    s$expr <- GetAssayData(s, layer = "data")[g, ]
    fit <- summary(glm(expr ~ condition + percent.hb + nCount_RNA, data = s@meta.data))$coefficients
    data.frame(gene = g, cell_type = ct,
               estimate = fit["conditionGC", "Estimate"], p = fit["conditionGC", "Pr(>|t|)"])
  }))
}))
down_hb$p_adj <- p.adjust(down_hb$p, method = "BH")
down_hb %>% arrange(gene, p_adj)


#save a checkpoint
saveRDS(combined, "combined_final_with_degs.rds")


#check to see if poster findings survive sequencing depth correction check
mac_ifn <- subset(combined, cell_type == "Macrophages")
snca_expr <- GetAssayData(mac_ifn, layer = "data")["Snca", ]
mac_ifn$snca_expr <- snca_expr

summary(glm(snca_expr ~ condition + percent.hb + nCount_RNA, data = mac_ifn@meta.data))
summary(glm(IFN_score1 ~ condition + percent.hb + nCount_RNA, data = mac_ifn@meta.data))


cell_types_snca_all8 <- c("B cells", "T cells", "B cells (atypical/memory)",
                          "Macrophages", "B cells (naive)", "Innate B cells",
                          "Monocytes (patrolling)", "NK cells")

snca_contamination_check_all8_depth <- lapply(cell_types_snca_all8, function(ct) {
  sub_obj <- subset(combined, cell_type == ct)

  snca_expr <- GetAssayData(sub_obj, layer = "data")["Snca", ]
  sub_obj$snca_expr <- snca_expr

  glm_fit <- glm(snca_expr ~ condition + percent.hb + nCount_RNA, data = sub_obj@meta.data)
  glm_summary <- summary(glm_fit)$coefficients

  data.frame(
    cell_type = ct,
    n_cells = ncol(sub_obj),
    condition_p_after_adjustment = glm_summary["conditionGC", "Pr(>|t|)"],
    condition_estimate = glm_summary["conditionGC", "Estimate"]
  )
})
snca_check_all8_depth_df <- bind_rows(snca_contamination_check_all8_depth)
print(snca_check_all8_depth_df)


b_cells_check <- subset(combined, cell_type == "B cells")
summary(glm(IFN_score1 ~ condition + percent.hb + nCount_RNA, data = b_cells_check@meta.data))

innate_b <- subset(combined, cell_type == "Innate B cells")
summary(glm(IFN_score1 ~ condition + percent.hb + nCount_RNA, data = innate_b@meta.data))




#repeats the Snca and IFN contamination tests using only cells with little hemoglobin contamination
low_hb <- subset(combined, percent.hb < 2)
ncol(low_hb)
table(low_hb$sample, low_hb$condition)

snca_low <- bind_rows(lapply(cell_types_snca_all8, function(ct) {
  s <- subset(low_hb, cell_type == ct)
  s$expr <- GetAssayData(s, layer = "data")["Snca", ]
  fit <- summary(glm(expr ~ condition + percent.hb + nCount_RNA, data = s@meta.data))$coefficients
  data.frame(cell_type = ct, n_FLT = sum(s$condition == "FLT"), n_GC = sum(s$condition == "GC"),
             estimate = fit["conditionGC", "Estimate"], p = fit["conditionGC", "Pr(>|t|)"])
}))
snca_low$p_adj <- p.adjust(snca_low$p, method = "BH")
snca_low

ifn_low <- bind_rows(lapply(c("B cells", "Innate B cells", "T cells", "Monocytes (patrolling)", "B cells (atypical/memory)"), function(ct) {
  s <- subset(low_hb, cell_type == ct)
  fit <- summary(glm(IFN_score1 ~ condition + percent.hb + nCount_RNA, data = s@meta.data))$coefficients
  data.frame(cell_type = ct, n_FLT = sum(s$condition == "FLT"), n_GC = sum(s$condition == "GC"),
             estimate = fit["conditionGC", "Estimate"], p = fit["conditionGC", "Pr(>|t|)"])
}))
ifn_low

#re-test Snca and the IFN score for a FLT vs GC difference after adjusting for both Hb contamination (percent.hb) and depth (nCount_RNA)
#Macrophages: Snca and IFN with both covariates
mac_ifn <- subset(combined, cell_type == "Macrophages")
mac_ifn$snca_expr <- GetAssayData(mac_ifn, layer = "data")["Snca", ]
summary(glm(snca_expr ~ condition + percent.hb + nCount_RNA, data = mac_ifn@meta.data))$coefficients
summary(glm(IFN_score1 ~ condition + percent.hb + nCount_RNA, data = mac_ifn@meta.data))$coefficients

#Snca across all 8 cell types, side by side with the hb-only result
snca_hb_only <- bind_rows(lapply(cell_types_snca_all8, function(ct) {
  s <- subset(combined, cell_type == ct)
  s$expr <- GetAssayData(s, layer = "data")["Snca", ]
  fit <- summary(glm(expr ~ condition + percent.hb, data = s@meta.data))$coefficients
  data.frame(cell_type = ct, p_hb_only = fit["conditionGC", "Pr(>|t|)"],
             est_hb_only = fit["conditionGC", "Estimate"])
}))
snca_both <- bind_rows(lapply(cell_types_snca_all8, function(ct) {
  s <- subset(combined, cell_type == ct)
  s$expr <- GetAssayData(s, layer = "data")["Snca", ]
  fit <- summary(glm(expr ~ condition + percent.hb + nCount_RNA, data = s@meta.data))$coefficients
  data.frame(cell_type = ct, p_hb_depth = fit["conditionGC", "Pr(>|t|)"],
             est_hb_depth = fit["conditionGC", "Estimate"])
}))
left_join(snca_hb_only, snca_both, by = "cell_type")

#IFN in the three populations on the poster/conclusions
for (ct in c("B cells", "Innate B cells", "T cells", "Monocytes (patrolling)", "B cells (atypical/memory)")) {
  s <- subset(combined, cell_type == ct)
  cat("\n", ct, "\n")
  print(summary(glm(IFN_score1 ~ condition + percent.hb + nCount_RNA, data = s@meta.data))$coefficients["conditionGC", ])
}


#CellChat setup and dependencies
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(c("BiocNeighbors", "ComplexHeatmap"))
#install.packages("NMF")
#install.packages("circlize")
#install.packages("ggalluvial")
#if (!requireNamespace("devtools", quietly = TRUE))
  #install.packages("devtools")
#devtools::install_github("jinworks/CellChat")

library(CellChat)
library(Seurat)
library(patchwork)
library(NMF)
library(foreach)
library(ggalluvial)

options(stringsAsFactors = FALSE)
registerDoSEQ()  #force sequential NMF - avoids parallel worker errors

dir.create("plots/comparison", recursive = TRUE, showWarnings = FALSE)

#create a CellChat-safe version of cell_type with parentheses removed
combined$cell_type_cc <- gsub("\\s*\\(([^)]+)\\)", " \\1", combined$cell_type)
#this converts "Monocytes (patrolling)" -> "Monocytes patrolling", etc.
#check it worked correctly across all your cell types:
unique(combined$cell_type_cc)

#Add broad cell type layer (built on parenthesis-free cell_type_cc, for CellChat)
#confirm exact values before mapping - critical since names must match exactly
unique(combined$cell_type_cc)

combined$celltype_fine <- combined$cell_type_cc

combined$celltype_broad <- plyr::mapvalues(combined$celltype_fine,
                                     from = c("B cells", "B cells atypical/memory", "B cells naive", "B cells progenitor",
                                              "Innate B cells", "Plasma cells",
                                              "Dendritic cells", "Plasmacytoid DC",
                                              "Monocytes", "Monocytes inflammatory", "Monocytes patrolling",
                                              "Macrophages", "Monocyte-derived Macrophages",
                                              "Neutrophils", "Neutrophils immature", "Neutrophils inflammatory",
                                              "NK cells", "T cells", "T cells CD8+",
                                              "HSCs/Progenitors", "Proerythroblasts", "Platelets"),
                                     to   = c("B Cells", "B Cells", "B Cells", "B Cells",
                                              "B Cells", "B Cells",
                                              "Dendritic Cells", "Dendritic Cells",
                                              "Monocytes", "Monocytes", "Monocytes",
                                              "Macrophages", "Macrophages",
                                              "Neutrophils", "Neutrophils", "Neutrophils",
                                              "NK Cells", "T Cells", "T Cells",
                                              "HSCs/Progenitors", "Erythroid Cells", "Platelets")
)


#Reorder Clusters
Idents(combined) <- "celltype_broad"
broad_order <- c("HSCs/Progenitors", "Monocytes", "Macrophages", "Dendritic Cells",
                 "T Cells", "NK Cells", "Neutrophils", "B Cells", "Erythroid Cells", "Platelets")
combined$celltype_broad <- factor(combined$celltype_broad, levels = broad_order)
DimPlot(combined, group.by = "celltype_broad", label = TRUE) + ggtitle("Broad Cell Types")

Idents(combined) <- "celltype_fine"
fine_order <- c("HSCs/Progenitors",
                "Monocytes", "Monocytes inflammatory", "Monocytes patrolling",
                "Macrophages", "Monocyte-derived Macrophages",
                "Dendritic cells", "Plasmacytoid DC",
                "T cells", "T cells CD8+", "NK cells",
                "Neutrophils", "Neutrophils immature", "Neutrophils inflammatory",
                "B cells", "B cells naive", "B cells progenitor",
                "B cells atypical/memory", "Innate B cells", "Plasma cells",
                "Proerythroblasts", "Platelets")
combined$celltype_fine <- factor(combined$celltype_fine, levels = fine_order)
DimPlot(combined, group.by = "celltype_fine", label = TRUE) + ggtitle("Fine Cell Types")

#Part 1: FLT group
seurat_flt <- subset(combined, condition == "FLT")
table(seurat_flt$cell_type_cc)

#map sample -> samples column CellChat v2 expects (must be a factor)
seurat_flt$samples <- factor(seurat_flt$sample)
Idents(seurat_flt) <- "cell_type_cc"

cellchat_flt <- createCellChat(
  object   = seurat_flt,
  group.by = "cell_type_cc",
  assay    = "RNA"
)

#MOUSE database - mouse PBMC data
cellchat_flt@DB <- CellChatDB.mouse

cellchat_flt <- subsetData(cellchat_flt)
cellchat_flt <- identifyOverExpressedGenes(cellchat_flt, do.fast = FALSE)
cellchat_flt <- identifyOverExpressedInteractions(cellchat_flt)

cellchat_flt <- computeCommunProb(cellchat_flt, type = "triMean", nboot = 100)
cellchat_flt <- filterCommunication(cellchat_flt, min.cells = 10)
cellchat_flt <- computeCommunProbPathway(cellchat_flt)
cellchat_flt <- aggregateNet(cellchat_flt)
cellchat_flt <- netAnalysis_computeCentrality(cellchat_flt, slot.name = "netP")

saveRDS(cellchat_flt, file = "cellchat_FLT.rds")


#Part 2: GC group - identical pipeline
seurat_gc <- subset(combined, condition == "GC")
seurat_gc$samples <- factor(seurat_gc$sample)
Idents(seurat_gc) <- "cell_type_cc"

cellchat_gc <- createCellChat(
  object   = seurat_gc,
  group.by = "cell_type_cc",
  assay    = "RNA"
)
cellchat_gc@DB <- CellChatDB.mouse
cellchat_gc <- subsetData(cellchat_gc)
cellchat_gc <- identifyOverExpressedGenes(cellchat_gc, do.fast = FALSE)
cellchat_gc <- identifyOverExpressedInteractions(cellchat_gc)
cellchat_gc <- computeCommunProb(cellchat_gc, type = "triMean", nboot = 100)
cellchat_gc <- filterCommunication(cellchat_gc, min.cells = 10)
cellchat_gc <- computeCommunProbPathway(cellchat_gc)
cellchat_gc <- aggregateNet(cellchat_gc)
cellchat_gc <- netAnalysis_computeCentrality(cellchat_gc, slot.name = "netP")

saveRDS(cellchat_gc, file = "cellchat_GC.rds")

#Part 3: Comparative analysis - FLT vs GC
cellchat_flt <- readRDS("cellchat_FLT.rds")
cellchat_gc  <- readRDS("cellchat_GC.rds")

#order determines direction of every differential result below.
#query (index 2) - reference (index 1) = GC - FLT
#red = stronger in GC; blue = stronger in FLT
object.list <- list(
  FLT = cellchat_flt,
  GC  = cellchat_gc
)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

#confirm identical cell type labels between the two objects before trusting the merge
setdiff(levels(cellchat_flt@idents), levels(cellchat_gc@idents))
setdiff(levels(cellchat_gc@idents), levels(cellchat_flt@idents))

#total interaction number/strength - absolute values, no directionality
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
png("plots/comparison/01_interaction_count_vs_strength.png", width = 1200, height = 750, res = 150)
print(gg1 + gg2)
dev.off()

#per-condition circle plots
for (i in 1:length(object.list)) {
  fname <- paste0("plots/comparison/02_circle_", tolower(names(object.list)[i]), ".png")
  png(fname, width = 700, height = 700, res = 150)
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE, label.edge = FALSE,
    title.name   = paste("Number of Interactions -", names(object.list)[i])
  )
  dev.off()
}

#differential interaction circle plots (red = more in GC, blue = more in FLT per list order above)
#FLIGHT TO GROUND CONTROL !!!!!!!!!!!!!!!! red is decrease in FLT compared to GC and blue is increase in FLT compared to GC
png("plots/comparison/03a_diff_interactions_count.png", width = 700, height = 700, res = 150)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
dev.off()

png("plots/comparison/03b_diff_interactions_strength.png", width = 700, height = 700, res = 150)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
dev.off()

#differential heatmap - more systematic view across all cell-type pairs
gg1 <- netVisual_heatmap(cellchat_merged)
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
png("plots/comparison/04_differential_heatmaps.png", width = 1800, height = 750, res = 150)
print(gg1 + gg2)
dev.off()

#pathway-level comparison
png("plots/comparison/05_rankNet_stacked.png", width = 1050, height = 1200, res = 150)
rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
dev.off()

png("plots/comparison/06_rankNet_sidebyside.png", width = 1050, height = 1200, res = 150)
rankNet(cellchat_merged, mode = "comparison", stacked = FALSE)
dev.off()

#outgoing/incoming signaling heatmaps per condition - who changed sending vs receiving
png("plots/comparison/07_outgoing_heatmap_FLT.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing",
                                  title = "FLT — Outgoing Signaling", height = 12, color.heatmap = "GnBu")
dev.off()

png("plots/comparison/08_outgoing_heatmap_GC.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing",
                                  title = "GC — Outgoing Signaling", height = 12, color.heatmap = "GnBu")
dev.off()

png("plots/comparison/09_incoming_heatmap_FLT.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming",
                                  title = "FLT — Incoming Signaling", height = 12, color.heatmap = "OrRd")
dev.off()

png("plots/comparison/10_incoming_heatmap_GC.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming",
                                  title = "GC — Incoming Signaling", height = 12, color.heatmap = "OrRd")
dev.off()

#bubble plot - Macrophages as sender to your validated IFN-elevated populations
png("plots/comparison/11_bubble_macrophages_to_targets.png", width = 1500, height = 1200, res = 150)
p <- netVisual_bubble(
  cellchat_merged,
  sources.use = c("Macrophages"),
  targets.use = c("Innate B cells", "Monocytes patrolling", "B cells", "T cells"),
  comparison  = c(1, 2),
  angle.x     = 45,
  remove.isolate = FALSE
)
print(p)
dev.off()

saveRDS(cellchat_merged, file = "cellchat_merged_FLT_GC.rds")

list.files("plots/comparison", full.names = TRUE)


#Same as above - for broad cell type
dir.create("plots/comparison_broad", recursive = TRUE, showWarnings = FALSE)

#Broad-resolution CellChat: FLT
seurat_flt_broad <- subset(combined, condition == "FLT")
seurat_flt_broad$samples <- factor(seurat_flt_broad$sample)
Idents(seurat_flt_broad) <- "celltype_broad"

cellchat_flt_broad <- createCellChat(object = seurat_flt_broad, group.by = "celltype_broad", assay = "RNA")
cellchat_flt_broad@DB <- CellChatDB.mouse
cellchat_flt_broad <- subsetData(cellchat_flt_broad)
cellchat_flt_broad <- identifyOverExpressedGenes(cellchat_flt_broad, do.fast = FALSE)
cellchat_flt_broad <- identifyOverExpressedInteractions(cellchat_flt_broad)
cellchat_flt_broad <- computeCommunProb(cellchat_flt_broad, type = "triMean", nboot = 100)
cellchat_flt_broad <- filterCommunication(cellchat_flt_broad, min.cells = 10)
cellchat_flt_broad <- computeCommunProbPathway(cellchat_flt_broad)
cellchat_flt_broad <- aggregateNet(cellchat_flt_broad)
cellchat_flt_broad <- netAnalysis_computeCentrality(cellchat_flt_broad, slot.name = "netP")
saveRDS(cellchat_flt_broad, file = "cellchat_FLT_broad.rds")

#Broad-resolution CellChat: GC
seurat_gc_broad <- subset(combined, condition == "GC")
seurat_gc_broad$samples <- factor(seurat_gc_broad$sample)
Idents(seurat_gc_broad) <- "celltype_broad"

cellchat_gc_broad <- createCellChat(object = seurat_gc_broad, group.by = "celltype_broad", assay = "RNA")
cellchat_gc_broad@DB <- CellChatDB.mouse
cellchat_gc_broad <- subsetData(cellchat_gc_broad)
cellchat_gc_broad <- identifyOverExpressedGenes(cellchat_gc_broad, do.fast = FALSE)
cellchat_gc_broad <- identifyOverExpressedInteractions(cellchat_gc_broad)
cellchat_gc_broad <- computeCommunProb(cellchat_gc_broad, type = "triMean", nboot = 100)
cellchat_gc_broad <- filterCommunication(cellchat_gc_broad, min.cells = 10)
cellchat_gc_broad <- computeCommunProbPathway(cellchat_gc_broad)
cellchat_gc_broad <- aggregateNet(cellchat_gc_broad)
cellchat_gc_broad <- netAnalysis_computeCentrality(cellchat_gc_broad, slot.name = "netP")
saveRDS(cellchat_gc_broad, file = "cellchat_GC_broad.rds")


#Part 3 (broad): Comparative analysis - FLT vs GC
cellchat_flt_broad <- readRDS("cellchat_FLT_broad.rds")
cellchat_gc_broad  <- readRDS("cellchat_GC_broad.rds")

object.list_broad <- list(FLT = cellchat_flt_broad, GC = cellchat_gc_broad)
cellchat_merged_broad <- mergeCellChat(object.list_broad, add.names = names(object.list_broad))

setdiff(levels(cellchat_flt_broad@idents), levels(cellchat_gc_broad@idents))
setdiff(levels(cellchat_gc_broad@idents), levels(cellchat_flt_broad@idents))

#total interaction number/strength
gg1 <- compareInteractions(cellchat_merged_broad, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_merged_broad, show.legend = FALSE, group = c(1, 2), measure = "weight")
png("plots/comparison_broad/01_interaction_count_vs_strength.png", width = 1200, height = 750, res = 150)
print(gg1 + gg2)
dev.off()

#per-condition circle plots
for (i in 1:length(object.list_broad)) {
  fname <- paste0("plots/comparison_broad/02_circle_", tolower(names(object.list_broad)[i]), ".png")
  png(fname, width = 700, height = 700, res = 150)
  netVisual_circle(
    object.list_broad[[i]]@net$count,
    weight.scale = TRUE, label.edge = FALSE,
    title.name   = paste("Number of Interactions -", names(object.list_broad)[i])
  )
  dev.off()
}

#differential interaction circle plots
png("plots/comparison_broad/03a_diff_interactions_count.png", width = 700, height = 700, res = 150)
netVisual_diffInteraction(cellchat_merged_broad, weight.scale = TRUE)
dev.off()

png("plots/comparison_broad/03b_diff_interactions_strength.png", width = 700, height = 700, res = 150)
netVisual_diffInteraction(cellchat_merged_broad, weight.scale = TRUE, measure = "weight")
dev.off()

#differential heatmap
gg1 <- netVisual_heatmap(cellchat_merged_broad)
gg2 <- netVisual_heatmap(cellchat_merged_broad, measure = "weight")
png("plots/comparison_broad/04_differential_heatmaps.png", width = 1800, height = 750, res = 150)
print(gg1 + gg2)
dev.off()

#pathway-level comparison
png("plots/comparison_broad/05_rankNet_stacked.png", width = 1050, height = 1200, res = 150)
rankNet(cellchat_merged_broad, mode = "comparison", stacked = TRUE, do.stat = TRUE)
dev.off()

png("plots/comparison_broad/06_rankNet_sidebyside.png", width = 1050, height = 1200, res = 150)
rankNet(cellchat_merged_broad, mode = "comparison", stacked = FALSE)
dev.off()

#outgoing/incoming signaling heatmaps per condition
png("plots/comparison_broad/07_outgoing_heatmap_FLT.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list_broad[[1]], pattern = "outgoing",
                                  title = "FLT — Outgoing Signaling", height = 12, color.heatmap = "GnBu")
dev.off()

png("plots/comparison_broad/08_outgoing_heatmap_GC.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list_broad[[2]], pattern = "outgoing",
                                  title = "GC — Outgoing Signaling", height = 12, color.heatmap = "GnBu")
dev.off()

png("plots/comparison_broad/09_incoming_heatmap_FLT.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list_broad[[1]], pattern = "incoming",
                                  title = "FLT — Incoming Signaling", height = 12, color.heatmap = "OrRd")
dev.off()

png("plots/comparison_broad/10_incoming_heatmap_GC.png", width = 1200, height = 1000, res = 150)
netAnalysis_signalingRole_heatmap(object.list_broad[[2]], pattern = "incoming",
                                  title = "GC — Incoming Signaling", height = 12, color.heatmap = "OrRd")
dev.off()

#bubble plot - Macrophages as sender to broad target categories
png("plots/comparison_broad/11_bubble_macrophages_to_targets.png", width = 1500, height = 1200, res = 150)
p_broad <- netVisual_bubble(
  cellchat_merged_broad,
  sources.use = c("Macrophages"),
  targets.use = c("B Cells", "T Cells", "NK Cells", "Neutrophils"),
  comparison  = c(1, 2),
  angle.x     = 45,
  remove.isolate = FALSE
)
print(p_broad)
dev.off()

saveRDS(cellchat_merged_broad, file = "cellchat_merged_FLT_GC_broad.rds")
list.files("plots/comparison_broad", full.names = TRUE)

Layers(combined[["RNA"]])

library(Matrix)

#export counts matrix, gene names, cell metadata as separate files
writeMM(GetAssayData(combined, layer = "counts"), "counts.mtx")
write.csv(rownames(combined), "genes.csv", row.names = FALSE)
write.csv(combined@meta.data, "metadata.csv", row.names = TRUE)

list.files(pattern = "counts.mtx|genes.csv|metadata.csv")
