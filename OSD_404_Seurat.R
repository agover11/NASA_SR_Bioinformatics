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


#filter low quality cells based on QC metrics
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

#remove cluster 23 ( has only 37 cells, very small)
combined <- subset(combined, seurat_clusters != "23")

#FEATUREPLOTS FOR AMBIGUOUS CLUSTERS
#cluster 11 — dendritic cells
FeaturePlot(combined,
            features = c("Flt3", "Zbtb46", "Siglech", "Irf8"),
            reduction = "umap", ncol = 2)

#cluster 1 — monocytes
FeaturePlot(combined,
            features = c("Csf1r", "Cd14", "F13a1", "Kit"),
            reduction = "umap", ncol = 2)

#cluster 13 — B cells
FeaturePlot(combined,
            features = c("Ms4a1", "Cd19", "Pax5", "Ly6d"),
            reduction = "umap", ncol = 2)

#cluster 17 — macrophages (red pulp macrophages - investigate this further)
FeaturePlot(combined,
            features = c("Gata6", "Mrc1", "Adgre1", "Cx3cr1"),
            reduction = "umap", ncol = 2)

#cluster 22 — monocyte-derived macrophages
FeaturePlot(combined,
            features = c("Mertk", "Mrc1", "Adgre1", "Spic"),
            reduction = "umap", ncol = 2)

#cluster 3 — mixed progenitor population
FeaturePlot(combined,
            features = c("Mcpt8", "Ms4a2", "Cpa3", "Itga2b"),
            reduction = "umap", ncol = 2)


#cell type annotation based on marker genes, SingleR, and FeaturePlot validation
cell_types <- c(
  "0"  = "Erythrocytes",
  "1"  = "Monocytes",
  "2"  = "Monocytes",
  "3"  = "HSCs/Progenitors",
  "4"  = "B cells (naive)",
  "5"  = "T cells",
  "6"  = "B cells",
  "7"  = "Monocytes (patrolling)",
  "8"  = "Neutrophils",
  "9"  = "T cells (CD8+)",
  "10" = "Plasmacytoid DC",
  "11" = "Dendritic cells",
  "12" = "B cells (progenitor)",
  "13" = "B cells",
  "14" = "Plasma cells",
  "15" = "Neutrophils",
  "16" = "NK cells",
  "17" = "Macrophages",
  "18" = "Erythrocytes",
  "19" = "Platelets",
  "20" = "Neutrophils (inflammatory)",
  "21" = "Proerythroblasts",
  "22" = "Monocyte-derived Macrophages"
)

#apply cell type labels to each cell based on cluster assignment
combined$cell_type <- as.character(cell_types[as.character(combined$seurat_clusters)])

#CLUSTER 3 SUBCLUSTERING 
#basophil markers (Mcpt8, Ms4a2, Cpa3) also occur with progenitor markers
#indicating a mixed population — subcluster to resolve
cluster3 <- subset(combined, seurat_clusters == "3")
cluster3 <- FindNeighbors(cluster3, reduction = "harmony", dims = 1:30)
cluster3 <- FindClusters(cluster3, resolution = 0.3)

DimPlot(cluster3, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(cluster3, features = c("Mcpt8", "Mpo", "Angpt1", "Cd34"), ncol = 2)
FeaturePlot(cluster3, features = c("Kit", "Cd34", "Mpo", "Elane", "Mcpt8", "Ms4a2", "Cpa3"),
            ncol = 3)
table(cluster3$seurat_clusters)

#find markers for cluster3 subclusters
cluster3 <- JoinLayers(cluster3)
cluster3_markers <- FindAllMarkers(cluster3,
                                   only.pos        = TRUE,
                                   min.pct         = 0.25,
                                   logfc.threshold = 0.25)

#top markers for HSCs/MPPs, Granulocyte-monocyte progenitors, and Progenitors intermediate
cluster3_markers %>%
  filter(cluster %in% c("0","1", "2")) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  top_n(10, avg_log2FC) %>%
  print(n = 30)

# map subcluster identities back onto the parent object using cell barcodes
combined$hsc_subcluster <- NA
combined$hsc_subcluster[colnames(cluster3)] <- as.character(cluster3$seurat_clusters)
combined$cell_type[combined$hsc_subcluster == "0"] <- "HSCs/MPPs"
combined$cell_type[combined$hsc_subcluster == "1"] <- "Granulocyte-monocyte progenitors"
combined$cell_type[combined$hsc_subcluster == "2"] <- "Progenitors (intermediate)"

#confirm final cell type count
length(unique(combined$cell_type))

#ANNOTATED UMAP 
#view final annotated UMAP
DimPlot(combined, reduction = "umap", group.by = "cell_type",
        label = TRUE, repel = TRUE)

#UMAP formatted for poster
p <- DimPlot(combined, reduction = "umap", group.by = "cell_type",
             label = TRUE, label.size = 8, repel = TRUE,
             pt.size = 0.5) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.35, "cm"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("umap_annotated.png",
       plot = p,
       width = 16, height = 15.81, units = "in",
       dpi = 300, bg = "white")

#CANONICAL MARKER DOTPLOT 
#canonical mouse immune marker panel, grouped by lineage
#check against data first so DotPlot doesn't choke on missing genes
canonical_markers <- c(
  # erythroid / progenitor
  "Hba-a1", "Hba-a2", "Hbb-bs", "Klf1", "Gata1", "Tfrc", "Kit", "Cd34", "Mpo", "Elane",
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
  "NK cells",
  "T cells (CD8+)",
  "T cells",
  "Macrophages",
  "Monocyte-derived Macrophages",
  "Plasmacytoid DC",
  "Dendritic cells",
  "Neutrophils (inflammatory)",
  "Neutrophils",
  "Monocytes (patrolling)",
  "Monocytes",
  "Granulocyte-monocyte progenitors",
  "Progenitors (intermediate)",
  "HSCs/MPPs",
  "Proerythroblasts",
  "Platelets",
  "Erythrocytes"
)

#one representative marker per cell type (roughly), ordered for diagonal expression pattern
poster_markers <- c(
  "Alas2",    # erythrocytes
  "Itga2b",   # platelets
  "Car2",     # proerythroblasts
  "Mecom",    # HSCs/MPPs
  "Mpo",      # granulocyte-monocyte progenitors
  "F13a1",    # monocytes
  "Nr4a1",    # monocytes (patrolling)
  "S100a8",   # neutrophils
  "Camp",     # neutrophils (inflammatory)
  "Zbtb46",   # dendritic cells
  "Siglech",  # plasmacytoid DC
  "Mertk",    # monocyte-derived macrophages
  "Gata6",    # macrophages
  "Lef1",     # T cells
  "Cd8b1",    # T cells (CD8+)
  "Nkg7",     # NK cells
  "Vpreb3",   # B cells (progenitor)
  "Ighd",     # B cells (naive)
  "Sox5",     # B cells
  "Jchain"    # plasma cells
)

p_dotplot <- DotPlot(combined,
                     features = poster_markers,
                     group.by = "cell_type") +
  scale_y_discrete(limits = cell_type_order) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 14))

p_dotplot

ggsave("dotplot_markers.png",
       plot = p_dotplot,
       width = 11.27, height = 7.75, units = "in",
       dpi = 300, bg = "white")

#M1/M2 MACROPHAGE POLARIZATION SCORING 
m1_genes <- c("Nos2", "Il1b", "Tnf", "Il6", "Cxcl9", "Cxcl10", "Il12b", "Ptgs2", "Cd86", "Socs3")
m2_genes <- c("Arg1", "Mrc1", "Chil3", "Retnla", "Il10", "Cd163", "Il1rn", "Tgfb1", "Vegfa")
m1_genes <- m1_genes[m1_genes %in% rownames(combined)]
m2_genes <- m2_genes[m2_genes %in% rownames(combined)]

combined <- AddModuleScore(combined, features = list(m1_genes), name = "M1_score")
combined <- AddModuleScore(combined, features = list(m2_genes), name = "M2_score")

VlnPlot(combined, features = c("M1_score1", "M2_score1"),
        group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = m1_genes, group.by = "cell_type", pt.size = 0, ncol = 3)
VlnPlot(combined, features = m2_genes, group.by = "cell_type", pt.size = 0, ncol = 3)

#QC AND DOUBLET CHECKS 
#also compare basic QC metrics -- doublets often run higher on both
VlnPlot(combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        group.by = "cell_type", pt.size = 0)

#check whether "macrophage" cells also carry strong platelet/RBC signal
#a common sign of phagocytosed material or doublets rather than a clean identity
mac_cells <- subset(combined, cell_type == "Monocyte-derived Macrophages")
VlnPlot(mac_cells, features = c("Adgre1", "Pf4", "Itga2b", "Ppbp", "Hba-a1"), pt.size = 0)

#QC metrics per erythrocyte cluster to confirm they are not low-quality debris
combined@meta.data %>%
  filter(seurat_clusters %in% c("0", "18")) %>%
  group_by(seurat_clusters) %>%
  summarise(mean_nCount   = mean(nCount_RNA),
            mean_nFeature = mean(nFeature_RNA),
            mean_pct_mt   = mean(percent.mt),
            n             = n())

table(combined$seurat_clusters[combined$seurat_clusters %in% c("0", "18")],
      combined$sample[combined$seurat_clusters %in% c("0", "18")])

#DOUBLET DETECTION
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
#FO20 notably lower than other FLT mice (708 vs 2655-3965)
#GO14 and GO19 higher than other GC mice
table(combined$sample, combined$condition)
table(combined$cell_type, combined$condition)

#what fraction of total cells is this?
round(prop.table(table(combined$cell_type)) * 100, 2)

#normalize by total cells per condition to get proportions
#necessary bc GC has significantly more cells than FLT (14750 vs 10562)
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

ggsave("cell_type_composition.png",
       plot = p_composition,
       width = 11.27, height = 7.75, units = "in",
       dpi = 300, bg = "white")

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
#GC dominates visually due to higher cell count
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
#load MitoCarta3.0 mouse mitochondrial gene set for module scoring
mitocarta <- read_excel("/Users/lexigover/Downloads/Mouse.MitoCarta3.0.xls",
                        sheet = "A Mouse MitoCarta3.0")
mito_genes <- mitocarta$Symbol

#keep only genes that are actually in the Seurat object
mito_genes <- mito_genes[mito_genes %in% rownames(combined)]
cat("MitoCarta genes found in dataset:", length(mito_genes), "\n")


#run module scoring
combined <- AddModuleScore(
  combined,
  features = list(mito_genes),
  name = "MitoCarta_score"
)

VlnPlot(combined, 
        features = "MitoCarta_score1",
        group.by = "cell_type",
        pt.size = 0)

FeaturePlot(combined,
            features = "MitoCarta_score1",
            reduction = "umap")

VlnPlot(combined,
        features = "MitoCarta_score1",
        group.by = "cell_type",
        split.by = "condition",
        pt.size = 0)

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
#function to read a single gene set from a .gmt file
read_gmt <- function(file) {
  lines <- readLines(file)
  genes <- strsplit(lines, "\t")[[1]]
  genes <- genes[-(1:2)]  # remove gene set name and description
  return(genes)
}


#load all 5 gene sets
mitophagy_genes  <- read_gmt("~/Desktop/GOBP_MITOPHAGY.v2026.1.Mm.gmt")
cgas_sting_genes <- read_gmt("~/Desktop/GOBP_CGAS_STING_SIGNALING_PATHWAY.v2026.1.Mm.gmt")
nlrp3_genes      <- read_gmt("~/Desktop/GOBP_NLRP3_INFLAMMASOME_COMPLEX_ASSEMBLY.v2026.1.Mm.gmt")
aim2_genes       <- read_gmt("~/Desktop/GOBP_AIM2_INFLAMMASOME_COMPLEX_ASSEMBLY.v2026.1.Mm.gmt")
ifn_genes        <- read_gmt("~/Desktop/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2026.1.Mm.gmt")



#check how many genes were found in each set
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


#violin plots split by condition for each score
VlnPlot(combined, features = "Mitophagy_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "cGAS_STING_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "NLRP3_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "AIM2_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "IFN_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)

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


#IFN Alpha Response score in Macrophages (Wilcoxon padj = 2.31e-20)
VlnPlot(mac_ifn,
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


#Snca upregulated in FLT Macrophages (Wilcoxon p = 1.18e-22)
VlnPlot(mac_ifn,
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

#CELL-LEVEL DEG ANALYSIS 
#overall DEG analysis: all FLT vs all GC cells regardless of cell type
#note: confounded by cell type composition differences - use as overview?
Idents(combined) <- "condition"
all_degs <- FindMarkers(combined,
                        ident.1 = "FLT",
                        ident.2 = "GC",
                        min.pct = 0.25,
                        logfc.threshold = 0.25)
all_degs$gene <- rownames(all_degs)

#classify genes as up in FLT, up in GC, or not significant
all_degs <- all_degs %>%
  mutate(direction = case_when(
    avg_log2FC > 0.25  & p_val_adj < 0.05 ~ "Up in FLT",
    avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Up in GC",
    TRUE ~ "NS"
  ))

#top 10 upregulated and downregulated genes
top10_up   <- all_degs %>% filter(direction == "Up in FLT") %>%
  arrange(desc(avg_log2FC)) %>% head(10)
top10_down <- all_degs %>% filter(direction == "Up in GC") %>%
  arrange(avg_log2FC) %>% head(10)
top10 <- bind_rows(top10_up, top10_down)

top10

#overall volcano plot
ggplot(all_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(data = top10, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Up in FLT" = "#2166AC", "Up in GC" = "#D6604D", "NS" = "grey70")) +
  theme_classic() +
  labs(title = "Overall DEGs: Flight vs Ground Control",
       x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

#DEG analysis within each cell type with volcano plot
#note: treats each cell as independent — pseudoreplication at n=4
Idents(combined) <- "cell_type"
cell_types_to_test <- c(
  "Monocytes", "Monocytes (patrolling)", "Neutrophils",
  "Neutrophils (inflammatory)", "NK cells", "T cells",
  "T cells (CD8+)", "B cells (naive)", "B cells",
  "Macrophages", "Monocyte-derived Macrophages",
  "Dendritic cells", "Plasmacytoid DC"
)

deg_results <- list()
for (ct in cell_types_to_test) {
  cat("Running DEG for", ct, "\n")
  tryCatch({
    degs <- FindMarkers(combined,
                        ident.1      = "FLT",
                        ident.2      = "GC",
                        group.by     = "condition",
                        subset.ident = ct,
                        min.pct      = 0.25,
                        logfc.threshold = 0.25)
    degs$gene      <- rownames(degs)
    degs$cell_type <- ct
    degs <- degs %>%
      mutate(direction = case_when(
        avg_log2FC > 0.25  & p_val_adj < 0.05 ~ "Up in FLT",
        avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Up in GC",
        TRUE ~ "NS"
      ))
    deg_results[[ct]] <- degs
    
    top10_up   <- degs %>% filter(direction == "Up in FLT") %>%
      arrange(desc(avg_log2FC)) %>% head(10)
    top10_down <- degs %>% filter(direction == "Up in GC") %>%
      arrange(avg_log2FC) %>% head(10)
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
  filter(!cell_type %in% c("Neutrophils", "NK cells", "T cells (CD8+)")) %>%
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
#aggregate raw counts per cell type per mouse since mouse is the statistical unit (n=4)
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
  
  #require at least 6 pseudobulk columns (3 FLT + 3 GC minimum)
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
    dds  <- DESeqDataSetFromMatrix(countData = cluster_counts, colData = meta,
                                   design = ~ condition)
    keep <- rowSums(counts(dds) >= 5) >= 3
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


#summary of DESeq2 significant hits
all_deseq2_degs <- bind_rows(deg_results_deseq2)
deg_summary_deseq2 <- all_deseq2_degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cell_type, direction) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  arrange(cell_type)
print(deg_summary_deseq2, n = 50)


#Hbb-bs and Hba-a1 upregulated in FLT Monocyte-derived Macrophages (DESeq2 padj = 0.009)
VlnPlot(mac_cells,
        features = c("Hbb-bs", "Hba-a1"),
        split.by = "condition",
        pt.size  = 0,
        ncol     = 2)

#LIMMA-VOOM
#pseudobulk limma-voom per cell type
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
    keep           <- rowSums(cluster_counts >= 5) >= 3
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



#summary of limma significant hits
all_limma_degs <- bind_rows(deg_results_limma)
deg_summary_limma <- all_limma_degs %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(cell_type, direction) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  arrange(cell_type)
print(deg_summary_limma, n = 50)

#Cxcl13 and Serpinb2 downregulated in FLT Macrophages (limma padj < 0.001 and 0.005)
VlnPlot(mac_ifn,
        features = c("Cxcl13", "Serpinb2"),
        split.by = "condition",
        pt.size  = 0,
        ncol     = 2)

#compare DESeq2 and limma hits 
deseq2_sig <- all_deseq2_degs %>%
  filter(p_val_adj < 0.05) %>%
  dplyr::select(gene, cell_type, avg_log2FC, p_val_adj, direction)

limma_sig <- all_limma_degs %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::select(gene, cell_type, logFC, adj.P.Val, direction)

cat("\nDESeq2 significant genes:\n"); print(deseq2_sig)
cat("\nlimma significant genes:\n");  print(limma_sig)


#erythroid genes in monocyte-derived macrophages might reflect erythrophagocytosis not macrophage-intrinsic transcriptional changes under spaceflight ?
cat("\nerythroid gene check — Monocyte-derived Macrophages vs Macrophages:\n")
bind_rows(deg_results_deseq2) %>%
  filter(cell_type %in% c("Monocyte-derived Macrophages", "Macrophages"),
         gene %in% c("Hbb-bs", "Hba-a1")) %>%
  dplyr::select(gene, cell_type, avg_log2FC, p_val_adj)

bind_rows(deg_results_limma) %>%
  filter(cell_type %in% c("Monocyte-derived Macrophages", "Macrophages"),
         gene %in% c("Hbb-bs", "Hba-a1")) %>%
  dplyr::select(gene, cell_type, logFC, adj.P.Val)






#compute percent.hb directly 
is_hb <- grepl("^Hb[ab]-", rownames(combined))
combined[["percent.hb"]] <- PercentageFeatureSet(combined, features = rownames(combined)[is_hb])

mac_ifn <- subset(combined, cell_type == "Macrophages")

#does percent.hb even differ between FLT and GC in Macrophages? (the real risk factor)
wilcox.test(percent.hb ~ condition, data = mac_ifn@meta.data)

#does Snca correlate with contamination level?
snca_expr <- GetAssayData(mac_ifn, layer = "data")["Snca", ]
cor.test(snca_expr, mac_ifn$percent.hb, method = "spearman")

#does Snca still differ by condition AFTER controlling for percent.hb?
mac_ifn$snca_expr <- snca_expr
summary(glm(snca_expr ~ condition + percent.hb, data = mac_ifn@meta.data))

#same for IFN score
summary(glm(IFN_score1 ~ condition + percent.hb, data = mac_ifn@meta.data))


# make sure percent.hb is computed on the full object first 
cell_types_snca <- c("T cells", "B cells", "Macrophages", "B cells (naive)", 
                     "Monocytes (patrolling)", "NK cells")

snca_contamination_check <- lapply(cell_types_snca, function(ct) {
  sub_obj <- subset(combined, cell_type == ct)
  
  hb_test <- wilcox.test(percent.hb ~ condition, data = sub_obj@meta.data)
  
  snca_expr <- GetAssayData(sub_obj, layer = "data")["Snca", ]
  sub_obj$snca_expr <- snca_expr
  
  cor_test <- cor.test(snca_expr, sub_obj$percent.hb, method = "spearman")
  
  glm_fit <- glm(snca_expr ~ condition + percent.hb, data = sub_obj@meta.data)
  glm_summary <- summary(glm_fit)$coefficients
  
  data.frame(
    cell_type = ct,
    hb_by_condition_p = hb_test$p.value,
    snca_hb_correlation = cor_test$estimate,
    condition_p_after_adjustment = glm_summary["conditionGC", "Pr(>|t|)"],
    condition_estimate = glm_summary["conditionGC", "Estimate"]
  )
})

snca_contamination_check_df <- bind_rows(snca_contamination_check)
print(snca_contamination_check_df)