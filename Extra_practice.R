# fgsea on the full ranked DESeq2 gene list for a given cell type
# extracts pathway-level signal even when no individual gene survives BH
library(fgsea)

# example for monocyte-derived macrophages
mac_res <- deg_results_deseq2[["Monocyte-derived Macrophages"]]

ranks <- mac_res %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%
  tibble::deframe()

#already-loaded gene sets
pathway_list <- list(
  MitoCarta   = mito_genes,
  Mitophagy   = mitophagy_genes,
  cGAS_STING  = cgas_sting_genes,
  NLRP3       = nlrp3_genes,
  AIM2        = aim2_genes,
  IFN         = ifn_genes
)

set.seed(42)
gsea_res <- fgsea(pathways = pathway_list,
                  stats    = ranks,
                  minSize  = 5,
                  maxSize  = 500)

gsea_res %>% arrange(pval) %>% print()


# build DEG-derived gene sets from pseudobulk hits and pathway overlaps
# these score cells on genes that actually changed in the data

all_mitocarta_degs <- c("Gpx1", "Alas2", "Prdx2", "Slc25a37", "mt-Cytb",
                        "mt-Nd1", "mt-Nd4", "Ucp2", "Fech", "Prxl2a",
                        "Cat", "mt-Nd2", "Slc25a39", "Ghitm", "Fth1",
                        "Bnip3l", "Hagh", "Fis1", "Slc25a51", "Bcl2", "Arf5") %>%
  intersect(rownames(combined))
cat("MitoCarta DEG gene set:", length(all_mitocarta_degs), "genes\n")

ifn_deg_genes <- c("Isg20", "Isg15", "Rsad2", "Cd74", "Ptpn22") %>%
  intersect(rownames(combined))
cat("IFN DEG gene set:", length(ifn_deg_genes), "genes\n")

combined <- AddModuleScore(combined, features = list(all_mitocarta_degs), name = "MitoCarta_DEG_score")
combined <- AddModuleScore(combined, features = list(ifn_deg_genes),      name = "IFN_DEG_score")

VlnPlot(combined, features = "MitoCarta_DEG_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)
VlnPlot(combined, features = "IFN_DEG_score1",
        group.by = "cell_type", split.by = "condition", pt.size = 0)

FeaturePlot(combined, features = "MitoCarta_DEG_score1", reduction = "umap")
FeaturePlot(combined, features = "IFN_DEG_score1",       reduction = "umap")

#classify mitoDAMP-responsive cells:
# above median on both mitochondrial stress AND innate immune activation scores
mito_threshold <- median(combined$MitoCarta_DEG_score1)
ifn_threshold  <- median(combined$IFN_DEG_score1)

combined$mitoDAMP_state <- case_when(
  combined$MitoCarta_DEG_score1 > mito_threshold &
    combined$IFN_DEG_score1 > ifn_threshold ~ "mitoDAMP-responsive",
  TRUE ~ "homeostatic"
)

table(combined$mitoDAMP_state)
table(combined$cell_type, combined$mitoDAMP_state)

#proportion mitoDAMP-responsive per condition
combined@meta.data %>%
  group_by(condition, mitoDAMP_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = round(n / sum(n) * 100, 2))

#proportion per cell type per condition
combined@meta.data %>%
  group_by(cell_type, condition, mitoDAMP_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type, condition) %>%
  mutate(proportion = round(n / sum(n) * 100, 2)) %>%
  filter(mitoDAMP_state == "mitoDAMP-responsive") %>%
  arrange(desc(proportion)) %>%
  print(n = 40)

DimPlot(combined,
        reduction = "umap",
        group.by  = "mitoDAMP_state",
        cols      = c("mitoDAMP-responsive" = "#D62728", "homeostatic" = "grey85"),
        pt.size   = 0.3) +
  ggtitle("mitoDAMP-responsive vs Homeostatic cells")

DimPlot(combined,
        reduction = "umap",
        group.by  = "mitoDAMP_state",
        cols      = c("mitoDAMP-responsive" = "#D62728", "homeostatic" = "grey85"),
        split.by  = "condition",
        pt.size   = 0.3)



