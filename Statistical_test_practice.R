
# statistical test of proportions between FLT and GC per cell type
# using two-proportion z-test with BH correction for multiple comparisons
# significant findings: NK cells, monocytes (patrolling), neutrophils (inflammatory) elevated in FLT
# erythrocytes, neutrophils (mature), B cells (naive) elevated in GC
prop_test <- combined@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = condition,
              values_from = c(n, total)) %>%
  rowwise() %>%
  mutate(p_value = prop.test(c(n_FLT, n_GC),
                             c(total_FLT, total_GC))$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
print(prop_test, n = 22)

table(combined$seurat_clusters)
