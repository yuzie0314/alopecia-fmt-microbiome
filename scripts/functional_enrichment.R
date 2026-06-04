# =============================================================================
#  functional_enrichment.R
#  Functional enrichment analysis
#
#  資料：
#    MetaCyc pathway（全名格式：e.g. "glycolysis I (from glucose 6-phosphate)"）
#    EC number（純數字格式：e.g. "1.1.1.1"）
#
#  分析：
#    Part 1 — MetaCyc SuperPathway rollup（MaAsLin2 only）
#    Part 2 — fgsea on MetaCyc（MaAsLin2 主要 + DESeq2 sensitivity，各自獨立）
#    Part 3 — EC → KEGG module ORA
#    Part 4 — EC → KEGG module fgsea
#    Part 5 — 整合視覺化
#
#  主要針對：FMT-specific 比較（FMT vs placebo @ 6M）
#
#  前置條件：
#    - maas_metacyc_full, maas_ec_full（來自 pathway_analysis.R）
#    - res_pw_A（來自 pairwise_timepoint_da.R，optional）
#    - metacyc_raw, ec_raw
#    - meta_raw, col_treatment, theme_fmt
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

library(tidyverse)
library(ggrepel)
library(patchwork)

if (!requireNamespace("fgsea", quietly = TRUE))
  BiocManager::install("fgsea")
library(fgsea)

if (!requireNamespace("KEGGREST", quietly = TRUE))
  BiocManager::install("KEGGREST")
library(KEGGREST)

OUT_DIR <- file.path(getwd(), "outputs", "functional_enrichment")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_enrich_fig <- function(plot, filename, w = 10, h = 7) {
  if (is.null(plot)) return(invisible())
  ggsave(file.path(OUT_DIR, filename), plot = plot,
         width = w, height = h, dpi = 300, bg = "white", device = cairo_pdf)
  cat(sprintf("  Saved: %s\n", filename))
}


# =============================================================================
#  PART 1 — MetaCyc SuperPathway rollup
#  使用官方 MetaCyc 層級結構
#
#  檔案：
#    map_metacyc-pwy_name_txt.gz  — no header, tab: PWY-ID \t full name
#    map_metacyc-pwy_lineage.tsv  — no header, tab: PWY-ID \t lineage（一對多）
#    lineage 格式：Biosynthesis|Cofactor-Biosynthesis|...（| 分隔層級）
# =============================================================================

cat("\n── Part 1: MetaCyc SuperPathway rollup ─────────────────────────────────\n")

# ── 1a. 讀取 mapping 檔案 ─────────────────────────────────────────────────────

pwy_name_map <- read_tsv(
  "../emulsion/full_mapping_v4_alpha/map_metacyc-pwy_name.txt.gz",
  col_names      = c("pwy_id", "pwy_name"),
  show_col_types = FALSE
) %>%
  mutate(pwy_name = str_trim(pwy_name),
         pwy_name = str_remove_all(pwy_name, "&"),
         pwy_name = str_remove_all(pwy_name, ";"))

pwy_lineage_raw <- read_tsv(
  "../emulsion/map_metacyc-pwy_lineage.tsv",
  col_names      = c("pwy_id", "lineage"),
  show_col_types = FALSE
)

cat(sprintf("Name map:    %d entries\n", nrow(pwy_name_map)))
cat(sprintf("Lineage map: %d entries, %d unique PWY IDs\n",
            nrow(pwy_lineage_raw), n_distinct(pwy_lineage_raw$pwy_id)))

# ── 1b. 解析層級 ──────────────────────────────────────────────────────────────
# top_level    = 第一層（最高，e.g. "Biosynthesis"、"Degradation"）
# second_level = 前兩層（e.g. "Biosynthesis > Cofactor-Biosynthesis"）

pwy_lineage <- pwy_lineage_raw %>%
  mutate(
    top_level    = str_extract(lineage, "^[^|]+"),
    second_level = str_extract(lineage, "^[^|]+\\|[^|]+") %>%
      str_replace("\\|", " > ")
  )

# 每個 PWY-ID 取第一筆 top_level（多筆時保留全部供後續用）/ second_level: 56
pwy_top <- pwy_lineage %>%
  group_by(pwy_id) %>%
  summarise(
    top_level    = first(top_level),
    all_lineages = paste(lineage, collapse = " || "),
    .groups      = "drop"
  )

cat(sprintf("Top-level categories: %s\n",
            paste(sort(unique(pwy_lineage$top_level)), collapse = ", ")))

# ── 1c. MaAsLin2 ranked list + 對應 category ─────────────────────────────────

ranked_pw <- maas_metacyc_full %>% mutate(Pathway = str_replace_all(Pathway,"\\.", "-"),
                                          Pathway = str_replace(Pathway, "^X(?=[0-9])", "")) %>%
  filter(metadata == "Treatment", !is.na(coef)) %>%
  transmute(
    pwy_id  = Pathway,
    coef    = coef,
    t       = coef/stderr,
    pval    = pval,
    qval    = qval
  ) %>%
  # pwy_id → top_level
  left_join(pwy_top %>% dplyr::select(pwy_id, top_level, all_lineages),
            by = "pwy_id") %>%
  mutate(category = if_else(is.na(top_level), "Unmapped", top_level)) %>%
  arrange(desc(abs(coef)))

n_mapped   <- sum(!is.na(ranked_pw$pwy_id))
n_unmapped <- nrow(ranked_pw) - n_mapped
cat(sprintf("Mapping: %d/%d (%.1f%%)  unmapped: %d\n",
            n_mapped, nrow(ranked_pw),
            100 * n_mapped / nrow(ranked_pw), n_unmapped))

if (n_unmapped > 0) {
  cat("First 5 unmapped:\n")
  print(ranked_pw %>% filter(is.na(pwy_id)) %>%
          head(5) %>% dplyr::select(feature_name))
}

write.csv(ranked_pw,
          file.path(OUT_DIR, "metacyc_ranked_with_category.csv"),
          row.names = FALSE)

# ── 1d. Category summary ──────────────────────────────────────────────────────

category_summary <- ranked_pw %>%
  filter(category != "Unmapped") %>%
  group_by(category) %>%
  summarise(
    n_pathways       = n(),
    n_up             = sum(coef > 0),
    n_down           = sum(coef < 0),
    n_sig_up         = sum(coef > 0 & pval < 0.05),
    n_sig_down       = sum(coef < 0 & pval < 0.05),
    mean_coef        = round(mean(coef),   4),
    median_coef      = round(median(coef), 4),
    enrichment_score = round(mean(coef) * log10(n() + 1), 4),
    .groups          = "drop"
  ) %>%
  arrange(desc(enrichment_score))

cat("\nCategory summary:\n")
print(category_summary, n = 25)

write.csv(category_summary,
          file.path(OUT_DIR, "metacyc_category_summary.csv"),
          row.names = FALSE)

# ── 1e. 視覺化 ────────────────────────────────────────────────────────────────

p_category_lollipop <- ggplot(
  category_summary %>% filter(n_pathways >= 3),
  aes(x = reorder(category, mean_coef), y = mean_coef,
      colour = mean_coef > 0)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(xend = reorder(category, mean_coef), yend = 0),
               linewidth = 0.8) +
  geom_point(aes(size = n_pathways)) +
  geom_text(aes(
    label = sprintf("n=%d \n(↑%d ↓%d)", n_pathways, n_sig_up, n_sig_down),
    hjust = ifelse(mean_coef > 0, -0.1, 1.1)
  ), size = 2.8, colour = "grey30") +
  coord_flip() +
  scale_colour_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "#2166AC"),
                      guide  = "none") +
  scale_size_continuous(range = c(3, 8), name = "# pathways") +
  labs(
    title    = "MetaCyc functional category rollup (official hierarchy)",
    subtitle = "FMT vs Placebo @ 6M  |  label: total (↑↓ p<0.05)",
    x        = NULL,
    y        = "Mean coefficient (positive = enriched in FMT)"
  ) +
  theme_minimal(base_family = "sans")+
  theme_fmt

top_per_cat <- ranked_pw %>%
  filter(category != "Unmapped") %>%
  group_by(category) %>%
  slice_max(abs(coef), n = 3) %>%
  ungroup() %>%
  mutate(
    pw_short = str_trunc(pwy_id, 50),
    sig_lab  = if_else(pval < 0.05, "*", "") # use pval instead！！！
  )

p_category_heatmap <- ggplot(
  top_per_cat,
  aes(x = "FMT vs Placebo", y = reorder(pw_short, coef), fill = coef)
) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_lab), colour = "white",
            fontface = "bold", size = 3.5) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, name = "Coefficient"
  ) +
  labs(
    title    = "Top 3 pathways per MetaCyc category (official hierarchy)",
    subtitle = "* p<0.05  |  FMT vs Placebo @ 6M",
    x = NULL, y = NULL
  ) +
  theme_fmt +
  theme(
    axis.text.y     = element_text(size = 7),
    strip.text.y    = element_text(size = 7, angle = 0, hjust = 0),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "right"
  )

save_enrich_fig(p_category_lollipop, "fig1a_category_lollipop.pdf", w = 11, h = 7)
save_enrich_fig(p_category_heatmap,  "fig1b_category_heatmap.pdf",  w = 13, h = 16)

metacyc_maaslin2_enrichm<- p_category_heatmap@data %>% select(pwy_id, coef, pval, qval, all_lineages) %>% inner_join(pwy_name_map, by = "pwy_id")
write.csv(metacyc_maaslin2_enrichm,
          file.path(OUT_DIR, "metacyc_maaslin2_enrichm.csv"),
          row.names = FALSE)
          
# ── 1f. Gene sets 供 Part 2 fgsea 使用 ───────────────────────────────────────

# Top-level gene sets
category_sets_official <- ranked_pw %>%
  filter(category != "Unmapped") %>%
  group_by(category) %>%
  summarise(pathways = list(pwy_id), .groups = "drop") %>%
  deframe()
category_sets_official <- category_sets_official[
  sapply(category_sets_official, length) >= 3
]

# Second-level gene sets（更細分類）
category_sets_level2 <- ranked_pw %>%
  filter(!is.na(pwy_id)) %>%
  left_join(
    pwy_lineage %>% dplyr::select(pwy_id, second_level) %>% distinct(),
    by = "pwy_id"
  ) %>%
  filter(!is.na(second_level)) %>%
  group_by(second_level) %>%
  summarise(pathways = list(pwy_id), .groups = "drop") %>%
  deframe()
category_sets_level2 <- category_sets_level2[
  sapply(category_sets_level2, length) >= 3
]

cat(sprintf("Gene sets — top-level: %d  second-level: %d\n",
            length(category_sets_official),
            length(category_sets_level2)))


# =============================================================================
#  PART 2 — fgsea on MetaCyc
#  MaAsLin2（主要）和 DESeq2（sensitivity）各自獨立跑，最後比較 NES 方向
# =============================================================================

cat("\n── Part 2: fgsea on MetaCyc ────────────────────────────────────────────\n")

# Gene sets：每個功能大類包含的 pathway 名稱（至少 3 個）
category_sets <- ranked_pw %>%
  filter(category != "Other") %>%
  group_by(category) %>%
  summarise(pathways = list(pwy_id), .groups = "drop") %>%
  deframe()
category_sets <- category_sets[sapply(category_sets, length) >= 3]

cat(sprintf("Gene sets: %d categories\n", length(category_sets)))

# ── 2a. fgsea — MaAsLin2（主要）─────────────────────────────────────────────

ranked_vec_maas <- ranked_pw %>%
  arrange(desc(t)) %>%
  { setNames(.$t, .$pwy_id) }

set.seed(42)
fgsea_maas <- fgsea(
  pathways = category_sets,
  stats    = ranked_vec_maas,
  minSize  = 3,
  maxSize  = 500,
  nperm    = 1000
) %>%
  as_tibble() %>%
  arrange(pval) %>%
  mutate(
    method    = "MaAsLin2",
    sig       = case_when(padj < 0.05 ~ "*", padj < 0.25 ~ "†", TRUE ~ ""),
    direction = if_else(NES > 0, "Enriched in FMT", "Depleted in FMT")
  )

cat("fgsea (MaAsLin2):\n")
print(fgsea_maas %>% dplyr::select(pathway, NES, pval, padj, sig, size))

# ── 2b. fgsea — DESeq2（sensitivity，若有）───────────────────────────────────

fgsea_deseq2 <- NULL

if (exists("res_pw_A") && !is.null(res_pw_A)) {
  # DESeq2 pathway result: 欄位為 Pathway, log2FoldChange, pvalue, padj
  ranked_vec_deseq2 <- res_pw_A %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(log2FoldChange)) %>%
    { setNames(.$log2FoldChange, .$Pathway) }
  
  # 只保留 category_sets 裡有的 pathway
  ranked_vec_deseq2 <- ranked_vec_deseq2[
    names(ranked_vec_deseq2) %in% unlist(category_sets)
  ]
  
  if (length(ranked_vec_deseq2) >= 10) {
    set.seed(42)
    fgsea_deseq2 <- fgsea(
      pathways = category_sets,
      stats    = ranked_vec_deseq2,
      minSize  = 3,
      maxSize  = 500,
      nperm    = 1000
    ) %>%
      as_tibble() %>%
      arrange(pval) %>%
      mutate(
        method    = "DESeq2",
        sig       = case_when(padj < 0.05 ~ "*", padj < 0.25 ~ "†", TRUE ~ ""),
        direction = if_else(NES > 0, "Enriched in FMT", "Depleted in FMT")
      )
    cat("fgsea (DESeq2 sensitivity):\n")
    print(fgsea_deseq2 %>% dplyr::select(pathway, NES, pval, padj, sig, size))
  } else {
    cat("DESeq2 ranked list too short for fgsea, skipping\n")
  }
}

# ── 2c. NES concordance（MaAsLin2 vs DESeq2）─────────────────────────────────

fgsea_concordance <- NULL
if (!is.null(fgsea_deseq2)) {
  fgsea_concordance <- fgsea_maas %>%
    dplyr::select(pathway, NES_maas = NES, padj_maas = padj, sig_maas = sig) %>%
    left_join(
      fgsea_deseq2 %>%
        dplyr::select(pathway, NES_deseq2 = NES, padj_deseq2 = padj, sig_deseq2 = sig),
      by = "pathway"
    ) %>%
    mutate(
      concordant = case_when(
        is.na(NES_deseq2)                    ~ "DESeq2 NA",
        sign(NES_maas) == sign(NES_deseq2)   ~ "Concordant",
        TRUE                                  ~ "Discordant"
      )
    )
  
  cat("\nNES concordance (MaAsLin2 vs DESeq2):\n")
  print(fgsea_concordance %>%
          dplyr::select(pathway, NES_maas, NES_deseq2, concordant))
}

# ── 2d. 儲存 ─────────────────────────────────────────────────────────────────

write.csv(fgsea_maas %>% dplyr::select(-leadingEdge),
          file.path(OUT_DIR, "fgsea_maaslin2_metacyc_tval.csv"),
          row.names = FALSE)
if (!is.null(fgsea_deseq2))
  write.csv(fgsea_deseq2 %>% dplyr::select(-leadingEdge),
            file.path(OUT_DIR, "fgsea_deseq2_metacyc.csv"),
            row.names = FALSE)
if (!is.null(fgsea_concordance))
  write.csv(fgsea_concordance,
            file.path(OUT_DIR, "fgsea_concordance_maas_deseq2.csv"),
            row.names = FALSE)

# Leading edge（MaAsLin2）
leading_edge_df <- fgsea_maas %>%
  dplyr::select(pathway, NES, padj, leadingEdge) %>%
  unnest(leadingEdge) %>%
  rename(pwy_id = leadingEdge) %>%
  left_join(
    ranked_pw %>% dplyr::select(pwy_id, coef, qval),
    by = "pwy_id"
  ) %>%
  arrange(padj, desc(abs(coef)))

write.csv(leading_edge_df,
          file.path(OUT_DIR, "fgsea_leading_edge_maaslin2_tval.csv"),
          row.names = FALSE)

# ── 2e. 視覺化 ────────────────────────────────────────────────────────────────

# MaAsLin2 NES bar
p_fgsea_maas <- ggplot(fgsea_maas,
                       aes(x = reorder(pathway, NES), y = NES,
                           fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(sig, sprintf(" (n=%d)", size))),
            hjust = if_else(fgsea_maas$NES > 0, -0.1, 1.1),
            size  = 3, colour = "grey30") +
  geom_hline(yintercept = 0, colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = c("Enriched in FMT" = "#D6604D",
                               "Depleted in FMT" = "#2166AC")) +
  labs(
    title    = "fgsea: MetaCyc functional categories (MaAsLin2)",
    subtitle = "FMT vs Placebo @ 6M  |  * padj<0.05  † padj<0.25",
    x = NULL, y = "NES", fill = NULL
  ) +
  theme_fmt + theme(legend.position = "bottom")

save_enrich_fig(p_fgsea_maas, "fig2a_fgsea_maaslin2_tval.pdf", w = 11, h = 7)

# Concordance scatter（MaAsLin2 NES vs DESeq2 NES）
if (!is.null(fgsea_concordance)) {
  p_fgsea_concordance <- ggplot(
    fgsea_concordance %>% filter(concordant != "DESeq2 NA"),
    aes(x = NES_maas, y = NES_deseq2, colour = concordant)
  ) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(aes(label = pathway), size = 3,
                    segment.colour = "grey60", max.overlaps = 12) +
    scale_colour_manual(
      values = c(Concordant = "#1D9E75", Discordant = "#D85A30")
    ) +
    labs(
      title    = "fgsea NES concordance: MaAsLin2 vs DESeq2",
      subtitle = "Categories above/below dashed line = same NES direction",
      x        = "NES (MaAsLin2)",
      y        = "NES (DESeq2)",
      colour   = NULL
    ) +
    theme_fmt
  
  save_enrich_fig(p_fgsea_concordance, "fig2b_fgsea_concordance.pdf",
                  w = 9, h = 7)
}

# ko_ec <- read_rds(file = "../emulsion/ko_ec_cache.rds") %>% head() # df: map ec to ko ontology ids
# read_rds(file = "../emulsion/ko_module_cache.rds") %>% head() # list: ko ontology ids to multiple module ids
# ko_pw <-read_rds(file = "../emulsion/ko_pathway_cache.rds") %>% head() # list: ko ontology ids to pathway ids
# =============================================================================
#  PART 3 — EC → KEGG module ORA
#  使用 ko_ec_cache.rds 和 ko_module_cache.rds
#
#  ko_ec_cache.rds    : data.frame, cols = ko（KO ID）+ ec（EC number）
#  ko_module_cache.rds: list, names = KO ID, values = module ID vector
#
#  流程：
#    EC（from maas_ec_full）→ KO（via ko_ec）→ module（via ko_module）
# =============================================================================

cat("\n── Part 3: EC → KEGG module ORA ────────────────────────────────────────\n")

# ── 3a. 載入 RDS cache ───────────────────────────────────────────────────────

ko_ec     <- readRDS("../emulsion/ko_ec_cache.rds")       # data.frame: ko, ec
ko_module <- readRDS("../emulsion/ko_module_cache.rds")   # list: KO → module vector

cat(sprintf("ko_ec:     %d rows, %d unique KO, %d unique EC\n",
            nrow(ko_ec), n_distinct(ko_ec$ko), n_distinct(ko_ec$ec)))
cat(sprintf("ko_module: %d KO entries\n", length(ko_module)))

# ── 3b. 建立 EC → module mapping ─────────────────────────────────────────────
# EC → KO（一對多）→ module（一對多）
# 結果：EC → unique module vector

ec_to_ko <- ko_ec %>%
  dplyr::select(ec, ko) %>%
  filter(!is.na(ec), !is.na(ko)) %>%
  distinct()

# KO → module（把 list 轉成 data.frame）
ko_to_module <- tibble(
  ko     = names(ko_module),
  module = ko_module
) %>%
  unnest(module) %>%
  filter(!is.na(module))

# EC → module（經由 KO bridge）
ec_to_module <- ec_to_ko %>%
  left_join(ko_to_module, by = "ko", relationship = "many-to-many") %>%
  filter(!is.na(module)) %>%
  dplyr::select(ec, module) %>%
  distinct()

cat(sprintf("EC→module links: %d  (unique EC: %d, unique module: %d)\n",
            nrow(ec_to_module),
            n_distinct(ec_to_module$ec),
            n_distinct(ec_to_module$module)))

# module → EC sets（for ORA / fgsea）
module_ec_sets <- ec_to_module %>%
  group_by(module) %>%
  summarise(ec_list = list(unique(ec)), .groups = "drop") %>%
  deframe()

# 至少 3 個 EC 才納入
module_ec_sets <- module_ec_sets[sapply(module_ec_sets, length) >= 3]
cat(sprintf("Modules with >= 3 EC: %d\n", length(module_ec_sets))) # Modules with >= 3 EC: 248

# ── 3c. 顯著 EC list（FMT，雙向）────────────────────────────────────────────

extract_ec <- function(maas_res, direction = c("up", "down", "both"),
                       pval_cut = 0.05) {
  direction <- match.arg(direction)
  df <- maas_res %>% filter(metadata == "Treatment", pval < pval_cut)
  if (direction == "up")   df <- df %>% filter(coef > 0)
  if (direction == "down") df <- df %>% filter(coef < 0)
  df %>%
    pull(Pathway) %>%
    str_extract("\\d+\\.\\d+\\.\\d+\\.\\d+") %>%
    na.omit() %>% unique()
}

if (exists("maas_ec_full")) {
  ec_sig_up   <- extract_ec(maas_ec_full, "up")
  ec_sig_down <- extract_ec(maas_ec_full, "down")
  ec_sig_both <- extract_ec(maas_ec_full, "both")
} else {
  message("maas_ec_full not found")
  ec_sig_up <- ec_sig_down <- ec_sig_both <- character(0)
}

# Background = 所有 EC 在 ko_ec 中有 mapping 的（有功能注解的）
ec_background <- intersect(
  rownames(ec_raw) %>% str_extract("\\d+\\.\\d+\\.\\d+\\.\\d+") %>% na.omit(),
  unique(ko_ec$ec)
)

cat(sprintf("Sig EC — up: %d  down: %d  both: %d\n",
            length(ec_sig_up), length(ec_sig_down), length(ec_sig_both)))
cat(sprintf("Background EC (with KO mapping): %d\n", length(ec_background)))

# ── 3d. ORA（Fisher's exact test）────────────────────────────────────────────

run_ora <- function(sig_list, background, gene_sets, min_overlap = 2) {
  purrr::map_dfr(names(gene_sets), function(mod) {
    mod_genes <- gene_sets[[mod]]
    n_overlap <- length(intersect(sig_list, mod_genes))
    if (n_overlap < min_overlap) return(NULL)
    mat <- matrix(c(
      n_overlap,
      length(sig_list)  - n_overlap,
      length(mod_genes) - n_overlap,
      length(background) - length(sig_list) - length(mod_genes) + n_overlap
    ), nrow = 2)
    ft <- fisher.test(mat, alternative = "greater")
    tibble(
      module     = mod,
      n_module   = length(mod_genes),
      n_overlap  = n_overlap,
      odds_ratio = ft$estimate,
      pval       = ft$p.value,
      overlap_ec = paste(intersect(sig_list, mod_genes), collapse = ";")
    )
  }) %>%
    { if (nrow(.) > 0) mutate(., padj = p.adjust(pval, method = "BH")) %>% arrange(padj)
      else tibble(module=character(), n_module=integer(), n_overlap=integer(),
                  odds_ratio=numeric(), pval=numeric(), overlap_ec=character(), padj=numeric()) }
}

# 雙向各跑一次
ora_up   <- run_ora(ec_sig_up,   ec_background, module_ec_sets)
ora_down <- run_ora(ec_sig_down, ec_background, module_ec_sets)

cat(sprintf("ORA Up   — padj<0.05: %d  padj<0.25: %d\n",
            sum(ora_up$padj   < 0.05, na.rm = TRUE),
            sum(ora_up$padj   < 0.25, na.rm = TRUE)))
cat(sprintf("ORA Down — padj<0.05: %d  padj<0.25: %d\n",
            sum(ora_down$padj < 0.05, na.rm = TRUE),
            sum(ora_down$padj < 0.25, na.rm = TRUE)))

ora_res <- bind_rows(
  ora_up   %>% mutate(direction = "Up"),
  ora_down %>% mutate(direction = "Down")
)

write.csv(ora_res, file.path(OUT_DIR, "ora_ec_kegg_module.csv"),
          row.names = FALSE)


# =============================================================================
#  PART 4 — EC → KEGG module + pathway fgsea
#  使用 ko_module_cache.rds 和 ko_pathway_cache.rds
# =============================================================================

cat("\n── Part 4: EC fgsea on KEGG modules + pathways ─────────────────────────\n")

ko_pathway <- readRDS("../emulsion/ko_pathway_cache.rds")  # list: KO → pathway vector

cat(sprintf("ko_pathway: %d KO entries\n", length(ko_pathway)))

# ── 4a. 建立 module + pathway gene sets（EC 為單位）─────────────────────────

# pathway → EC sets（同 ec_to_module 的邏輯）
ko_to_pathway <- tibble(
  ko      = names(ko_pathway),
  pathway = ko_pathway
) %>%
  unnest(pathway) %>%
  filter(!is.na(pathway))

pathway_ec_sets <- ec_to_ko %>%
  left_join(ko_to_pathway, by = "ko", relationship = "many-to-many") %>%
  filter(!is.na(pathway)) %>%
  dplyr::select(ec, pathway) %>%
  distinct() %>%
  group_by(pathway) %>%
  summarise(ec_list = list(unique(ec)), .groups = "drop") %>%
  deframe()

pathway_ec_sets <- pathway_ec_sets[sapply(pathway_ec_sets, length) >= 3]
cat(sprintf("Pathways with >= 3 EC: %d\n", length(pathway_ec_sets)))

# ── 4b. Ranked EC list（MaAsLin2 coefficient）────────────────────────────────

fgsea_ec_res         <- NULL
fgsea_ec_pathway_res <- NULL

if (exists("maas_ec_full")) {
  ec_ranked_vec <- maas_ec_full %>%
    filter(metadata == "Treatment", !is.na(coef)) %>%
    mutate(ec_clean = str_extract(Pathway, "\\d+\\.\\d+\\.\\d+\\.\\d+"),
           t = coef/stderr) %>%
    filter(!is.na(ec_clean)) %>%
    group_by(ec_clean) %>%
    summarise(coef = mean(coef),tval = mean(t), .groups = "drop") %>%
    arrange(desc(tval)) %>%
    { setNames(.$tval, .$ec_clean) }
  
  cat(sprintf("Ranked EC: %d entries\n", length(ec_ranked_vec)))
  
  # ── 4c. fgsea on modules ─────────────────────────────────────────────────
  
  if (length(module_ec_sets) > 0 && length(ec_ranked_vec) >= 10) {
    set.seed(42)
    fgsea_ec_res <- tryCatch(
      fgsea(
        pathways = module_ec_sets,
        stats    = ec_ranked_vec,
        minSize  = 3,
        maxSize  = 500,
        nperm    = 1000
      ) %>%
        as_tibble() %>%
        arrange(pval) %>%
        mutate(
          sig          = case_when(padj < 0.05 ~ "*", padj < 0.25 ~ "†", TRUE ~ ""),
          direction    = if_else(NES > 0, "Enriched in FMT", "Depleted in FMT"),
          module_label = str_replace(pathway, "^M", "M") %>%
            str_replace_all("_", " ")
        ),
      error = function(e) { message("EC fgsea (module) error: ", e$message); NULL }
    )
    
    if (!is.null(fgsea_ec_res)) {
      cat(sprintf("fgsea modules — padj<0.25: %d\n",
                  sum(fgsea_ec_res$padj < 0.25, na.rm = TRUE)))
      write.csv(fgsea_ec_res %>% dplyr::select(-leadingEdge),
                file.path(OUT_DIR, "fgsea_ec_kegg_modules.csv"),
                row.names = FALSE)
    }
  }
  
  # ── 4d. fgsea on pathways ────────────────────────────────────────────────
  
  if (length(pathway_ec_sets) > 0 && length(ec_ranked_vec) >= 10) {
    set.seed(42)
    fgsea_ec_pathway_res <- tryCatch(
      fgsea(
        pathways = pathway_ec_sets,
        stats    = ec_ranked_vec,
        minSize  = 3,
        maxSize  = 500,
        nperm    = 1000
      ) %>%
        as_tibble() %>%
        arrange(pval) %>%
        mutate(
          sig          = case_when(padj < 0.05 ~ "*", padj < 0.25 ~ "†", TRUE ~ ""),
          direction    = if_else(NES > 0, "Enriched in FMT", "Depleted in FMT"),
          pathway_label = str_replace(pathway, "^map", "")
        ),
      error = function(e) { message("EC fgsea (pathway) error: ", e$message); NULL }
    )
    
    if (!is.null(fgsea_ec_pathway_res)) {
      cat(sprintf("fgsea pathways — padj<0.25: %d\n",
                  sum(fgsea_ec_pathway_res$padj < 0.25, na.rm = TRUE)))
      write.csv(fgsea_ec_pathway_res %>% dplyr::select(-leadingEdge),
                file.path(OUT_DIR, "fgsea_ec_kegg_pathways.csv"),
                row.names = FALSE)
    }
  }
} else {
  cat("maas_ec_full not found, skipping EC fgsea\n")
}


# =============================================================================
#  PART 5 — 整合視覺化
# =============================================================================

cat("\n── Part 5: Combined visualisation ──────────────────────────────────────\n")

# ORA bubble plot
p_ora <- NULL
if (exists("ora_res") && nrow(ora_res) > 0) {
  top_ora <- ora_res %>%
    slice_min(padj, n = 20) %>%
    mutate(
      module_label = str_replace(module, "^M\\d+_", "") %>%
        str_replace_all("_", " "),
      sig = case_when(padj < 0.05 ~ "**", padj < 0.25 ~ "*", TRUE ~ "")
    )
  
  p_ora <- ggplot(top_ora,
                  aes(x = log2(odds_ratio + 0.01),
                      y = reorder(module_label, -log10(padj)),
                      size = n_overlap, colour = -log10(padj))) +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = sig), hjust = -0.5, size = 4, colour = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_gradient(low = "grey70", high = "#D6604D",
                          name = "-log₁₀(padj)") +
    scale_size_continuous(range = c(3, 9), name = "# overlap EC") +
    labs(title    = "KEGG module ORA: EC enriched in FMT",
         subtitle = "FMT vs Placebo @ 6M  |  Fisher's exact (BH corrected)",
         x = "log₂(odds ratio)", y = NULL) +
    theme_fmt
  
  save_enrich_fig(p_ora, "fig3_ora_kegg_module.pdf", w = 11, h = 7)
}

# EC fgsea bar
p_fgsea_ec <- NULL
if (!is.null(fgsea_ec_res) && nrow(fgsea_ec_res) > 0) { # fgsea_ec_pathway_res or fgsea_ec_res
  plot_data <- fgsea_ec_res %>% 
    filter(!is.na(NES)) %>% 
    group_by(direction) %>% 
    slice_max(order_by = abs(NES), n =20, with_ties = FALSE) %>% 
    ungroup()
  p_fgsea_ec <- ggplot(plot_data, # 同樣使用過濾後的資料
                                        aes(x = reorder(module_label, NES), y = NES, fill = direction)) + # module_label or pathway
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(sig, sprintf(" (n=%d)", size)),
                  hjust = if_else(NES > 0, -0.1, 1.1)), size = 3, colour = "grey30") +
    geom_hline(yintercept = 0, colour = "grey40") +
    coord_flip() +
    scale_fill_manual(values = c("Enriched in FMT" = "#D6604D",
                                 "Depleted in FMT" = "#2166AC")) +
    labs(title    = "fgsea: KEGG modules (EC-based, MaAsLin2)", # modules or pathways
         subtitle = "FMT vs Placebo @ 6M", x = NULL, y = "NES", fill = NULL) +
    theme_fmt + theme(legend.position = "bottom")
  
  save_enrich_fig(p_fgsea_ec, "fig4_fgsea_ec_kegg_module.pdf", w = 10, h = 7)
}

# Summary panel
p_summary <- p_category_lollipop + p_fgsea_maas +
  plot_annotation(
    title    = "Functional enrichment summary: FMT vs Placebo @ 6M",
    subtitle = "Left: MetaCyc category rollup (MaAsLin2)  |  Right: fgsea NES (MaAsLin2)",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )

save_enrich_fig(p_summary, "fig0_summary_panel.pdf", w = 18, h = 8)


# =============================================================================
#  Summary
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("✓ Functional enrichment analysis complete\n")
cat(strrep("=", 60), "\n\n")
cat("Output:\n")
cat("  fig0_summary_panel.pdf\n")
cat("  fig1a_category_lollipop.pdf   MetaCyc rollup (MaAsLin2 only)\n")
cat("  fig1b_category_heatmap.pdf    Top 3 pathways per category\n")
cat("  fig2a_fgsea_maaslin2.pdf      fgsea NES (MaAsLin2)\n")
cat("  fig2b_fgsea_concordance.pdf   NES concordance: MaAsLin2 vs DESeq2\n")
cat("  fig3_ora_kegg_module.pdf      KEGG module ORA bubble\n")
cat("  fig4_fgsea_ec_kegg.pdf        KEGG module fgsea (EC)\n")
cat("\n  CSVs in:", OUT_DIR, "\n")