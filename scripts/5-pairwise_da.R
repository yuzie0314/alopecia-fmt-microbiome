# =============================================================================
#  5-pairwise_da.R
#  Pairwise timepoint DA — Alopecia Oral FMT Study
#
#  3 transitions × 2 groups × 3 data types = 18 comparisons per method
#    Transitions : BL→2M | 2M→6M | BL→6M
#    Groups      : FMT | placebo
#    Data types  : taxonomy (genus) | metacyc | ec
#    Methods     : MaAsLin2 (primary, TSS) + DESeq2 (sensitivity)
#
#  前置條件：
#    - ps_filt, meta_raw, metacyc_raw, ec_raw, theme_fmt 已在 0_1-loadData.R 定義
#    - GM03 無 Baseline → BL 相關比較自動排除
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

library(DESeq2)
library(phyloseq)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(Maaslin2)

select <- dplyr::select

# =============================================================================
#  Configuration
# =============================================================================

TRANSITIONS <- list(
  BL_2M = list(
    tp       = c("Baseline", "2M post-FMT"),
    ref      = "Timepoint,Baseline",
    excl_gm03 = TRUE
  ),
  M2_6M = list(
    tp       = c("2M post-FMT", "6M post-FMT"),
    ref      = "Timepoint,2M post-FMT",
    excl_gm03 = FALSE
  ),
  BL_6M = list(
    tp       = c("Baseline", "6M post-FMT"),
    ref      = "Timepoint,Baseline",
    excl_gm03 = TRUE
  )
)
GROUPS     <- c("FMT", "placebo")
DATA_TYPES <- c("taxonomy", "metacyc", "ec")
OUT_ROOT   <- "pairwise_results"

for (dt in DATA_TYPES)
  for (tr in names(TRANSITIONS))
    dir.create(file.path(OUT_ROOT, dt, tr), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "combined"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
#  Helper functions
# =============================================================================

# Subset phyloseq by SampleID vector; update Timepoint factor levels
subset_ps_ordered <- function(ps_obj, sids, tp_levels) {
  ps_sub <- subset_samples(ps_obj, SampleID %in% sids)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  sample_data(ps_sub)$Timepoint <- factor(
    as.character(sample_data(ps_sub)$Timepoint), levels = tp_levels)
  ps_sub
}

# Build phyloseq from pathway matrix (rows=pathways, cols=SampleID)
make_pathway_ps <- function(mat, meta_sub, tp_levels) {
  sids   <- meta_sub$SampleID
  sids   <- intersect(sids, colnames(mat))
  mat_s  <- mat[, sids, drop = FALSE]
  meta_p <- meta_sub %>% filter(SampleID %in% sids) %>%
    remove_rownames() %>%
    column_to_rownames("SampleID")
  meta_p$Timepoint <- factor(as.character(meta_p$Timepoint), levels = tp_levels)
  phyloseq(
    otu_table(mat_s, taxa_are_rows = TRUE),
    sample_data(meta_p)
  )
}

# ── DESeq2 ────────────────────────────────────────────────────────────────────
run_deseq2_pair <- function(ps_obj, trans, grp, dtype) {
  # Aggregate to genus for taxonomy
  if (dtype == "taxonomy") ps_obj <- tax_glom(ps_obj, "Genus", NArm = FALSE)

  otu_table(ps_obj) <- otu_table(
    round(otu_table(ps_obj)), taxa_are_rows = taxa_are_rows(ps_obj))

  dds <- phyloseq_to_deseq2(ps_obj, ~ PatientID + Timepoint)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)

  coef_nm <- resultsNames(dds)[grepl("Timepoint", resultsNames(dds))][1]
  res <- tryCatch(
    lfcShrink(dds, coef = coef_nm, type = "ashr", quiet = TRUE),
    error = function(e) results(dds, name = coef_nm)
  )

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("feature_id") %>%
    filter(!is.na(pvalue))

  if (dtype == "taxonomy") {
    tax_df <- as.data.frame(as(tax_table(ps_obj), "matrix")) %>%
      tibble::rownames_to_column("feature_id")
    res_df <- left_join(res_df, tax_df, by = "feature_id")
  }

  res_df %>%
    mutate(
      sig = case_when(
        pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**",
        pvalue < 0.05  ~ "*",   pvalue < 0.1  ~ "†",
        TRUE ~ "ns"
      ),
      direction = case_when(
        pvalue < 0.05 & log2FoldChange > 0 ~ "Up",
        pvalue < 0.05 & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      method = "DESeq2", data_type = dtype,
      transition = trans, group = grp
    ) %>%
    arrange(pvalue)
}

# ── MaAsLin2 ─────────────────────────────────────────────────────────────────
run_maaslin2_pair <- function(ps_obj, meta_sub, trans, grp, ref_str,
                               dtype, out_dir) {
  if (dtype == "taxonomy") {
    ps_g <- transform_sample_counts(
      tax_glom(ps_obj, "Genus", NArm = FALSE),
      function(x) x / sum(x))
    mat  <- t(as.data.frame(otu_table(ps_g)))
    norm <- "NONE"
  } else {
    mat  <- t(as.data.frame(otu_table(ps_obj)))
    norm <- "TSS"
  }

  meta_m <- meta_sub %>%
    dplyr::select(SampleID, PatientID, Timepoint, Age) %>%
    remove_rownames() %>%
    column_to_rownames("SampleID")

  fit <- tryCatch(
    Maaslin2(
      input_data       = mat,
      input_metadata   = meta_m,
      output           = out_dir,
      fixed_effects    = c("Timepoint", "Age"),
      random_effects   = c("PatientID"),
      reference        = ref_str,
      normalization    = norm,
      transform        = "LOG",
      analysis_method  = "LM",
      min_prevalence   = 0.1,
      min_abundance    = 0.0,
      max_significance = 0.25,
      plot_heatmap     = FALSE,
      plot_scatter     = FALSE
    ),
    error = function(e) { message("MaAsLin2 error [", trans, "|", grp, "|", dtype, "]: ", e$message); NULL }
  )

  if (is.null(fit)) return(NULL)

  fit$results %>%
    as_tibble() %>%
    rename(feature_name = feature) %>%
    filter(metadata == "Timepoint") %>%
    mutate(method = "MaAsLin2", data_type = dtype,
           transition = trans, group = grp) %>%
    arrange(qval)
}

# =============================================================================
#  Main loop
# =============================================================================

all_d2   <- list()
all_maas <- list()

for (trans in names(TRANSITIONS)) {
  cfg     <- TRANSITIONS[[trans]]
  tp_lvls <- cfg$tp
  ref_str <- cfg$ref

  for (grp in GROUPS) {
    cat(sprintf("\n── %s | %s ─────────────────────────────────────────────\n", trans, grp))

    meta_sub <- meta_raw %>%
      filter(Treatment == grp, Timepoint %in% tp_lvls) %>%
      mutate(Timepoint = factor(Timepoint, levels = tp_lvls))
    if (cfg$excl_gm03) meta_sub <- meta_sub %>% filter(PatientID != "GM03")

    sids <- meta_sub$SampleID

    for (dtype in DATA_TYPES) {
      key <- paste(trans, grp, dtype, sep = "_")

      # Build phyloseq
      ps_sub <- tryCatch({
        if (dtype == "taxonomy") {
          subset_ps_ordered(ps_filt, sids, tp_lvls)
        } else {
          raw_mat <- if (dtype == "metacyc") metacyc_raw else ec_raw
          make_pathway_ps(raw_mat, meta_sub, tp_lvls)
        }
      }, error = function(e) { message("Build PS error: ", e$message); NULL })

      if (is.null(ps_sub)) next

      # ── DESeq2
      d2 <- tryCatch(
        run_deseq2_pair(ps_sub, trans, grp, dtype),
        error = function(e) { message("DESeq2 error [", key, "]: ", e$message); NULL }
      )
      if (!is.null(d2)) {
        fpath <- file.path(OUT_ROOT, dtype, trans,
                           sprintf("deseq2_%s_%s.csv", grp, trans))
        write.csv(d2, fpath, row.names = FALSE)
        all_d2[[key]] <- d2
        cat(sprintf("  DESeq2  %s: %d sig (p<0.05)\n", dtype,
                    sum(d2$pvalue < 0.05, na.rm = TRUE)))
      }

      # ── MaAsLin2
      out_dir_m <- file.path(OUT_ROOT, dtype, trans,
                             sprintf("maaslin2_%s_%s", grp, trans))
      m2 <- run_maaslin2_pair(ps_sub, meta_sub, trans, grp, ref_str, dtype, out_dir_m)
      if (!is.null(m2)) {
        write.csv(m2, paste0(out_dir_m, ".csv"), row.names = FALSE)
        all_maas[[key]] <- m2
        cat(sprintf("  MaAsLin2 %s: %d sig (p<0.05)\n", dtype,
                    sum(m2$pval < 0.05, na.rm = TRUE)))
      }
    }
  }
}

# =============================================================================
#  Combine & summarise
# =============================================================================

d2_combined   <- bind_rows(all_d2)
maas_combined <- bind_rows(all_maas)

write.csv(d2_combined,   file.path(OUT_ROOT, "combined", "deseq2_all_pairwise.csv"),   row.names = FALSE)
write.csv(maas_combined, file.path(OUT_ROOT, "combined", "maaslin2_all_pairwise.csv"), row.names = FALSE)

summary_counts <- maas_combined %>%
  group_by(data_type, transition, group) %>%
  summarise(
    total    = n(),
    sig_0.05 = sum(pval < 0.05,  na.rm = TRUE),
    sig_0.25 = sum(qval < 0.25,  na.rm = TRUE),
    .groups  = "drop"
  )
write.csv(summary_counts, file.path(OUT_ROOT, "combined", "summary_sig_counts.csv"), row.names = FALSE)

# =============================================================================
#  Visualisation
# =============================================================================

save_fig <- function(p, path, w = 12, h = 7) {
  ggsave(path, plot = p, width = w, height = h, dpi = 300, bg = "white")
}

# ── LFC heatmap (MaAsLin2 coef across 3 transitions)
make_lfc_heatmap <- function(df, grp_filter, dtype_filter, n_top = 30, title = "") {
  top_feats <- df %>%
    filter(group == grp_filter, data_type == dtype_filter, pval < 0.05) %>%
    count(feature_name) %>%
    slice_max(n, n = n_top) %>%
    pull(feature_name)

  if (length(top_feats) == 0) return(NULL)

  df %>%
    filter(group == grp_filter, data_type == dtype_filter,
           feature_name %in% top_feats) %>%
    mutate(transition = factor(transition, levels = c("BL_2M", "M2_6M", "BL_6M")),
           star = ifelse(pval < 0.05, "*", "")) %>%
    ggplot(aes(x = transition, y = feature_name, fill = coef)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = star), size = 3.5, vjust = 0.75) +
    scale_fill_gradient2(low = "#D85A30", mid = "white", high = "#1D9E75",
                         midpoint = 0, name = "coef") +
    labs(title = title, x = "Transition", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 7),
          plot.title  = element_text(size = 11, face = "bold"))
}

# ── Volcano panel (3 transitions)
make_volcano_panel <- function(df, grp_filter, dtype_filter) {
  df_sub <- df %>%
    filter(group == grp_filter, data_type == dtype_filter) %>%
    mutate(
      direction  = case_when(
        pval < 0.05 & coef > 0 ~ "Up",
        pval < 0.05 & coef < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      transition = factor(transition, levels = c("BL_2M", "M2_6M", "BL_6M"))
    )

  lbl_df <- df_sub %>% filter(pval < 0.05) %>%
    group_by(transition) %>% slice_max(abs(coef), n = 10)

  ggplot(df_sub, aes(x = coef, y = -log10(pval), colour = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = lbl_df, aes(label = feature_name),
                    size = 2.4, max.overlaps = 15, segment.colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    scale_colour_manual(
      values = c(Up = "#1D9E75", Down = "#D85A30", NS = "grey70"), guide = "none") +
    facet_wrap(~ transition, scales = "free_x") +
    labs(x = "Coefficient", y = "-log10(p)") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey95"))
}

for (dtype in DATA_TYPES) {
  for (grp in GROUPS) {
    tag <- sprintf("%s_%s", dtype, grp)

    p_heat <- make_lfc_heatmap(
      maas_combined, grp, dtype,
      title = sprintf("MaAsLin2 LFC — %s | %s", dtype, grp))
    if (!is.null(p_heat))
      save_fig(p_heat,
               file.path(OUT_ROOT, "combined", sprintf("heatmap_maaslin2_%s.pdf", tag)),
               w = 9, h = 10)

    p_volc <- make_volcano_panel(maas_combined, grp, dtype)
    if (!is.null(p_volc))
      save_fig(p_volc,
               file.path(OUT_ROOT, "combined", sprintf("volcano_panel_%s.pdf", tag)),
               w = 12, h = 5)
  }
}

cat("\n✓ pairwise_timepoint_da complete\n")
cat(sprintf("  DESeq2 rows   : %d\n", nrow(d2_combined)))
cat(sprintf("  MaAsLin2 rows : %d\n", nrow(maas_combined)))
print(summary_counts)
