# =============================================================================
#  ── 修改 Section E：加入 |LFC| > 1 門檻 ──────────────────────────────────
#  把原本的 fmt_specific 定義換成以下版本
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

fmt_specific <- res_A %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > 1) %>%   # 加 |LFC|>1
  semi_join(
    res_B %>% filter(pvalue < 0.05, abs(log2FoldChange) > 1),
    by = "Genus"
  ) %>%
  anti_join(
    res_C %>% filter(pvalue < 0.05),
    by = "Genus"
  )

cat(sprintf("FMT-specific genera (p<0.05 & |LFC|>1): %d\n", nrow(fmt_specific)))
print(fmt_specific %>% select(Genus, log2FoldChange, pvalue, sig))
write.csv(fmt_specific, "deseq2_E_FMT_specific_taxa.csv", row.names = FALSE)


# =============================================================================
#  SECTION F — LEfSe on Δabundance (6M − Baseline)
#
#  核心邏輯：
#    每個病人只有一個 Δ 值 → 獨立觀測 → LEfSe 可正確比較 FMT vs placebo
#    Δ = log2(rel.abund 6M + ε) − log2(rel.abund BL + ε)
#      = log2 ratio，正值代表 6M 後增加
# =============================================================================

cat("\n── Section F: LEfSe on Δabundance ─────────────────────────────────────\n")

if (!requireNamespace("microbiomeMarker", quietly = TRUE))
  BiocManager::install("microbiomeMarker")
library(microbiomeMarker)

# ── F0. 準備相對豐度（genus level）────────────────────────────────────────────

rel_genus <- tax_glom(ps_filt, taxrank = "Genus") %>%
  transform_sample_counts(function(x) x / sum(x))

meta_df <- as.data.frame(sample_data(rel_genus))

# 有完整 BL + 6M 的病人（排除 GM03）
patients_complete <- meta_df %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT"),
         PatientID != "GM03") %>%
  group_by(PatientID) %>%
  filter(n() == 2) %>%
  pull(PatientID) %>%
  unique()

cat(sprintf("Patients with complete BL + 6M: %d\n", length(patients_complete)))

# SampleID 按 PatientID 排序，確保 BL 和 6M 一一對應
sid_bl <- meta_df %>%
  filter(Timepoint == "Baseline", PatientID %in% patients_complete) %>%
  arrange(PatientID) %>% pull(SampleID)

sid_6m <- meta_df %>%
  filter(Timepoint == "6M post-FMT", PatientID %in% patients_complete) %>%
  arrange(PatientID) %>% pull(SampleID)

# ── F1. Δ matrix（log2 ratio）────────────────────────────────────────────────

abund_mat <- as.matrix(otu_table(rel_genus))   # rows = taxa, cols = samples

delta_mat <- log2(abund_mat[, sid_6m] + 1e-6) -
             log2(abund_mat[, sid_bl] + 1e-6)

colnames(delta_mat) <- patients_complete   # cols = PatientID

cat(sprintf("Delta matrix: %d genera × %d patients\n",
            nrow(delta_mat), ncol(delta_mat)))

# ── F2. 建立 phyloseq for LEfSe ──────────────────────────────────────────────

meta_delta <- meta_df %>%
  filter(Timepoint == "Baseline", PatientID %in% patients_complete) %>%
  arrange(PatientID) %>%
  select(PatientID, Treatment, Age) %>%
  column_to_rownames("PatientID")

tax_delta <- tax_table(rel_genus)[rownames(delta_mat), ]

ps_delta <- phyloseq(
  otu_table(delta_mat, taxa_are_rows = TRUE),
  tax_delta,
  sample_data(meta_delta)
)

# ── F3. LEfSe ────────────────────────────────────────────────────────────────

set.seed(42)
lefse_res <- run_lefse(
  ps_delta,
  group           = "Treatment",
  taxa_rank       = "none",    # 已是 genus level
  wilcoxon_cutoff = 0.05,
  kw_cutoff       = 0.05,
  lda_cutoff      = 2.0,
  norm            = "none"     # delta 不需要 normalize
)

lefse_tab <- marker_table(lefse_res) %>%
  as.data.frame() %>%
  rename(Genus = feature) %>%
  arrange(desc(abs(ef_lda)))

cat(sprintf("LEfSe significant genera (LDA>2): %d\n", nrow(lefse_tab)))
print(lefse_tab)

write.csv(lefse_tab, "lefse_delta_abundance.csv", row.names = FALSE)

# ── F4. LEfSe LDA bar plot ────────────────────────────────────────────────────

p_lefse <- plot_ef_bar(lefse_res) +
  labs(
    title    = "LEfSe on Δabundance (6M − Baseline)",
    subtitle = "log₂(rel.abund 6M) − log₂(rel.abund BL)  |  LDA>2, p<0.05",
    x        = "LDA score (log₁₀)"
  ) +
  theme_fmt +
  theme(axis.text.y = element_text(face = "italic", size = 9))

# ── F5. DESeq2 vs LEfSe concordance ──────────────────────────────────────────

cat("\n── DESeq2 (Section A) vs LEfSe concordance ─────────────────────────────\n")

concordance <- lefse_tab %>%
  left_join(
    res_A %>% select(Genus, log2FoldChange, pvalue, sig),
    by = "Genus"
  ) %>%
  mutate(
    lefse_dir  = if_else(enrich_group == "FMT", "Up", "Down"),
    deseq2_dir = case_when(
      is.na(log2FoldChange) ~ "Not tested",
      log2FoldChange > 0    ~ "Up",
      TRUE                  ~ "Down"
    ),
    concordant = lefse_dir == deseq2_dir
  ) %>%
  select(Genus, enrich_group, ef_lda, lefse_dir,
         log2FoldChange, pvalue, deseq2_dir, concordant)

cat("Direction concordance:\n")
print(concordance)

# 三方一致：LEfSe + DESeq2_A + DESeq2_B
triple_sig <- concordance %>%
  filter(concordant == TRUE, pvalue < 0.05) %>%
  semi_join(
    res_B %>% filter(pvalue < 0.05, abs(log2FoldChange) > 1),
    by = "Genus"
  )

cat(sprintf("\nTriple concordant (LEfSe + DESeq2_A + DESeq2_B): %d\n",
            nrow(triple_sig)))
print(triple_sig %>% select(Genus, enrich_group, ef_lda, log2FoldChange, pvalue))

write.csv(concordance, "lefse_deseq2_concordance.csv",  row.names = FALSE)
write.csv(triple_sig,  "lefse_deseq2_triple_sig.csv",   row.names = FALSE)

# ── F6. Concordance scatter：DESeq2 LFC vs LEfSe LDA ─────────────────────────

if (nrow(concordance %>% filter(!is.na(log2FoldChange))) > 0) {
  p_concordance <- ggplot(
    concordance %>% filter(!is.na(log2FoldChange)),
    aes(x = log2FoldChange, y = ef_lda,
        colour = concordant, shape = concordant)
  ) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60") +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(aes(label = Genus), size = 2.8,
                    fontface = "italic", segment.colour = "grey60",
                    max.overlaps = 15) +
    scale_colour_manual(
      values = c(`TRUE`  = "#1D9E75", `FALSE` = "#D85A30"),
      labels = c(`TRUE`  = "Concordant", `FALSE` = "Discordant")
    ) +
    scale_shape_manual(
      values = c(`TRUE`  = 16, `FALSE` = 17),
      labels = c(`TRUE`  = "Concordant", `FALSE` = "Discordant")
    ) +
    labs(
      title    = "DESeq2 LFC vs LEfSe LDA score",
      subtitle = "Dashed lines: |LFC|=1, |LDA|=2  |  Concordant = same direction",
      x        = "DESeq2 log₂FC (FMT vs placebo @ 6M)",
      y        = "LEfSe LDA score (Δabundance)",
      colour   = NULL, shape = NULL
    ) +
    theme_fmt
}

# ── F7. 新增 save ─────────────────────────────────────────────────────────────
# 加到你原本的 save_fig block 裡

save_fig(p_lefse,       "deseq2_fig6_lefse_delta.pdf",              w = 10, h = 8)
if (exists("p_concordance"))
  save_fig(p_concordance, "deseq2_fig7_deseq2_lefse_concordance.pdf", w = 9,  h = 7)

cat("\n── Summary ─────────────────────────────────────────────────────────────\n")
cat(sprintf("  DESeq2 FMT-specific taxa (p<0.05 & |LFC|>1): %d\n", nrow(fmt_specific)))
cat(sprintf("  LEfSe significant genera (LDA>2):             %d\n", nrow(lefse_tab)))
cat(sprintf("  Triple concordant:                            %d\n", nrow(triple_sig)))
cat("\n  triple_sig = most reliable FMT-driven taxa\n")
