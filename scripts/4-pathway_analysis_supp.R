# =============================================================================
#  SECTION 2b — Differential pathway（DESeq2 sensitivity analysis）
#  注意：ec 值 round 成整數，統計假設不完全成立
#  作為 MaAsLin2 結果的 sensitivity check
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

cat("\n── Section 2b: DESeq2 pathway (sensitivity analysis) ───────────────────\n")

make_pw_phyloseq <- function(mat, meta, timepoint_filter = NULL,
                             treatment_filter = NULL,
                             exclude_patient  = NULL) {
  # 篩選樣本
  meta_sub <- meta
  if (!is.null(timepoint_filter))
    meta_sub <- meta_sub %>% filter(Timepoint %in% timepoint_filter)
  if (!is.null(treatment_filter))
    meta_sub <- meta_sub %>% filter(Treatment %in% treatment_filter)
  if (!is.null(exclude_patient))
    meta_sub <- meta_sub %>% filter(!PatientID %in% exclude_patient)
  
  common   <- intersect(colnames(mat), meta_sub$SampleID)
  meta_sub <- meta_sub %>%
    filter(SampleID %in% common) %>%
    remove_rownames() %>% 
    tibble::column_to_rownames("SampleID")
  
  # round → integer count matrix
  count_mat <- round(mat[, common])
  mode(count_mat) <- "integer"
  
  # 移除全零 pathway
  count_mat <- count_mat[rowSums(count_mat) > 0, ]
  
  # 建立簡易 phyloseq（只需要 otu_table + sample_data）
  otu <- otu_table(count_mat, taxa_are_rows = TRUE)
  smd <- sample_data(meta_sub)
  phyloseq(otu, smd)
}

tidy_pw_deseq2 <- function(dds, contrast, alpha = 0.05) {
  res <- results(dds, contrast = contrast, alpha = alpha,
                 independentFiltering = TRUE)
  
  res_shrink <- tryCatch(
    lfcShrink(dds, contrast = contrast, res = res,
              type = "ashr", quiet = TRUE),
    error = function(e) { message("lfcShrink skipped: ", e$message); res }
  )
  
  as.data.frame(res_shrink) %>%
    tibble::rownames_to_column("Pathway") %>%
    filter(!is.na(pvalue)) %>%
    mutate(
      sig = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01  ~ "**",
        pvalue < 0.05  ~ "*",
        pvalue < 0.1   ~ "†",
        TRUE           ~ "ns"
      ),
      direction = case_when(
        pvalue < 0.05 & log2FoldChange > 0 ~ "Up",
        pvalue < 0.05 & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(pvalue)
}

# ── A. FMT vs Placebo @ 6M ───────────────────────────────────────────────────

ps_pw_6m <- make_pw_phyloseq(ec, meta_raw,
                             timepoint_filter = "6M post-FMT")
sample_data(ps_pw_6m)$Treatment <- relevel(
  factor(sample_data(ps_pw_6m)$Treatment), ref = "placebo"
)

dds_pw_6m <- phyloseq_to_deseq2(ps_pw_6m, ~ Treatment)
dds_pw_6m <- estimateSizeFactors(dds_pw_6m, type = "poscounts")
dds_pw_6m <- DESeq(dds_pw_6m, test = "Wald",
                   fitType = "parametric", quiet = TRUE)

res_pw_A <- tidy_pw_deseq2(dds_pw_6m,
                           contrast = c("Treatment", "FMT", "placebo"))

cat(sprintf("DESeq2 [A] FMT vs placebo @ 6M — sig pathways (p<0.05): %d\n",
            sum(res_pw_A$pvalue < 0.05, na.rm = TRUE)))

# ── B. FMT 組內：Baseline vs 6M（paired）─────────────────────────────────────

ps_pw_fmt <- make_pw_phyloseq(ec,
                              meta_raw,
                              timepoint_filter = c("Baseline", "6M post-FMT"),
                              treatment_filter = "FMT",
                              exclude_patient  = "GM03")
sample_data(ps_pw_fmt)$Timepoint <- relevel(
  factor(sample_data(ps_pw_fmt)$Timepoint,
         levels = c("Baseline", "6M post-FMT")),
  ref = "Baseline"
)

dds_pw_fmt <- phyloseq_to_deseq2(ps_pw_fmt, ~ PatientID + Timepoint)
dds_pw_fmt <- estimateSizeFactors(dds_pw_fmt, type = "poscounts")
dds_pw_fmt <- DESeq(dds_pw_fmt, test = "Wald",
                    fitType = "parametric", quiet = TRUE)

res_pw_B <- tidy_pw_deseq2(dds_pw_fmt,
                           contrast = c("Timepoint", "6M post-FMT", "Baseline"))

cat(sprintf("DESeq2 [B] FMT Baseline→6M — sig pathways (p<0.05): %d\n",
            sum(res_pw_B$pvalue < 0.05, na.rm = TRUE)))

# ── C. Placebo 組內（對照）───────────────────────────────────────────────────

ps_pw_plc <- make_pw_phyloseq(ec,meta_raw,
                              timepoint_filter = c("Baseline", "6M post-FMT"),
                              treatment_filter = "placebo")
sample_data(ps_pw_plc)$Timepoint <- relevel(
  factor(sample_data(ps_pw_plc)$Timepoint,
         levels = c("Baseline", "6M post-FMT")),
  ref = "Baseline"
)

dds_pw_plc <- phyloseq_to_deseq2(ps_pw_plc, ~ PatientID + Timepoint)
dds_pw_plc <- estimateSizeFactors(dds_pw_plc, type = "poscounts")
dds_pw_plc <- DESeq(dds_pw_plc, test = "Wald",
                    fitType = "parametric", quiet = TRUE)

res_pw_C <- tidy_pw_deseq2(dds_pw_plc,
                           contrast = c("Timepoint", "6M post-FMT", "Baseline"))

cat(sprintf("DESeq2 [C] Placebo Baseline→6M — sig pathways (p<0.05): %d\n",
            sum(res_pw_C$pvalue < 0.05, na.rm = TRUE)))

# ── D. FMT-specific（A ∩ B, not C）──────────────────────────────────────────
# plot_volcano in 3-da.R at row 86
p_volcano_A_ec <- plot_volcano(res_df = res_pw_A, lab = "Pathway", title = "ec DESeq2 [A] FMT vs Placebo @ 6M")
p_volcano_B_ec <- plot_volcano(res_df = res_pw_B, lab = "Pathway", title = "ec DESeq2 [B] FMT Baseline→6M")
p_volcano_C_ec <- plot_volcano(res_df = res_pw_C, lab = "Pathway", title = "ec DESeq2 [C] Placebo Baseline→6M")

fmt_pw_deseq2 <- res_pw_A %>%
  filter(pvalue < 0.05, log2FoldChange > 0) %>%
  semi_join(res_pw_B %>% filter(pvalue < 0.05, log2FoldChange > 0),
            by = "Pathway") %>%
  anti_join(res_pw_C %>% filter(pvalue < 0.05),
            by = "Pathway")

cat(sprintf("DESeq2 FMT-specific pathways (A∩B not C): %d\n",
            nrow(fmt_pw_deseq2)))
print(fmt_pw_deseq2 %>% dplyr::select(Pathway, log2FoldChange, pvalue, sig))

# ── E. MaAsLin2 vs DESeq2 一致性比較 ────────────────────────────────────────

overlap_pw <- intersect(fmt_specific_pw, fmt_pw_deseq2$Pathway)

cat(sprintf("\nPathways significant in BOTH MaAsLin2 and DESeq2: %d\n",
            length(overlap_pw)))
print(overlap_pw)

# 輸出
write.csv(res_pw_A,       "outputs/pathway/deseq2_ec_pw_A_FMT_vs_placebo_6M.csv",    row.names = FALSE)
write.csv(res_pw_B,       "outputs/pathway/deseq2_ec_pw_B_FMT_baseline_6M.csv",      row.names = FALSE)
write.csv(res_pw_C,       "outputs/pathway/deseq2_ec_pw_C_placebo_baseline_6M.csv",  row.names = FALSE)
write.csv(fmt_pw_deseq2,  "outputs/pathway/deseq2_ec_pw_D_FMT_specific.csv",         row.names = FALSE)
write.csv(data.frame(Pathway = overlap_pw),
          "outputs/pathway/ec_MaAsLin2_DESeq2_overlap.csv", row.names = FALSE)

save_fig(p_volcano_A_ec,     "outputs/pathway/deseq2_ec_fig1a_volcano_FMT_vs_placebo_6M.pdf")
save_fig(p_volcano_B_ec,     "outputs/pathway/deseq2_ec_fig1b_volcano_FMT_baseline_6M.pdf")
save_fig(p_volcano_C_ec,     "outputs/pathway/deseq2_ec_fig1c_volcano_placebo_baseline_6M.pdf")

cat("\nNote: DESeq2 on ec pathway data uses rounded integer values.\n")
cat("Treat as sensitivity analysis; primary results from MaAsLin2.\n")