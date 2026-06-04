# =============================================================================
#  deseq2_analysis.R
#  Differential abundance analysis — Alopecia Oral FMT Study
#
#  比較設計：
#    A. FMT vs Placebo @ 6M（主要比較）
#    B. FMT 組內：Baseline vs 6M（within-group temporal change, paired）
#    C. Placebo 組內：Baseline vs 6M（對照確認）
#    D. 三時間點整體趨勢（LRT）
#    E. FMT-specific taxa（A ∩ B, not C）+ |LFC| > 1
#    F. LEfSe on Δabundance（paired 結構在資料層面處理）
#    G. 視覺化
#
#  前置條件：
#    - ps_filt, meta_raw, col_treatment, col_timepoint, tp_labels, theme_fmt
#      已在主 script 定義
#    - GM03 無 Baseline，相關比較自動排除
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

library(DESeq2)
library(phyloseq)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(Maaslin2)


# =============================================================================
#  輔助 functions
# =============================================================================

# ── phyloseq → DESeq2 dataset ─────────────────────────────────────────────────
make_dds <- function(ps_obj, formula) {
  otu_table(ps_obj) <- otu_table(
    round(otu_table(ps_obj)),
    taxa_are_rows = taxa_are_rows(ps_obj)
  )
  dds <- phyloseq_to_deseq2(ps_obj, formula)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  return(dds)
}

# ── DESeq2 結果整理 ────────────────────────────────────────────────────────────
# pvalue < 0.05  & log2FoldChange > 0 ~ "Up"
# pvalue < 0.05  & log2FoldChange < 0 ~ "Down"
tidy_deseq2 <- function(dds, ps_obj, contrast, alpha = 0.05, lfc_threshold = 0) {
  res <- results(dds, contrast = contrast, alpha = alpha,
                 lfcThreshold = lfc_threshold, independentFiltering = TRUE)
  
  res_shrink <- tryCatch(
    lfcShrink(dds, contrast = contrast, res = res, type = "ashr", quiet = TRUE),
    error = function(e) {
      message("lfcShrink failed, using unshrunken LFC: ", e$message)
      res
    }
  )
  
  tax_df <- as.data.frame(as(tax_table(ps_obj), "matrix")) %>%
    tibble::rownames_to_column("ASV_ID")
  
  as.data.frame(res_shrink) %>%
    tibble::rownames_to_column("ASV_ID") %>%
    left_join(tax_df, by = "ASV_ID") %>%
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
        pvalue < 0.05  & log2FoldChange > 0 ~ "Up", # & abs(log2FoldChange) > 1
        pvalue < 0.05  & log2FoldChange < 0 ~ "Down", # & abs(log2FoldChange) > 1
        TRUE ~ "NS"
      )
    ) %>%
    arrange(pvalue)
}

# ── Volcano plot ──────────────────────────────────────────────────────────────
# pval <0.05 & > lfc_cut label text 
# pval <0.05 show up and down
plot_volcano <- function(res_df, title,lab = "Genus",
                         lfc_cut = 1, pvalue_cut = 0.05, n_label = 15) {
  label_df <- res_df %>%
    filter(pvalue < pvalue_cut, abs(log2FoldChange) > lfc_cut) %>%
    slice_max(abs(log2FoldChange), n = n_label)
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue),
                     colour = direction)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_text_repel(data = label_df, aes(label = !!sym(lab)),
                    size = 2.8, fontface = "italic",
                    max.overlaps = 20, segment.colour = "grey60",
                    box.padding  = 0.4) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = -log10(pvalue_cut),
               linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    scale_colour_manual(
      values = c(Up = "#D6604D", Down = "#2166AC", NS = "grey70"),
      labels = c(Up = "Enriched", Down = "Depleted", NS = "NS")
    ) +
    labs(title    = title,
         subtitle = sprintf("Up: %d  |  Down: %d  |  p<%.2f & |LFC|>%d",
                            sum(res_df$direction == "Up"),
                            sum(res_df$direction == "Down"),
                            pvalue_cut, lfc_cut),
         x = "Log2 fold change", y = "-log10(p-value)", colour = NULL) +
    theme_fmt
}

# ── LFC bar chart ─────────────────────────────────────────────────────────────
plot_lfc_bar <- function(res_df, title, n = 20, lfc_cut = 1) {
  top <- res_df %>%
    filter(pvalue < 0.05, abs(log2FoldChange) > lfc_cut) %>%
    slice_max(abs(log2FoldChange), n = n) %>%
    mutate(label = coalesce(Genus, paste0("ASV_", str_sub(ASV_ID, 1, 6))),
           label = if_else(is.na(label) | label == "NA", ASV_ID, label))
  
  if (nrow(top) == 0) {
    message(sprintf("[%s] No significant taxa (p<0.05 & |LFC|>0)", title))
    return(NULL)
  }
  
  ggplot(top, aes(x = reorder(label, log2FoldChange),
                  y = log2FoldChange, fill = direction)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                      ymax = log2FoldChange + lfcSE),
                  width = 0.3, colour = "grey40") +
    geom_text(aes(label = sig),
              hjust = if_else(top$log2FoldChange > 0, -0.3, 1.3),
              size  = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c(Up = "#D6604D", Down = "#2166AC")) +
    labs(title = title, subtitle = "Error bar = ±SE (shrunken LFC)",
         x = NULL, y = "Log2 fold change", fill = NULL) +
    theme_fmt +
    theme(axis.text.y = element_text(face = "italic", size = 9))
}

# ── MA plot ───────────────────────────────────────────────────────────────────
# pval <0.05 & > lfc_cut label text 
# pval <0.05 show up and down
plot_ma <- function(res_df, title, lfc_cut = 1) {
  ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange,
                     colour = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text_repel(data = filter(res_df, pvalue < 0.05, abs(log2FoldChange) > lfc_cut),
                    aes(label = Genus), size = 2.5, fontface = "italic",
                    max.overlaps = 15, segment.colour = "grey60") +
    scale_colour_manual(values = c(Up = "#D6604D", Down = "#2166AC", NS = "grey75")) +
    labs(title = title, x = "log10(mean abundance + 1)",
         y = "Log2 fold change", colour = NULL) +
    theme_fmt
}


# =============================================================================
#  SECTION A — FMT vs Placebo @ 6M
# =============================================================================

cat("\n── Section A: FMT vs Placebo @ 6M ─────────────────────────────────────\n")

ps_6m       <- subset_samples(ps_filt, Timepoint == "6M post-FMT")
ps_6m_genus <- tax_glom(ps_6m, taxrank = "Genus", NArm = FALSE) # 144 tax
ps_6m_genus <- filter_ps_by_prev(ps_obj = ps_6m_genus) # 139 tax

cat(sprintf("Samples: %d  |  Genera: %d\n",
            nsamples(ps_6m_genus), ntaxa(ps_6m_genus)))

sample_data(ps_6m_genus)$Treatment <- relevel(
  factor(sample_data(ps_6m_genus)$Treatment), ref = "placebo" ### REF
)

dds_6m <- make_dds(ps_6m_genus, ~ Treatment)
min_samples <- ceiling(0.10 * ncol(dds_6m))
keep_sub <- rowSums(counts(dds_6m) >= 2) >= min_samples
dds_6m <- dds_6m[keep_sub, ]
# plotDispEsts(dds_6m)

dds_6m <- DESeq(dds_6m, test = "Wald", fitType = "parametric", quiet = TRUE)
plotDispEsts(dds_6m)

counts_data <- plotCounts(dds = dds_6m, gene = "ASV0295", intgroup = "Treatment", returnData = TRUE) # check Klebsiella


ggplot(data = counts_data, aes(x = Treatment, y = count, color = Treatment))+
  geom_point(size = 6, alpha = 0.6)+
  scale_color_manual(values = col_treatment)+
  labs(y="Normalized Count", x = "")+
  theme_fmt

res_A <- tidy_deseq2(dds_6m, ps_6m_genus,
                     contrast = c("Treatment", "FMT", "placebo")) # (FMT 減去 placebo) 的配對差異

cat(sprintf("Significant (p<0.05 & |LFC|>0): Up=%d  Down=%d\n",
            sum(res_A$direction == "Up"), sum(res_A$direction == "Down")))
print(filter(res_A, pvalue < 0.05) %>%
        dplyr::select(Genus, log2FoldChange, lfcSE, pvalue, sig))

p_volcano_A <- plot_volcano(res_A, "FMT vs Placebo @ 6M post-FMT")
p_lfc_A     <- plot_lfc_bar(res_A, "Top diff. genera: FMT vs Placebo @ 6M")
p_ma_A      <- plot_ma(res_A,      "MA plot: FMT vs Placebo @ 6M")

write.csv(res_A, "outputs/tax_da/deseq2_A_FMT_vs_placebo_6M.csv", row.names = FALSE)

# check expected samples size in the future [FAILED]
# 1. 提取所有菌屬的平均讀數 (baseMean)
current_means <- mcols(dds_6m)$baseMean

# 2. 提取 DESeq2 計算出的離散度 (Dispersion)
current_disps <- dispersions(dds_6m)
valid_idx <- !is.na(current_disps) & current_disps > 0
biomarker_means <- current_means[valid_idx]
biomarker_disps <- current_disps[valid_idx]

# install.packages("ssizeRNA")
library(ssizeRNA)

size_summary <- ssizeRNA_single(
  m = length(biomarker_disps),   # 總菌屬數量 (Features)
  mu = median(biomarker_means),  # 目前菌屬豐度的中位數
  disp = median(biomarker_disps),# 目前菌屬變異度(離散度)的中位數
  fc = 2,                        # 你未來預期想抓到的最低Fold Change (2代表兩倍差異，即LFC=1)
  fdr = 0.05,                    # 預期控制的 BH 校正後 FDR 門檻
  power = 0.8,                   # 期望達到的統計檢定力 (80%)
  pi0 = 0.8                      # 預期有多少比例的菌「沒有差異」（0.8代表猜測有20%的菌有差）
)

# 查看建議的每組樣本數
print(size_summary$n) 



# =============================================================================
#  SECTION B — FMT 組內：Baseline vs 6M（paired）
# =============================================================================

cat("\n── Section B: FMT group — Baseline vs 6M ──────────────────────────────\n")

ps_fmt <- subset_samples(ps_filt,
                         Treatment == "FMT" &
                           Timepoint %in% c("Baseline", "6M post-FMT") &
                           PatientID != "GM03")

ps_fmt_genus <- tax_glom(ps_fmt, taxrank = "Genus", NArm = FALSE)
ps_fmt_genus <- filter_ps_by_prev(ps_obj = ps_fmt_genus)

cat(sprintf("Samples: %d  |  Genera: %d\n",
            nsamples(ps_fmt_genus), ntaxa(ps_fmt_genus)))

sample_data(ps_fmt_genus)$Timepoint <- relevel(
  factor(sample_data(ps_fmt_genus)$Timepoint,
         levels = c("Baseline", "6M post-FMT")),
  ref = "Baseline"
)

dds_fmt <- make_dds(ps_fmt_genus, ~ PatientID + Timepoint)
dds_fmt <- DESeq(dds_fmt, test = "Wald", fitType = "parametric", quiet = TRUE)

res_B <- tidy_deseq2(dds_fmt, ps_fmt_genus,
                     contrast = c("Timepoint", "6M post-FMT", "Baseline"))

cat(sprintf("Significant (p<0.05 & |LFC|>0): %d\n",
            sum(res_B$direction != "NS")))
print(filter(res_B, pvalue < 0.05) %>%
        dplyr::select(Genus, log2FoldChange, lfcSE, pvalue, sig))

p_volcano_B <- plot_volcano(res_B, "FMT group: Baseline → 6M post-FMT")
p_lfc_B     <- plot_lfc_bar(res_B, "Top diff. genera: FMT Baseline→6M")
p_ma_B      <- plot_ma(res_B,      "MA plot: FMT Baseline→6M")

write.csv(res_B, "outputs/tax_da/deseq2_B_FMT_baseline_vs_6M.csv", row.names = FALSE)

# don't have any bac significant in result B
# ps_6m_genus %>% 
# =============================================================================
#  SECTION C — Placebo 組內：Baseline vs 6M（自然波動對照）
# =============================================================================

cat("\n── Section C: Placebo group — Baseline vs 6M ──────────────────────────\n")

ps_plc <- subset_samples(ps_filt,
                         Treatment == "placebo" &
                           Timepoint %in% c("Baseline", "6M post-FMT"))

ps_plc_genus <- tax_glom(ps_plc, taxrank = "Genus", NArm = FALSE)
ps_plc_genus <- filter_ps_by_prev(ps_obj = ps_plc_genus)

sample_data(ps_plc_genus)$Timepoint <- relevel(
  factor(sample_data(ps_plc_genus)$Timepoint,
         levels = c("Baseline", "6M post-FMT")),
  ref = "Baseline"
)

dds_plc <- make_dds(ps_plc_genus, ~ PatientID + Timepoint)
dds_plc <- DESeq(dds_plc, test = "Wald", fitType = "parametric", quiet = TRUE) # Wald 檢定（逐一兩兩比較）

res_C <- tidy_deseq2(dds_plc, ps_plc_genus,
                     contrast = c("Timepoint", "6M post-FMT", "Baseline"))

cat(sprintf("Significant (p<0.05 & |LFC|>0): %d\n",
            sum(res_C$direction != "NS")))
print(filter(res_C, pvalue < 0.05) %>%
        dplyr::select(Genus, log2FoldChange, lfcSE, pvalue, sig))

p_volcano_C <- plot_volcano(res_C, "Placebo group: Baseline → 6M post-FMT")
p_lfc_C     <- plot_lfc_bar(res_C, "Top diff. genera: Placebo Baseline→6M")

write.csv(res_C, "outputs/tax_da/deseq2_C_placebo_baseline_vs_6M.csv", row.names = FALSE)


# =============================================================================
#  SECTION D — 三時間點整體趨勢（LRT）Likelihood Ratio Test
# =============================================================================
# DESeq2 透過 LRT 去計算：這個菌在任何一個時間點有發生顯著波動嗎？
cat("\n── Section D: Timepoint trend — LRT ───────────────────────────────────\n")

ps_all       <- subset_samples(ps_filt, PatientID != "GM03")
ps_all_genus <- tax_glom(ps_all, taxrank = "Genus", NArm = FALSE)
ps_all_genus <- filter_ps_by_prev(ps_obj = ps_all_genus)


sample_data(ps_all_genus)$Timepoint <- factor(
  sample_data(ps_all_genus)$Timepoint,
  levels = c("Baseline", "2M post-FMT", "6M post-FMT")
)

run_lrt <- function(ps_obj, trt_label) {
  keep  <- sample_names(ps_obj)[sample_data(ps_obj)$Treatment == trt_label] # WHY: subset_samples()
  ps_sub <- prune_samples(keep, ps_obj)
  ps_g   <- tax_glom(ps_sub, taxrank = "Genus", NArm = FALSE)
  
  dds <- make_dds(ps_g, ~ PatientID + Timepoint) 
  dds <- DESeq(dds, test = "LRT", reduced = ~ PatientID,# Timepoint／多組別的綜合檢定（Omnibus Test）。
               fitType = "parametric", quiet = TRUE)
  
  tax_df <- as.data.frame(as(tax_table(ps_g), "matrix")) %>%
    tibble::rownames_to_column("ASV_ID")
  
  results(dds, alpha = 0.05) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ASV_ID") %>%
    left_join(tax_df, by = "ASV_ID") %>%
    filter(!is.na(pvalue)) %>%
    mutate(Treatment = trt_label,
           sig = case_when(
             pvalue < 0.001 ~ "***",
             pvalue < 0.01  ~ "**",
             pvalue < 0.05  ~ "*",
             TRUE           ~ "ns"
           )) %>%
    arrange(pvalue)
}

res_D_fmt <- run_lrt(ps_all_genus, "FMT")
res_D_plc <- run_lrt(ps_all_genus, "placebo")

cat(sprintf("LRT sig — FMT: %d  |  placebo: %d\n",
            sum(res_D_fmt$pvalue < 0.05, na.rm = TRUE),
            sum(res_D_plc$pvalue < 0.05, na.rm = TRUE)))

write.csv(res_D_fmt, "outputs/tax_da/deseq2_D_LRT_FMT_timepoint.csv",     row.names = FALSE)
write.csv(res_D_plc, "outputs/tax_da/deseq2_D_LRT_placebo_timepoint.csv", row.names = FALSE)


# =============================================================================
#  SECTION E — FMT-specific taxa（A ∩ B, not C）+ |LFC| > 1
# =============================================================================

cat("\n── Section E: FMT-specific taxa ────────────────────────────────────────\n")

fmt_specific <- res_A %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > 0) %>%
  semi_join(
    res_B %>% filter(pvalue < 0.05, abs(log2FoldChange) > 0), 
    by = "Genus" # 這些菌不僅在 6M 時讓 FMT 組跟對照組看起來不一樣，而且它們確實是在 FMT 治療後才發生改變的。
  ) %>%
  anti_join(
    res_C %>% filter(pvalue < 0.05),
    by = "Genus" # 排除掉那些「不需要灌糞便，一般人隨著時間也會自然上下波動」的背景雜訊菌。
  )

cat(sprintf("FMT-specific genera (p<0.05 & |LFC|>0): %d\n", nrow(fmt_specific)))
print(fmt_specific %>% dplyr::select(Genus, log2FoldChange, pvalue, sig))

write.csv(fmt_specific, "outputs/tax_da/deseq2_E_FMT_specific_taxa.csv", row.names = FALSE)


# =============================================================================
#  SECTION F — MaAsLin2
#  主要優勢：支援 random effect (PatientID)、continuous values、repeated measures
# =============================================================================

cat("\n── Section F: MaAsLin2 ─────────────────────────────────────────────────\n")

# genus relative abundance matrix（rows = samples for MaAsLin2）
ps_genus_rel <- tax_glom(ps_filt, taxrank = "Genus") %>%
  transform_sample_counts(function(x) x / sum(x))
ps_genus_rel <- filter_ps_by_prev(ps_obj = ps_genus_rel)

genus_mat_full <- as.data.frame(t(as.matrix(otu_table(ps_genus_rel))))
colnames(genus_mat_full) <- as.data.frame(tax_table(ps_genus_rel))$Genus

run_maaslin2 <- function(ps_rel, meta_sub, label,
                         fixed_effects, random_effects,
                         reference, output_dir) {
  common   <- intersect(rownames(meta_sub), sample_names(ps_rel))
  mat_sub  <- as.data.frame(t(as.matrix(otu_table(ps_rel))))[common, ]
  colnames(mat_sub) <- as.data.frame(tax_table(ps_rel))$Genus
  meta_sub <- meta_sub[common, , drop = FALSE]
  
  fit <- Maaslin2(
    input_data      = mat_sub,
    input_metadata  = meta_sub,
    output          = output_dir,
    fixed_effects   = fixed_effects,
    random_effects  = random_effects,
    reference       = reference,
    normalization   = "NONE",   # already relative abundance
    transform       = "LOG",
    analysis_method = "LM",
    min_prevalence  = 0.1,
    min_abundance   = 0.0,
    max_significance = 0.25,
    plot_heatmap    = FALSE,
    plot_scatter    = FALSE
  )
  
  fit$results %>%
    as_tibble() %>%
    rename(Genus = feature) %>%
    mutate(method = "MaAsLin2", label = label) %>%
    arrange(qval)
}

# ── F-A. FMT vs Placebo @ 6M ─────────────────────────────────────────────────

meta_6m <- remove_rownames(meta_raw) %>%
  filter(Timepoint == "6M post-FMT") %>%
  mutate(Treatment = factor(Treatment, levels = c("placebo", "FMT"))) %>%
  column_to_rownames("SampleID")

ps_6m_rel <- subset_samples(ps_genus_rel, Timepoint == "6M post-FMT")

maas_A <- run_maaslin2(
  ps_6m_rel, meta_6m, "FMT_vs_placebo_6M",
  fixed_effects  = c("Treatment", "Age"),
  random_effects = character(0),          # 6M only，無 repeated measures
  reference      = c("Treatment,placebo"),
  output_dir     = file.path(getwd(), "outputs/tax_da/", "maaslin2_A")
)

cat(sprintf("MaAsLin2 [A] q<0.25: %d  q<0.05: %d\n",
            sum(maas_A$qval < 0.25, na.rm = TRUE),
            sum(maas_A$qval < 0.05, na.rm = TRUE)))

# ── F-B. FMT 組內：Baseline vs 6M（paired）────────────────────────────────────

meta_fmt_b <- remove_rownames(meta_raw) %>%
  filter(Treatment == "FMT",
         Timepoint %in% c("Baseline", "6M post-FMT"),
         PatientID != "GM03") %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Baseline", "6M post-FMT"))) %>%
  column_to_rownames("SampleID")

ps_fmt_rel <- subset_samples(ps_genus_rel,
                             Treatment == "FMT" &
                               Timepoint %in% c("Baseline", "6M post-FMT") &
                               PatientID != "GM03")

maas_B <- run_maaslin2(
  ps_fmt_rel, meta_fmt_b, "FMT_BL_vs_6M",
  fixed_effects  = c("Timepoint", "Age"),
  random_effects = c("PatientID"),
  reference      = c("Timepoint,Baseline"),
  output_dir     = file.path(getwd(), "outputs/tax_da/", "maaslin2_B")
)

cat(sprintf("MaAsLin2 [B] q<0.25: %d  q<0.05: %d\n",
            sum(maas_B$qval < 0.25, na.rm = TRUE),
            sum(maas_B$qval < 0.05, na.rm = TRUE)))

# ── F-C. Placebo 組內（對照）──────────────────────────────────────────────────

meta_plc_c <- remove_rownames(meta_raw) %>%
  filter(Treatment == "placebo",
         Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("Baseline", "6M post-FMT"))) %>%
  column_to_rownames("SampleID")

ps_plc_rel <- subset_samples(ps_genus_rel,
                             Treatment == "placebo" &
                               Timepoint %in% c("Baseline", "6M post-FMT"))

maas_C <- run_maaslin2(
  ps_plc_rel, meta_plc_c, "placebo_BL_vs_6M",
  fixed_effects  = c("Timepoint", "Age"),
  random_effects = c("PatientID"),
  reference      = c("Timepoint,Baseline"),
  output_dir     = file.path(getwd(), "outputs/tax_da/", "maaslin2_C")
)

cat(sprintf("MaAsLin2 [C] q<0.25: %d  q<0.05: %d\n",
            sum(maas_C$qval < 0.25, na.rm = TRUE),
            sum(maas_C$qval < 0.05, na.rm = TRUE)))

# ── F-D. MaAsLin2 FMT-specific（A ∩ B, not C）────────────────────────────────
maas_A_sig <- maas_A %>% filter(metadata == "Treatment", coef > 0, qval < 0.25)
maas_B_sig <- maas_B %>% filter(metadata == "Timepoint", coef > 0, qval < 0.25)
maas_C_sig <- maas_C %>% filter(metadata == "Timepoint", qval < 0.25)

fmt_specific_maas <- intersect(maas_A_sig$Genus, maas_B_sig$Genus) %>%
  setdiff(maas_C_sig$Genus)

cat(sprintf("MaAsLin2 FMT-specific genera: %d\n", length(fmt_specific_maas)))
print(fmt_specific_maas)

write.csv(maas_A, "outputs/tax_da/maaslin2_A_FMT_vs_placebo_6M.csv", row.names = FALSE)
write.csv(maas_B, "outputs/tax_da/maaslin2_B_FMT_BL_vs_6M.csv",      row.names = FALSE)
write.csv(maas_C, "outputs/tax_da/maaslin2_C_placebo_BL_vs_6M.csv",  row.names = FALSE)

# Safe defaults for objects defined in downstream scripts (must precede all viz code)
if (!exists("lefse_tab")     || !is.data.frame(lefse_tab))    lefse_tab    <- data.frame()
if (!exists("triple_sig")    || !is.data.frame(triple_sig))   triple_sig   <- data.frame()
if (!exists("p_lefse"))      p_lefse      <- ggplot() + labs(title = "LEfSe (not run)")
if (!exists("p_concordance")) p_concordance <- ggplot() + labs(title = "Concordance (not run)")

# MaAsLin2 volcano (defined here; same function in 4-pathway_analysis.R)
plot_maaslin_volcano <- function(res_df, coef_var = "Treatment", title = "") {
  df <- res_df %>%
    filter(metadata == coef_var) %>%
    mutate(
      direction = case_when(
        qval < 0.25 & coef > 0 ~ "Up",
        qval < 0.25 & coef < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label_pw = if_else(qval < 0.1, Pathway, NA_character_)
    )
  ggplot(df, aes(x = coef, y = -log10(qval), colour = direction)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_text_repel(aes(label = label_pw), size = 2.5,
                    max.overlaps = 15, segment.colour = "grey60",
                    fontface = "bold", na.rm = TRUE) +
    geom_hline(yintercept = -log10(0.25), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "#D85A30", linewidth = 0.4) +
    scale_colour_manual(values = c(Up = "#D6604D", Down = "#2166AC", NS = "grey70")) +
    labs(title = title,
         subtitle = sprintf("q<0.25: %d up, %d down",
                            sum(df$direction == "Up"), sum(df$direction == "Down")),
         x = "Coefficient", y = "-log10(q-value)", colour = NULL) +
    theme_fmt
}

p_vol_maas_A <- plot_maaslin_volcano(maas_A %>% rename(Pathway="Genus"),
                                      title = "TAXONOMY: FMT v.s. Placebo @ 6M Post")

p_vol_maas_B <- plot_maaslin_volcano(maas_B %>% rename(Pathway="Genus"),
                                     coef_var = "Timepoint",
                                     title = "TAXONOMY: FMT (Baseline -> 6M Post)")

p_vol_maas_C <- plot_maaslin_volcano(maas_C %>% rename(Pathway="Genus"),
                                     coef_var = "Timepoint",
                                     title = "TAXONOMY: Placebo (Baseline -> 6M Post)")

# =============================================================================
#  SECTION G — 視覺化
# =============================================================================

# ── G1. Volcano panel ─────────────────────────────────────────────────────────

p_volcano_panel <- (p_volcano_A | p_volcano_B) /
  (p_volcano_C | plot_spacer()) +
  plot_annotation(
    title    = "DESeq2 differential abundance — all comparisons",
    subtitle = "Dashed lines: |LFC|=1, p=0.05",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# ── G2. LFC bar panel ────────────────────────────────────────────────────────

if (!is.null(p_lfc_A) && !is.null(p_lfc_B)) {
  p_lfc_panel <- (p_lfc_A | p_lfc_B) +
    plot_annotation(
      title = "Top differential genera: FMT vs placebo @ 6M vs FMT temporal",
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
}

# ── G3. FMT-specific taxa trajectory ─────────────────────────────────────────

p_fmt_specific <- NULL
#res_A_sig <- filter(res_A, sig != 'ns')
if (nrow(fmt_specific) > 0) {
  ps_genus_rel <- tax_glom(ps_filt, taxrank = "Genus") %>%
    transform_sample_counts(function(x) x / sum(x) * 100)
  
  abund_df <- psmelt(ps_genus_rel) %>%
    filter(Genus %in% fmt_specific$Genus) %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_labels)) %>%
    group_by(Genus, Treatment, Timepoint) %>%
    summarise(mean_abund = mean(Abundance),
              se_abund   = sd(Abundance) / sqrt(n()),
              .groups = "drop")
  
  p_fmt_specific <- ggplot(abund_df,
                           aes(x = Timepoint, y = mean_abund,
                               colour = Treatment, group = Treatment)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_abund - se_abund,
                      ymax = mean_abund + se_abund),
                  width = 0.2) +
   # facet_wrap(~ Genus, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = col_treatment) +
    labs(title    = "FMT-specific taxa: relative abundance trajectory",
         subtitle = "Mean ± SE  |  DESeq2 A∩B not C  |  p<0.05 & |LFC|>0",
         x = NULL, y = "Relative abundance (%)") +
    theme_fmt +
    theme(axis.text.x    = element_text(angle = 25, hjust = 1),
          strip.text      = element_text(face = "italic", size = 9),
          legend.position = "bottom")
}

# ── G4. Heatmap @ 6M ─────────────────────────────────────────────────────────

p_heatmap_A <- NULL
sig_genera_A <- res_A %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > 1) %>%
  pull(Genus) %>% na.omit() %>% unique()

if (length(sig_genera_A) >= 2) {
  ps_6m_rel <- transform_sample_counts(ps_6m_genus, function(x) x / sum(x) * 100)
  
  hm_df <- psmelt(ps_6m_rel) %>%
    filter(Genus %in% sig_genera_A) %>%
    left_join(res_A %>% dplyr::select(Genus, log2FoldChange, pvalue, sig),
              by = "Genus") %>%
    mutate(log_abund = log10(Abundance + 0.001),
           Genus     = reorder(Genus, log2FoldChange))
  
  p_heatmap_A <- ggplot(hm_df, aes(x = SampleID, y = Genus, fill = log_abund)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_grid(. ~ Treatment, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = median(hm_df$log_abund),
      name = "log₁₀(rel. abund + 0.001)"
    ) +
    labs(title    = "Significant genera @ 6M (p<0.05 & |LFC|>0)",
         subtitle = "Sorted by log₂FC (FMT/placebo)",
         x = NULL, y = NULL) +
    theme_fmt +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y     = element_text(face = "italic", size = 9),
          legend.position = "right")
}

# ── Safe defaults for LEfSe objects (defined in deseq2_lefse_addon.R if run) ─
if (!exists("lefse_tab"))    lefse_tab    <- data.frame()
if (!exists("triple_sig"))   triple_sig   <- data.frame()
if (!exists("p_lefse"))      p_lefse      <- ggplot() + labs(title = "LEfSe (not run)")
if (!exists("p_concordance")) p_concordance <- ggplot() + labs(title = "Concordance (not run)")

# ── G5. Triple concordant summary table ──────────────────────────────────────

if (nrow(triple_sig) > 0) {
  p_triple <- ggplot(triple_sig,
                     aes(x = reorder(Genus, ef_lda), y = ef_lda,
                         fill = enrich_group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("LFC=%.1f", log2FoldChange)),
              hjust = if_else(triple_sig$ef_lda > 0, -0.1, 1.1),
              size  = 3) +
    coord_flip() +
    scale_fill_manual(values = col_treatment) +
    labs(title    = "Triple-concordant taxa (LEfSe + DESeq2_A + DESeq2_B)",
         subtitle = "Most reliable FMT-driven genera",
         x = NULL, y = "LEfSe LDA score", fill = "Enriched in") +
    theme_fmt +
    theme(axis.text.y = element_text(face = "italic"))
}


# =============================================================================
#  Save
# =============================================================================
dir.create(path = "outputs/tax_da/")
save_fig <- function(plot, filename, w = 10, h = 6) {
  if (is.null(plot)) { message("Skipping NULL: ", filename); return(invisible()) }
  ggsave(filename, plot = plot, width = w, height = h, dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", filename))
}

save_fig(p_volcano_A,     "outputs/tax_da/deseq2_fig1a_volcano_FMT_vs_placebo_6M.pdf")
save_fig(p_volcano_B,     "outputs/tax_da/deseq2_fig1b_volcano_FMT_baseline_6M.pdf")
save_fig(p_volcano_C,     "deseq2_fig1c_volcano_placebo_baseline_6M.pdf")
save_fig(p_volcano_panel, "outputs/tax_da/deseq2_fig1_volcano_panel.pdf",           w = 14, h = 10)
save_fig(p_ma_A,          "outputs/tax_da/deseq2_fig2a_MA_FMT_vs_placebo.pdf")
save_fig(p_ma_B,          "outputs/tax_da/deseq2_fig2b_MA_FMT_temporal.pdf")
save_fig(p_lfc_A,         "outputs/tax_da/deseq2_fig3a_lfc_FMT_vs_placebo.pdf",    w = 10, h = 8)
save_fig(p_lfc_B,         "outputs/tax_da/deseq2_fig3b_lfc_FMT_temporal.pdf",      w = 10, h = 8)
if (exists("p_lfc_panel"))
  save_fig(p_lfc_panel,   "outputs/tax_da/deseq2_fig3_lfc_panel.pdf",              w = 16, h = 8)
save_fig(p_fmt_specific,  "outputs/tax_da/deseq2_fig4_FMT_specific_taxa.pdf",      w = 14, h = 8)
save_fig(p_heatmap_A,     "outputs/tax_da/deseq2_fig5_heatmap_6M.pdf",             w = 12, h = 7)
save_fig(p_lefse,         "outputs/tax_da/deseq2_fig6_lefse_delta.pdf",            w = 10, h = 8)
save_fig(p_concordance,   "outputs/tax_da/deseq2_fig7_deseq2_lefse_concordance.pdf", w = 9, h = 7)
if (exists("p_triple"))
  save_fig(p_triple,      "outputs/tax_da/deseq2_fig8_triple_concordant.pdf",      w = 9,  h = 6)

cat("\n✓ DESeq2 + LEfSe analysis complete.\n")
cat("  Key outputs:\n")
cat(sprintf("    fmt_specific : DESeq2 FMT-specific taxa        n=%d\n", nrow(fmt_specific)))
cat(sprintf("    lefse_tab    : LEfSe significant genera         n=%d\n", nrow(lefse_tab)))
cat(sprintf("    triple_sig   : Triple-concordant (most reliable) n=%d\n", nrow(triple_sig)))
cat("  CSV: deseq2_A/B/C/D/E_*.csv | lefse_*.csv\n")