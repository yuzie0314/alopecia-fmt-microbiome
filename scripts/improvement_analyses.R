# =============================================================================
#  improvement_analyses.R
#  下一步改進分析 — Alopecia Oral FMT Study
#
#  PART A — 關鍵發表圖（publication-ready figures）
#  PART B — 統計強化（ALDEx2 mc.samples, LEfSe on Δabundance）
#  PART C — SALT-pathway 深度分析
#
#  前置條件：run_pipeline.R 已成功執行 Step 0 + Step 2
#           （ps_filt, meta_raw, metacyc_raw, res_ait 已載入）
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

library(tidyverse)
library(ggrepel)
library(patchwork)
library(vegan)
library(phyloseq)

select <- dplyr::select
OUT_DIR <- "outputs/improvements"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
#  PART A — 關鍵發表圖
# =============================================================================

cat("\n── PART A: 發表圖製作 ─────────────────────────────────────────────────\n")

# ── A1. ΔPC1 (Aitchison) vs ΔSALT — 核心臨床圖 ─────────────────────────────
# res_ait$p_salt 已由 2-beta.R 計算完成；重存為 publication 版本

tryCatch({
  if (!exists("res_ait")) stop("res_ait 不存在，請先執行 2-beta.R")

  p_A1 <- res_ait$p_salt +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40"),
      legend.position = "top"
    ) +
    labs(
      title    = "Community shift correlates with hair regrowth",
      subtitle = "Aitchison ΔPC1 ~ ΔSALT  |  Spearman ρ = 0.63, p = 0.067",
      caption  = "Positive ΔSALT = worsening; negative ΔSALT = improvement. GM03 excluded."
    )

  ggsave(file.path(OUT_DIR, "figA1_dPC1_vs_dSALT.pdf"),
         p_A1, width = 6, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "figA1_dPC1_vs_dSALT.png"),
         p_A1, width = 6, height = 5.5, dpi = 300, bg = "white")
  cat("  Saved: figA1_dPC1_vs_dSALT\n")
}, error = function(e) message("  [SKIP figA1] ", e$message))

# ── A2. Klebsiella 逐病人追蹤圖 ────────────────────────────────────────────

tryCatch({
  if (!exists("ps_filt")) stop("ps_filt 不存在，請先執行 0_1-loadData.R")

  ps_genus <- tax_glom(ps_filt, "Genus", NArm = FALSE)
  ps_rel   <- transform_sample_counts(ps_genus, function(x) x / sum(x) * 100)

  kleb_asv <- tax_table(ps_rel) %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    filter(Genus == "Klebsiella") %>%
    pull(ASV)

  if (length(kleb_asv) == 0) stop("Klebsiella ASV not found in ps_filt")

  tp_levels <- c("Baseline", "2M post-FMT", "6M post-FMT")
  kleb_df <- psmelt(prune_taxa(kleb_asv, ps_rel)) %>%
    group_by(PatientID, Treatment, Timepoint) %>%
    summarise(rel_abund = sum(Abundance), .groups = "drop") %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_levels))

  p_A2 <- ggplot(kleb_df, aes(x = Timepoint, y = rel_abund,
                                group = PatientID, colour = Treatment)) +
    geom_line(aes(linetype = Treatment), linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 3) +
    geom_text_repel(
      data = filter(kleb_df, Timepoint == "6M post-FMT"),
      aes(label = PatientID), size = 2.8, nudge_x = 0.15,
      segment.colour = "grey70", show.legend = FALSE
    ) +
    scale_colour_manual(values = c(FMT = "#D6604D", placebo = "#4393C3")) +
    scale_y_log10(
      name   = "Relative abundance (%, log scale)",
      labels = function(x) sprintf("%.2f%%", x)
    ) +
    scale_linetype_manual(values = c(FMT = "solid", placebo = "dashed")) +
    labs(
      title    = "Klebsiella relative abundance — per patient trajectory",
      subtitle = "DESeq2 Comparison A: LFC = −2.21, p = 0.00075 ***",
      x = NULL, colour = NULL, linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "top",
      axis.text.x     = element_text(angle = 20, hjust = 1)
    )

  ggsave(file.path(OUT_DIR, "figA2_Klebsiella_trajectory.pdf"),
         p_A2, width = 6, height = 4.5, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "figA2_Klebsiella_trajectory.png"),
         p_A2, width = 6, height = 4.5, dpi = 300, bg = "white")
  cat("  Saved: figA2_Klebsiella_trajectory\n")
}, error = function(e) message("  [SKIP figA2] ", e$message))

# ── A3. SALT-correlated pathways — volcano plot ──────────────────────────────
# pathway_salt_correlation_metacyc.csv 用 MetaCyc ID（PWY-xxx）格式
# 直接用 rho vs -log10(padj) 呈現，不需要 cross-join MaAsLin2 資料

tryCatch({
  salt_metacyc_path <- "outputs/pathway/pathway_salt_correlation_metacyc.csv"
  if (!file.exists(salt_metacyc_path)) stop("pathway_salt_correlation_metacyc.csv not found")

  salt_res <- read_csv(salt_metacyc_path, show_col_types = FALSE) %>%
    mutate(
      direction = case_when(
        rho < -0.7 & padj < 0.05 ~ "Beneficial (↑ pathway → ↓ SALT)",
        rho >  0.7 & padj < 0.05 ~ "Adverse (↑ pathway → ↑ SALT)",
        padj < 0.05               ~ "Significant",
        TRUE                      ~ "ns"
      ),
      sig_label = ifelse(padj < 0.01 & abs(rho) > 0.9,
                         str_trunc(Pathway, 35), NA_character_)
    )

  n_sig <- sum(salt_res$padj < 0.05, na.rm = TRUE)
  cat(sprintf("  figA3: %d pathways total, %d significant (padj<0.05)\n",
              nrow(salt_res), n_sig))

  p_A3 <- ggplot(salt_res, aes(x = rho, y = -log10(padj), colour = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.7, 0.7), linetype = "dashed", colour = "grey80", linewidth = 0.3) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = sig_label), size = 2.2, max.overlaps = 20,
                    na.rm = TRUE, segment.colour = "grey60", segment.size = 0.3) +
    scale_colour_manual(
      values = c(
        "Beneficial (↑ pathway → ↓ SALT)" = "#D6604D",
        "Adverse (↑ pathway → ↑ SALT)"    = "#2166AC",
        "Significant"                       = "#F4A582",
        "ns"                                = "grey75"
      )
    ) +
    annotate("text", x = -0.9, y = 0.2, label = sprintf("n = %d sig", n_sig),
             size = 3, colour = "grey40") +
    labs(
      title    = "MetaCyc pathways: SALT score correlation",
      subtitle = sprintf("Spearman ρ (pathway Δabundance ~ ΔSALT)  |  %d pathways significant (BH padj<0.05)", n_sig),
      x        = "Spearman ρ with ΔSALT (negative = pathway ↑ → improvement)",
      y        = "-log10(BH-adjusted p)",
      colour   = NULL,
      caption  = "Dashed: padj=0.05 threshold; |ρ|=0.7 reference lines"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text     = element_text(size = 8)
    )

  ggsave(file.path(OUT_DIR, "figA3_pathway_SALT_volcano.pdf"),
         p_A3, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "figA3_pathway_SALT_volcano.png"),
         p_A3, width = 8, height = 6, dpi = 300, bg = "white")
  cat("  Saved: figA3_pathway_SALT_volcano\n")
}, error = function(e) message("  [SKIP figA3] ", e$message))

# ── A4. fgsea MetaCyc category bar (publication) ────────────────────────────

tryCatch({
  fgsea_path <- "outputs/functional_enrichment/fgsea_maaslin2_metacyc_tval.csv"
  if (!file.exists(fgsea_path)) stop("fgsea_maaslin2_metacyc_tval.csv not found")

  fgsea_df <- read_csv(fgsea_path, show_col_types = FALSE) %>%
    mutate(
      pathway   = str_replace_all(pathway, "-", " "),
      sig_label = ifelse(padj < 0.05, "*", ifelse(padj < 0.2, "†", ""))
    )

  p_A4 <- ggplot(fgsea_df, aes(x = NES, y = reorder(pathway, NES), fill = NES > 0)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_text(aes(label = sig_label),
              hjust = ifelse(fgsea_df$NES > 0, -0.3, 1.3),
              size = 4.5, fontface = "bold") +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "#4393C3"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.2))) +
    labs(
      title    = "MetaCyc pathway categories — fgsea (FMT BL→6M)",
      subtitle = "NES: normalised enrichment score  |  * padj<0.05  † padj<0.2",
      x        = "NES (Normalised Enrichment Score)",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_DIR, "figA4_fgsea_metacyc_bar.pdf"),
         p_A4, width = 7, height = 4, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "figA4_fgsea_metacyc_bar.png"),
         p_A4, width = 7, height = 4, dpi = 300, bg = "white")
  cat("  Saved: figA4_fgsea_metacyc_bar\n")
}, error = function(e) message("  [SKIP figA4] ", e$message))

# =============================================================================
#  PART B — 統計強化
# =============================================================================

cat("\n── PART B: 統計強化 ───────────────────────────────────────────────────\n")

# ── B1. LEfSe on Δabundance（新設計）──────────────────────────────────────
# 每個病人 Baseline→6M 的 Δlog2(rel.abund) 作為一個獨立觀察值
# 消除 repeated-measures 問題

# B1. Δabundance comparison: Wilcoxon FMT vs Placebo (per-patient paired design)
# Note: lefser 1.22 incompatible with log2-FC values (negative) → use Wilcoxon instead

tryCatch({
  if (!exists("ps_filt")) stop("ps_filt 不存在")
  if (!exists("meta_raw")) stop("meta_raw 不存在")

  cat("B1. Δabundance Wilcoxon: FMT vs Placebo (paired at data level)...\n")

  ps_genus <- tax_glom(ps_filt, "Genus", NArm = FALSE)
  eps      <- 1e-5

  bl_meta <- meta_raw %>% filter(Timepoint == "Baseline",    PatientID != "GM03") %>% arrange(PatientID)
  m6_meta <- meta_raw %>% filter(Timepoint == "6M post-FMT", PatientID != "GM03") %>% arrange(PatientID)

  avail_bl <- intersect(bl_meta$SampleID, sample_names(ps_genus))
  avail_6m <- intersect(m6_meta$SampleID, sample_names(ps_genus))

  ps_bl <- prune_samples(avail_bl, ps_genus)
  ps_6m <- prune_samples(avail_6m, ps_genus)

  rel_bl <- as(otu_table(transform_sample_counts(ps_bl, function(x) x / sum(x))), "matrix")
  rel_6m <- as(otu_table(transform_sample_counts(ps_6m, function(x) x / sum(x))), "matrix")

  sid_to_pt <- setNames(c(bl_meta$PatientID, m6_meta$PatientID),
                        c(bl_meta$SampleID,  m6_meta$SampleID))
  colnames(rel_bl) <- sid_to_pt[colnames(rel_bl)]
  colnames(rel_6m) <- sid_to_pt[colnames(rel_6m)]

  shared_pts <- intersect(colnames(rel_bl), colnames(rel_6m))
  shared_pts <- shared_pts[!is.na(shared_pts)]
  if (length(shared_pts) < 4) stop(sprintf("Only %d shared patients", length(shared_pts)))

  delta_mat <- log2(rel_6m[, shared_pts] + eps) - log2(rel_bl[, shared_pts] + eps)

  # Treatment per patient
  trt_vec <- meta_raw %>%
    filter(PatientID %in% shared_pts, Timepoint == "Baseline") %>%
    distinct(PatientID, Treatment) %>%
    arrange(match(PatientID, shared_pts)) %>%
    pull(Treatment)

  fmt_pts <- shared_pts[trt_vec == "FMT"]
  plc_pts <- shared_pts[trt_vec == "placebo"]

  # Wilcoxon per genus: FMT Δ vs placebo Δ
  wx_res <- apply(delta_mat, 1, function(dv) {
    w <- tryCatch(wilcox.test(dv[fmt_pts], dv[plc_pts], exact = FALSE),
                  error = function(e) list(statistic = NA, p.value = 1))
    c(W = unname(w$statistic),
      pval = w$p.value,
      mean_fmt = mean(dv[fmt_pts], na.rm = TRUE),
      mean_plc = mean(dv[plc_pts], na.rm = TRUE))
  }) %>% t() %>% as.data.frame() %>%
    rownames_to_column("ASV") %>%
    mutate(
      padj = p.adjust(pval, method = "BH"),
      LDA_proxy = mean_fmt - mean_plc,  # Δ(FMT) − Δ(placebo)
      Genus = tax_table(ps_genus)[ASV, "Genus"]
    ) %>%
    arrange(pval)

  n_sig <- sum(wx_res$padj < 0.25, na.rm = TRUE)
  cat(sprintf("  Wilcoxon Δabundance: %d genera with padj<0.25, %d with padj<0.05\n",
              n_sig, sum(wx_res$padj < 0.05, na.rm = TRUE)))

  write.csv(wx_res, file.path(OUT_DIR, "wilcoxon_delta_abundance_FMTvsPLC.csv"), row.names = FALSE)

  # Plot top 20 by |LDA_proxy|, regardless of significance
  top_wx <- wx_res %>%
    filter(!is.na(Genus), Genus != "NA") %>%
    slice_max(abs(LDA_proxy), n = 20)

  p_B1 <- ggplot(top_wx,
                 aes(x = reorder(Genus, LDA_proxy), y = LDA_proxy,
                     fill = LDA_proxy > 0)) +
    geom_col(width = 0.75, alpha = 0.85) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "#4393C3"),
                      labels = c("TRUE" = "FMT ↑", "FALSE" = "Placebo ↑"),
                      name = NULL) +
    labs(
      title    = "Δabundance (BL→6M): FMT vs Placebo per-patient",
      subtitle = "Each patient = 1 data point (paired design removes repeated-measures issue)\nBar = mean(ΔFMT) − mean(ΔPlacebo)",
      x = NULL, y = "Δlog2(rel.abund) FMT − Placebo"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  ggsave(file.path(OUT_DIR, "figB1_delta_abundance_FMTvsPLC.pdf"),
         p_B1, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "figB1_delta_abundance_FMTvsPLC.png"),
         p_B1, width = 7, height = 6, dpi = 300, bg = "white")
  cat("  Saved: figB1_delta_abundance_FMTvsPLC\n")
}, error = function(e) message("  [SKIP B1] ", e$message))

# ── B2. ALDEx2 mc.samples = 512 ─────────────────────────────────────────────

tryCatch({
  if (!exists("ps_filt")) stop("ps_filt 不存在")
  if (!requireNamespace("ALDEx2", quietly = TRUE)) stop("ALDEx2 未安裝")
  library(ALDEx2)

  cat("B2. ALDEx2 mc.samples=512 (FMT vs Placebo @6M)...\n")

  ps_6m_all <- subset_samples(ps_filt, Timepoint == "6M post-FMT")
  ps_6m_g   <- tax_glom(ps_6m_all, "Genus", NArm = FALSE)
  count_mat  <- as(otu_table(ps_6m_g), "matrix")
  conditions <- as.character(sample_data(ps_6m_g)$Treatment)

  clr_512 <- aldex.clr(round(count_mat), conds = conditions,
                        mc.samples = 512, denom = "all", verbose = FALSE)
  tt_512  <- aldex.ttest(clr_512, paired.test = FALSE, verbose = FALSE)
  eff_512 <- aldex.effect(clr_512, CI = TRUE, verbose = FALSE)

  aldex_res <- cbind(tt_512, eff_512) %>%
    rownames_to_column("ASV_ID") %>%
    left_join(
      as.data.frame(as(tax_table(ps_6m_g), "matrix")) %>% rownames_to_column("ASV_ID"),
      by = "ASV_ID"
    ) %>%
    arrange(wi.eBH)

  write.csv(aldex_res, file.path(OUT_DIR, "aldex2_mc512_FMT_vs_placebo_6M.csv"),
            row.names = FALSE)
  n_sig <- sum(aldex_res$wi.eBH < 0.05, na.rm = TRUE)
  cat(sprintf("  ALDEx2 mc=512: %d genera significant (wi.eBH<0.05)\n", n_sig))
  cat("  Saved: aldex2_mc512_FMT_vs_placebo_6M.csv\n")
}, error = function(e) message("  [SKIP B2] ", e$message))

# =============================================================================
#  PART C — SALT-pathway 深度分析
# =============================================================================

cat("\n── PART C: SALT-pathway 深度分析 ───────────────────────────────────────\n")

# ── C1. Top SALT-correlated pathways heatmap ────────────────────────────────
# metacyc_raw 的 rownames 為全名格式（非 MetaCyc ID）
# 直接從 metacyc_raw 計算 SALT 相關，不依賴 pathway_salt_correlation_metacyc.csv

tryCatch({
  if (!exists("metacyc_raw")) stop("metacyc_raw 不存在，請先執行 0_1-loadData.R")
  if (!exists("meta_raw"))    stop("meta_raw 不存在")

  # ── 1. 計算 ΔSALT per patient（BL → 6M）
  salt_delta <- meta_raw %>%
    filter(Timepoint %in% c("Baseline", "6M post-FMT"), PatientID != "GM03",
           !is.na(SALT_score)) %>%
    select(PatientID, Timepoint, SALT_score) %>%
    pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
    mutate(dSALT = `6M post-FMT` - Baseline) %>%
    filter(!is.na(dSALT))

  # ── 2. FMT 組 6M 樣本的 pathway Δabundance (BL → 6M)
  fmt_bl_meta <- meta_raw %>% filter(Treatment == "FMT", Timepoint == "Baseline",    PatientID != "GM03")
  fmt_6m_meta <- meta_raw %>% filter(Treatment == "FMT", Timepoint == "6M post-FMT", PatientID != "GM03")

  shared_pts <- intersect(fmt_bl_meta$PatientID, salt_delta$PatientID)
  fmt_bl_meta <- filter(fmt_bl_meta, PatientID %in% shared_pts) %>% arrange(PatientID)
  fmt_6m_meta <- filter(fmt_6m_meta, PatientID %in% shared_pts) %>% arrange(PatientID)

  bl_sids <- fmt_bl_meta$SampleID
  m6_sids <- fmt_6m_meta$SampleID

  avail_bl <- intersect(bl_sids, colnames(metacyc_raw))
  avail_6m <- intersect(m6_sids, colnames(metacyc_raw))
  if (length(avail_bl) < 3 || length(avail_6m) < 3)
    stop(sprintf("Too few samples: BL=%d, 6M=%d", length(avail_bl), length(avail_6m)))

  mat_bl <- metacyc_raw[, avail_bl]
  mat_6m <- metacyc_raw[, avail_6m]

  # 相對豐度
  mat_bl_rel <- sweep(mat_bl, 2, colSums(mat_bl), "/")
  mat_6m_rel <- sweep(mat_6m, 2, colSums(mat_6m), "/")

  # 用 PatientID 對齊
  pts_bl <- setNames(fmt_bl_meta$PatientID, fmt_bl_meta$SampleID)[avail_bl]
  pts_6m <- setNames(fmt_6m_meta$PatientID, fmt_6m_meta$SampleID)[avail_6m]
  shared_pts2 <- intersect(pts_bl, pts_6m)

  mat_bl_rel <- mat_bl_rel[, names(pts_bl)[pts_bl %in% shared_pts2], drop = FALSE]
  mat_6m_rel <- mat_6m_rel[, names(pts_6m)[pts_6m %in% shared_pts2], drop = FALSE]
  colnames(mat_bl_rel) <- pts_bl[pts_bl %in% shared_pts2]
  colnames(mat_6m_rel) <- pts_6m[pts_6m %in% shared_pts2]
  mat_bl_rel <- mat_bl_rel[, shared_pts2]
  mat_6m_rel <- mat_6m_rel[, shared_pts2]

  eps <- 1e-6
  delta_mat <- log2(mat_6m_rel + eps) - log2(mat_bl_rel + eps)

  # ── 3. Spearman correlation per pathway vs ΔSALT
  dSALT_vec <- salt_delta$dSALT[match(shared_pts2, salt_delta$PatientID)]

  cor_res <- apply(delta_mat, 1, function(pw) {
    ct <- cor.test(pw, dSALT_vec, method = "spearman", exact = FALSE)
    c(rho = unname(ct$estimate), pval = ct$p.value)
  }) %>% t() %>% as.data.frame() %>%
    rownames_to_column("Pathway") %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(padj)

  cat(sprintf("  C1: %d pathways computed, %d sig (BH padj<0.05)\n",
              nrow(cor_res), sum(cor_res$padj < 0.05, na.rm = TRUE)))

  # ── 4. Top 20 neg + top 10 pos correlation
  top_neg <- cor_res %>% filter(rho < 0) %>% slice_min(rho, n = 20)
  top_pos <- cor_res %>% filter(rho > 0) %>% slice_max(rho, n = 10)
  top_pwy <- bind_rows(top_neg, top_pos) %>% pull(Pathway)
  if (length(top_pwy) < 3) stop("Not enough pathways for heatmap")

  # ── 5. Heatmap: FMT 組 6M 樣本豐度，按 ΔSALT 排序
  mat_hm <- log10(mat_6m_rel[top_pwy, shared_pts2, drop = FALSE] + eps + 1)
  sample_order <- shared_pts2[order(dSALT_vec)]

  hm_df <- mat_hm %>%
    as.data.frame() %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "PatientID", values_to = "log10_rel") %>%
    left_join(cor_res %>% select(Pathway, rho), by = "Pathway") %>%
    left_join(salt_delta %>% select(PatientID, dSALT), by = "PatientID") %>%
    mutate(
      Pathway   = str_trunc(Pathway, 50),
      PatientID = factor(PatientID, levels = sample_order)
    )

  p_C1 <- ggplot(hm_df, aes(x = PatientID, y = reorder(Pathway, -rho), fill = log10_rel)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                         midpoint = median(hm_df$log10_rel, na.rm = TRUE),
                         name = "log10(rel.abund)") +
    labs(
      title    = "SALT-correlated MetaCyc pathways (FMT group, sorted by ΔSALT)",
      subtitle = sprintf("Top 20 neg-corr + top 10 pos-corr  |  n=%d pathways  |  cols sorted by ΔSALT",
                         length(top_pwy)),
      x = "Patient ID (sorted by ΔSALT: left=most improved)", y = NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 6),
      plot.title  = element_text(face = "bold", size = 9)
    )

  ggsave(file.path(OUT_DIR, "figC1_SALT_pathway_heatmap.pdf"),
         p_C1, width = 10, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "figC1_SALT_pathway_heatmap.png"),
         p_C1, width = 10, height = 10, dpi = 300, bg = "white")

  write.csv(cor_res, file.path(OUT_DIR, "salt_pathway_cor_fullnames.csv"), row.names = FALSE)
  cat("  Saved: figC1_SALT_pathway_heatmap + salt_pathway_cor_fullnames.csv\n")
}, error = function(e) message("  [SKIP C1] ", e$message))

# ── C2. Antibiotic sensitivity ───────────────────────────────────────────────

tryCatch({
  if (!exists("meta_raw")) stop("meta_raw 不存在")
  abx_col <- grep("antibiotic|Antibiotic|abx|ABX", colnames(meta_raw),
                  value = TRUE, ignore.case = TRUE)
  if (length(abx_col) == 0) {
    cat("  C2: No antibiotic column found in meta_raw\n")
    cat("  Available columns:", paste(colnames(meta_raw), collapse = ", "), "\n")
  } else {
    abx_pts <- meta_raw %>%
      filter(Timepoint == "6M post-FMT",
             .data[[abx_col[1]]] %in% c(TRUE, "Yes", "yes", "Y", 1)) %>%
      pull(PatientID) %>% unique()
    cat(sprintf("  C2: Patients with antibiotics: %d (%s)\n",
                length(abx_pts), paste(abx_pts, collapse = ", ")))
    cat("  Sensitivity analysis: re-run DESeq2 Comp A excluding these patients\n")
  }
}, error = function(e) message("  [SKIP C2] ", e$message))

# =============================================================================
cat("\n✓ improvement_analyses.R complete\n")
cat(sprintf("  Outputs in: %s/\n", OUT_DIR))
cat("  figA1: ΔPC1 vs ΔSALT (publication scatter)\n")
cat("  figA2: Klebsiella per-patient trajectory\n")
cat("  figA3: Pathway FMT-effect vs SALT correlation scatter\n")
cat("  figA4: fgsea MetaCyc NES bar\n")
cat("  figB1: LEfSe on Δabundance\n")
cat("  B2:    ALDEx2 mc=512 results CSV\n")
cat("  figC1: SALT-correlated pathway heatmap\n")
