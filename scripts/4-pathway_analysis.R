# =============================================================================
#  pathway_analysis.R
#  Functional pathway analysis — Alopecia Oral FMT Study
#
#  輸入：
#    metacyc.csv  — rows = pathways (MetaCyc), cols = samples
#    ec.csv       — rows = EC numbers,         cols = samples
#
#  分析流程：
#    1. 資料載入與前處理
#    2. Beta diversity（Bray-Curtis + Aitchison）
#    3. Differential pathway（MaAsLin2 + Wilcoxon）—— 非整數資料不能用 DESeq2
#    4. Taxa–Pathway 連結（與 DESeq2 結果對照）
#    5. Pathway × SALT correlation
#
#  前置條件：
#    - meta_raw, col_treatment, col_timepoint, tp_labels, theme_fmt 已定義
#    - fmt_specific（DESeq2 Section E 結果）已存在
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(scales)

# MaAsLin2（Bioconductor）
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  BiocManager::install("Maaslin2")
}
library(Maaslin2)


# =============================================================================
#  SECTION 0 — 資料載入與前處理
# =============================================================================

cat("\n── Section 0: Load pathway data ────────────────────────────────────────\n")

# ── 0a. pwy_name_map（MetaCyc ID ↔ full name）────────────────────────────────
if (!exists("pwy_name_map")) {
  pwy_name_map <- read_tsv(
    "../emulsion/full_mapping_v4_alpha/map_metacyc-pwy_name.txt.gz",
    col_names      = c("pwy_id", "pwy_name"),
    show_col_types = FALSE
  ) %>%
    mutate(pwy_name = str_trim(pwy_name),
           pwy_name = str_remove_all(pwy_name, "&"),
           pwy_name = str_remove_all(pwy_name, ";"))
  cat(sprintf("pwy_name_map loaded: %d entries\n", nrow(pwy_name_map)))
}

# ── 0b. 讀取 ──────────────────────────────────────────────────────────────────
load_pathway <- function(path, type = "MetaCyc") {
  df <- read_csv(path, show_col_types = FALSE)
  if (type == "MetaCyc"){
    # switch full names to ids and put pwy_id to the first column
    df <- left_join(df, pwy_name_map, by = c('pathway' = 'pwy_name')) %>% 
      select(-pathway) %>% 
      filter(!is.na(pwy_id))
    id_col <- "pwy_id"
  } else {
    # 第一欄是 pathway/EC ID
    id_col <- names(df)[1]
  }

  mat <- df %>%
    column_to_rownames(id_col) %>%
    as.matrix()

  cat(sprintf("[%s] %d pathways × %d samples\n", type, nrow(mat), ncol(mat)))
  return(mat)
}

metacyc_raw <- load_pathway("inputs/metacyc_pathway_table.csv", "MetaCyc") # 58 unmapped
ec_raw      <- load_pathway("inputs/ec_pathway_table.csv",      "EC")

# ── 0b. Sample name → metadata 對應 ──────────────────────────────────────────
# 欄位名稱格式：GM01_V1, GM02_V2 等
# 對應到 metadata 的 SampleID（你的 lab_id 格式）

# 確認 metacyc 欄位名稱和 metadata SampleID 是否一致
cat("MetaCyc sample names (first 5):", head(colnames(metacyc_raw), 5), "\n")
cat("Metadata SampleID (first 5):",    head(meta_raw$SampleID, 5), "\n")

# 如果不一致，在這裡做轉換
# 例如 V1→Baseline, V2→2M post-FMT, V3→6M post-FMT
# colnames(metacyc_raw) <- colnames(metacyc_raw) %>%
#   str_replace("_V1$", "_Baseline") %>% ...

# ── 0c. 過濾低豐度 pathway ───────────────────────────────────────────────────

filter_pathway_mat <- function(mat, min_prevalence = 0.2, min_mean = 1) {
  # 移除：超過 (1-min_prevalence) 比例樣本為 0 的 pathway
  prevalence <- rowMeans(mat > 0)
  keep_prev  <- prevalence >= min_prevalence

  # 移除：平均值過低的 pathway
  keep_mean  <- rowMeans(mat) >= min_mean

  mat_filt <- mat[keep_prev & keep_mean, ]
  cat(sprintf("  After filtering: %d → %d pathways\n", nrow(mat), nrow(mat_filt)))
  return(mat_filt)
}

metacyc <- filter_pathway_mat(metacyc_raw)# 606  29 # filter_pathway_mat(metacyc_raw, min_prevalence = 0.1) %>% dim(): 607/29
ec      <- filter_pathway_mat(ec_raw)# 1477   29 # filter_pathway_mat(ec_raw, min_prevalence = 0.1) %>% dim(): 1480/29



# =============================================================================
#  SECTION 1 — Beta diversity（Pathway level）
# =============================================================================

cat("\n── Section 1: Pathway beta diversity ───────────────────────────────────\n")

run_pathway_beta <- function(mat, meta, label = "MetaCyc",
                             end_tp = "6M post-FMT", seed = 42) {
  
  # 確保樣本順序一致
  common_samples <- intersect(colnames(mat), meta$SampleID)
  mat_sub  <- mat[, common_samples]
  meta_sub <- meta %>%
    filter(SampleID %in% common_samples) %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_labels)) %>%
    arrange(match(SampleID, common_samples))
  
  # ── 1a. Bray-Curtis ──────────────────────────────────────────────────────
  
  bray_dist <- vegdist(t(mat_sub), method = "bray")
  
  pcoa_bc  <- cmdscale(bray_dist, k = 3, eig = TRUE)
  var_exp  <- pcoa_bc$eig / sum(pcoa_bc$eig[pcoa_bc$eig > 0]) * 100
  
  scores <- data.frame(
    PC1 = pcoa_bc$points[, 1],
    PC2 = pcoa_bc$points[, 2],
    PC3 = pcoa_bc$points[, 3]
  ) %>%
    tibble::rownames_to_column("SampleID") %>%
    left_join(meta_sub %>% dplyr::select(SampleID, PatientID, Timepoint, Treatment),
              by = "SampleID") %>%
    arrange(PatientID, Timepoint)
  
  cat(sprintf("[%s | Bray-Curtis] PC1: %.1f%%  PC2: %.1f%%\n",
              label, var_exp[1], var_exp[2]))
  
  perm_ctrl <- how(blocks = meta_sub$PatientID, nperm = 999)
  set.seed(seed)
  perm_res <- adonis2(
    bray_dist ~ Treatment * Timepoint + Age,
    data         = meta_sub,
    permutations = perm_ctrl,
    by           = "terms"
  )
  cat(sprintf("[%s | Bray-Curtis] PERMANOVA:\n", label))
  print(perm_res)
  
  # ── 1b. CLR + Aitchison ──────────────────────────────────────────────────
  
  mat_clr  <- apply(mat_sub + 0.001, 2, function(x) log(x) - mean(log(x)))
  ait_dist <- dist(t(mat_clr), method = "euclidean")
  
  pcoa_ait  <- cmdscale(ait_dist, k = 3, eig = TRUE)
  var_exp_a <- pcoa_ait$eig / sum(pcoa_ait$eig[pcoa_ait$eig > 0]) * 100
  
  scores_ait <- data.frame(
    PC1 = pcoa_ait$points[, 1],
    PC2 = pcoa_ait$points[, 2],
    PC3 = pcoa_ait$points[, 3]
  ) %>%
    tibble::rownames_to_column("SampleID") %>%
    left_join(meta_sub %>% dplyr::select(SampleID, PatientID, Timepoint, Treatment),
              by = "SampleID") %>%
    arrange(PatientID, Timepoint)
  
  cat(sprintf("[%s | Aitchison] PC1: %.1f%%  PC2: %.1f%%\n",
              label, var_exp_a[1], var_exp_a[2]))
  
  perm_ctrl_ait <- how(blocks = meta_sub$PatientID, nperm = 999)
  set.seed(seed)
  perm_res_ait <- adonis2(
    ait_dist ~ Treatment * Timepoint + Age,
    data         = meta_sub,
    permutations = perm_ctrl_ait,
    by           = "terms"
  )
  cat(sprintf("[%s | Aitchison] PERMANOVA:\n", label))
  print(perm_res_ait)
  
  # ── 1c. PCoA plots ───────────────────────────────────────────────────────
  
  make_pcoa_plot <- function(sc, vexp, dist_name) {
    ggplot(sc, aes(x = PC1, y = PC2, colour = Treatment)) +
      geom_path(aes(group = PatientID),
                colour = "gray40", linewidth = 0.4, linetype = "dashed",
                arrow  = arrow(length = unit(0.2, "cm"), type = "closed"),
                alpha  = 0.6) +
      geom_point(aes(shape = Timepoint), size = 3.5, alpha = 0.85) +
      stat_ellipse(aes(group = Treatment), level = 0.75, linewidth = 0.8) +
      geom_text(data = filter(sc, Timepoint == "Baseline"),
                aes(label = PatientID),
                nudge_x = 0.01, nudge_y = -0.01, size = 2.8) +
      scale_colour_manual(values = col_treatment) +
      scale_shape_manual(values  = c(16, 17, 15)) +
      labs(title  = sprintf("[%s | %s] PCoA", label, dist_name),
           x      = sprintf("PC1 (%.1f%%)", vexp[1]),
           y      = sprintf("PC2 (%.1f%%)", vexp[2]),
           colour = "Treatment", shape = "Time point") +
      theme_fmt
  }
  
  p_bray <- make_pcoa_plot(scores,     var_exp,   "Bray-Curtis")
  p_ait  <- make_pcoa_plot(scores_ait, var_exp_a, "Aitchison")
  
  # ── 1d + 1e. Loop over both distances：shift + envfit + arrow ─────────────
  
  dist_configs <- list(
    list(
      dist_name  = "Bray-Curtis",
      scores_df  = scores,
      var_e      = var_exp,
      env_mat    = t(mat_sub)          # Bray-Curtis envfit: relative abundance
    ),
    list(
      dist_name  = "Aitchison",
      scores_df  = scores_ait,
      var_e      = var_exp_a,
      env_mat    = t(mat_clr)          # Aitchison envfit: CLR matrix
    )
  )
  
  shift_results  <- list()
  envfit_results <- list()
  arrow_plots    <- list()
  
  for (cfg in dist_configs) {
    
    dist_name <- cfg$dist_name
    sc        <- cfg$scores_df
    vexp      <- cfg$var_e
    env_mat   <- cfg$env_mat
    
    cat(sprintf("\n[%s | %s] Directional shift (Baseline → %s)\n",
                label, dist_name, end_tp))
    
    # ── Shift df ────────────────────────────────────────────────────────────
    
    shift_df <- sc %>%
      filter(Timepoint %in% c("Baseline", end_tp)) %>%
      dplyr::select(PatientID, Treatment, Timepoint, PC1, PC2) %>%
      pivot_wider(names_from  = Timepoint,
                  values_from = c(PC1, PC2)) %>%
      rename(
        PC1_bl = `PC1_Baseline`,
        PC1_et = !!sym(paste0("PC1_", end_tp)),
        PC2_bl = `PC2_Baseline`,
        PC2_et = !!sym(paste0("PC2_", end_tp))
      ) %>%
      mutate(dPC1 = PC1_et - PC1_bl,
             dPC2 = PC2_et - PC2_bl) %>%
      filter(!is.na(dPC1))
    
    shift_df %>%
      group_by(Treatment) %>%
      summarise(
        mean_dPC1  = round(mean(dPC1), 4),
        n_pos      = sum(dPC1 > 0),
        n_neg      = sum(dPC1 < 0),
        consistent = all(dPC1 > 0) | all(dPC1 < 0),
        .groups    = "drop"
      ) %>%
      print()
    
    wx <- wilcox.test(dPC1 ~ Treatment, data = shift_df, exact = FALSE)
    cat(sprintf("[%s | %s] Wilcoxon dPC1: W=%.1f, p=%.4f\n",
                label, dist_name, wx$statistic, wx$p.value))
    
    shift_results[[dist_name]] <- shift_df
    
    # ── envfit ───────────────────────────────────────────────────────────────
    
    samples_sub <- filter(sc, Timepoint %in% c("Baseline", end_tp)) %>%
      pull(SampleID)
    
    pcoa_mat <- sc %>%
      filter(SampleID %in% samples_sub) %>%
      column_to_rownames("SampleID") %>%
      dplyr::select(PC1, PC2) %>%
      as.matrix()
    
    # env_mat rows = samples，需要對齊 pcoa_mat 的 rownames
    env_sub <- env_mat[rownames(pcoa_mat), , drop = FALSE]
    
    set.seed(seed)
    fit <- envfit(pcoa_mat, env_sub, permutations = 999, na.rm = TRUE)
    
    fit_scores <- as.data.frame(scores(fit, display = "vectors")) %>%
      tibble::rownames_to_column("Pathway") %>%
      mutate(
        r2   = fit$vectors$r,
        pval = fit$vectors$pvals
      ) %>%
      filter(pval < 0.05) %>%
      arrange(desc(r2))
    
    cat(sprintf("[%s | %s] envfit significant pathways (p<0.05): %d\n",
                label, dist_name, nrow(fit_scores)))
    
    envfit_results[[dist_name]] <- fit_scores
    
    top_pw <- slice_max(fit_scores, r2, n = 10)
    
    # auto scale_factor
    
    # 1. 找出 PCoA 圖中點的實際最大邊界半徑
    pcoa_max_radius <- max(sqrt(sc$PC1^2 + sc$PC2^2), na.rm = TRUE)
    
    # 2. 如果有顯著pathway，先將 envfit 的向量「長度歸一化（變成長度1）」，再乘以想呈現的長度
    if (nrow(top_pw) > 0) {
      top_pw <- top_pw %>%
        mutate(
          # 計算原始向量長度
          orig_len = sqrt(PC1^2 + PC2^2),
          # 歸一化後，乘以 PCoA 最大半徑的 35% (可以微調 0.35 控制箭頭長短)
          PC1_scaled = (PC1 / orig_len) * (pcoa_max_radius * 0.35),
          PC2_scaled = (PC2 / orig_len) * (pcoa_max_radius * 0.35)
        )
    }
    
    # ── Arrow plot ────────────────────────────────────────────────────────────
    
    p_arrow <- ggplot() +
      geom_point(
        data = filter(sc, Timepoint == "Baseline"),
        aes(x = PC1, y = PC2, colour = Treatment),
        size = 3.5, shape = 1, stroke = 1.2
      ) +
      geom_segment(
        data = shift_df,
        aes(x = PC1_bl, y = PC2_bl,
            xend = PC1_et, yend = PC2_et,
            colour = Treatment),
        arrow     = arrow(length = unit(0.015, "npc"), type = "closed"),
        linewidth = 0.9, alpha = 0.85
      ) +
      geom_text_repel(
        data = shift_df,
        aes(x = PC1_et, y = PC2_et, label = PatientID, colour = Treatment),
        size = 3, segment.colour = NA
      ) +
      scale_colour_manual(values = col_treatment) +
      labs(
        title    = sprintf("[%s | %s] Displacement + envfit", label, dist_name),
        subtitle = sprintf("Circle = Baseline  |  Arrowhead = %s", end_tp),
        x        = sprintf("PC1 (%.1f%%)", vexp[1]),
        y        = sprintf("PC2 (%.1f%%)", vexp[2])
      ) +
      theme_fmt
    
    # ── 修正後的疊加 envfit 向量（改用已縮放的欄位） ──────────────────────
    if (nrow(top_pw) > 0) {
      p_arrow <- p_arrow +
        geom_segment(
          data = top_pw,
          aes(x = 0, y = 0,
              xend = PC1_scaled,   # 修正：改用明確縮放後的欄位
              yend = PC2_scaled),  # 修正：改用明確縮放後的欄位
          arrow       = arrow(length = unit(0.012, "npc"), type = "closed"),
          colour      = "grey20",
          linewidth   = 0.65,
          inherit.aes = FALSE
        ) +
        geom_text_repel(
          data = top_pw,
          aes(x    = PC1_scaled,   # 修正：改用明確縮放後的欄位
              y    = PC2_scaled,   # 修正：改用明確縮放後的欄位
              label = Pathway),
          colour         = "grey15",
          size           = 2.5,
          segment.colour = NA,
          inherit.aes    = FALSE
        ) +
        labs(subtitle = sprintf(
          "Circle = Baseline  |  Arrowhead = %s  |  Grey arrows = envfit (p<0.05, top 10 r²)",
          end_tp
        ))
    }
    
    arrow_plots[[dist_name]] <- p_arrow
    
  }   # end loop over dist_configs
  
  # ── 1f. 並排比較兩種距離的 arrow + envfit ────────────────────────────────
  
  p_arrow_compare <- (arrow_plots[["Bray-Curtis"]] | arrow_plots[["Aitchison"]]) +
    plot_annotation(
      title    = sprintf("[%s] Displacement + envfit: Bray-Curtis vs Aitchison", label),
      subtitle = sprintf("Baseline → %s  |  Grey arrows = driving pathways", end_tp),
      theme    = theme(plot.title = element_text(size = 12, face = "bold"))
    )
  
  # ── Return ────────────────────────────────────────────────────────────────
  
  list(
    # distances
    bray_dist      = bray_dist,
    ait_dist       = ait_dist,
    # scores
    scores         = scores,
    scores_ait     = scores_ait,
    var_exp        = var_exp,
    var_exp_a      = var_exp_a,
    # PERMANOVA
    perm           = perm_res,
    perm_ait       = perm_res_ait,
    # shift
    shift_bray     = shift_results[["Bray-Curtis"]],
    shift_ait      = shift_results[["Aitchison"]],
    # envfit
    envfit_bray    = envfit_results[["Bray-Curtis"]],
    envfit_ait     = envfit_results[["Aitchison"]],
    # plots
    p_bray         = p_bray,
    p_ait          = p_ait,
    p_arrow_bray   = arrow_plots[["Bray-Curtis"]],
    p_arrow_ait    = arrow_plots[["Aitchison"]],
    p_arrow_compare = p_arrow_compare
  )
}
beta_metacyc <- run_pathway_beta(metacyc, meta_raw, label = "MetaCyc")
beta_ec      <- run_pathway_beta(ec,      meta_raw, label = "EC")

# beta_ec$p_arrow_compare
# =============================================================================
#  SECTION 2 — Differential pathway（MaAsLin2）
#  原因：pathway 值是 continuous，不是 integer count
#  MaAsLin2 支援 continuous 值、random effect、multiple comparisons
# =============================================================================

cat("\n── Section 2: Differential pathway (MaAsLin2) ──────────────────────────\n")

run_maaslin2 <- function(mat, meta, label,
                          fixed_effects  = c("Treatment", "Timepoint", "Age"),
                          random_effects = c("PatientID"),
                          reference      = c("Treatment,placebo",
                                             "Timepoint,Baseline"),
                          output_dir     = NULL,
                          subset_expr    = NULL) {

  # 樣本篩選
  common <- intersect(colnames(mat), meta$SampleID)
  mat_sub  <- t(mat[, common])   # MaAsLin2 要求 rows = samples
  meta_sub <- meta %>%
    filter(SampleID %in% common) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = tp_labels),
      Treatment = factor(Treatment, levels = c("placebo", "FMT"))
    ) 
  

  if (!is.null(subset_expr)) {
    keep <- eval(parse(text = subset_expr), envir = meta_sub)
    mat_sub  <- mat_sub[keep, ]
    meta_sub <- meta_sub[keep, ]
  }

  if (is.null(output_dir))
    output_dir <- tempdir()
  
  dir.create(output_dir)
  fit <- Maaslin2(
    input_data     = mat_sub,
    input_metadata = meta_sub,
    output         = output_dir,
    fixed_effects  = fixed_effects,
    random_effects = random_effects,
    reference      = reference,
    normalization  = "TSS",      # 已經是相對豐度
    transform      = "LOG",       # log 轉換處理偏態
    analysis_method = "LM",
    min_prevalence = 0.1,
    min_abundance  = 0.0001,
    max_significance = 0.25,
    plot_heatmap   = FALSE,
    plot_scatter   = FALSE
  )

  res <- fit$results %>%
    as_tibble() %>%
    rename(Pathway = feature) %>%
    mutate(label = label) %>%
    arrange(pval)

  cat(sprintf("[%s] pval<0.05: %d  pval<0.1: %d  q<0.25: %d\n", label,
              sum(res$pval < 0.05, na.rm = TRUE),
              sum(res$pval < 0.1,  na.rm = TRUE),
              sum(res$qval < 0.25, na.rm = TRUE)))

  return(res)
}

# ── 2a. 全模型：Treatment + Timepoint（所有樣本）────────────────────────────

dir_mc_full <- file.path(getwd(), "outputs/pathway/maaslin2_metacyc_full")
dir_ec_full <- file.path(getwd(), "outputs/pathway/maaslin2_ec_full")

maas_metacyc_full <- run_maaslin2(metacyc, meta_raw, "MetaCyc",
                                    output_dir = dir_mc_full)
maas_ec_full      <- run_maaslin2(ec,      meta_raw, "EC",
                                    output_dir = dir_ec_full)

# ── 2b. FMT 組內：Baseline vs 6M（paired）────────────────────────────────────

meta_fmt <- meta_raw %>%
  filter(Treatment == "FMT",
         Timepoint %in% c("Baseline", "6M post-FMT"),
         PatientID != "GM03")

dir_mc_fmt <- file.path(getwd(), "outputs/pathway/maaslin2_metacyc_FMT")

maas_metacyc_fmt <- run_maaslin2(
  metacyc, meta_fmt, "MetaCyc_FMT",
  fixed_effects  = c("Timepoint"),
  random_effects = c("PatientID"),
  reference      = c("Timepoint,Baseline"),
  output_dir     = dir_mc_fmt
)

# ── 2c. Placebo 組內：Baseline vs 6M（對照）──────────────────────────────────

meta_plc <- meta_raw %>%
  filter(Treatment == "placebo",
         Timepoint %in% c("Baseline", "6M post-FMT"))

dir_mc_plc <- file.path(getwd(), "outputs/pathway/maaslin2_metacyc_placebo")

maas_metacyc_plc <- run_maaslin2(
  metacyc, meta_plc, "MetaCyc_placebo",
  fixed_effects  = c("Timepoint"),
  random_effects = c("PatientID"),
  reference      = c("Timepoint,Baseline"),
  output_dir     = dir_mc_plc
)

# ── 2d. FMT-specific pathways（A ∩ B, not C）─────────────────────────────────

fmt_pw_A <- maas_metacyc_full %>%
  filter(metadata == "Treatment", coef > 0, pval < 0.05) %>%
  pull(Pathway)

fmt_pw_B <- maas_metacyc_fmt %>%
  filter(metadata == "Timepoint", coef > 0, pval < 0.05) %>%
  pull(Pathway)

fmt_pw_C <- maas_metacyc_plc %>%
  filter(metadata == "Timepoint", pval < 0.05) %>%
  pull(Pathway)

fmt_specific_pw <- intersect(fmt_pw_A, fmt_pw_B) %>% # 沒有任何重疊的pathway....
  setdiff(fmt_pw_C)

cat(sprintf("\nFMT-specific pathways (A∩B not C): %d\n",
            length(fmt_specific_pw)))
print(fmt_specific_pw)


# =============================================================================
#  SECTION 3 — Differential pathway 視覺化
# =============================================================================

cat("\n── Section 3: Visualisation ─────────────────────────────────────────────\n")

# ── 3a. Volcano plot（MaAsLin2 結果）─────────────────────────────────────────

plot_maaslin_volcano <- function(res_df, coef_var = "Treatment",
                                  title = "") {
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
    geom_hline(yintercept = -log10(0.25), linetype = "dashed",
               colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "#D85A30", linewidth = 0.4) +
    scale_colour_manual(
      values = c("Up" = "#D6604D", "Down" = "#2166AC", NS = "grey70")
    ) +
    labs(title    = title,
         subtitle = sprintf("q<0.25: %d up, %d down",
                            sum(df$direction == "Up"),
                            sum(df$direction == "Down")),
         x = "Coefficient (log scale)", y = "-log10(q-value)",
         colour = NULL) +
    theme_fmt
}

p_vol_metacyc <- plot_maaslin_volcano(maas_metacyc_full,
                                       title = "MetaCyc: FMT vs Placebo (all timepoints)")
p_vol_metacyc_fmt <- plot_maaslin_volcano(maas_metacyc_fmt,
                                           coef_var = "Timepoint",
                                           title = "MetaCyc: FMT Baseline → 6M")

# ── 3b. Top pathway trajectory plot ──────────────────────────────────────────

plot_pw_trajectory <- function(mat, meta, pathways, label) {
  pw_long <- as.data.frame(t(mat[pathways, , drop = FALSE])) %>%
    tibble::rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "Pathway", values_to = "abund") %>%
    left_join(meta %>% dplyr::select(SampleID, PatientID, Timepoint, Treatment),
              by = "SampleID") %>%
    filter(!is.na(Treatment)) %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_labels)) %>%
    group_by(Pathway, Treatment, Timepoint) %>%
    summarise(mean_abund = mean(abund),
              se         = sd(abund) / sqrt(n()),
              .groups = "drop")

  ggplot(pw_long, aes(x = Timepoint, y = mean_abund,
                       colour = Treatment, group = Treatment)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_abund - se, ymax = mean_abund + se),
                  width = 0.2) +
    facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = col_treatment) +
    labs(title    = sprintf("[%s] FMT-specific pathway trajectories", label),
         subtitle = "Mean ± SE",
         x = NULL, y = "Pathway abundance") +
    theme_fmt +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          strip.text  = element_text(size = 7),
          legend.position = "bottom")
}

if (length(fmt_specific_pw) > 0) {
  n_show <- min(12, length(fmt_specific_pw))
  p_pw_traj <- plot_pw_trajectory(
    metacyc, meta_raw,
    pathways = fmt_specific_pw[1:n_show],
    label    = "MetaCyc"
  )
}


# =============================================================================
#  SECTION 4 — Pathway × SALT correlation
# =============================================================================

cat("\n── Section 4: Pathway × SALT score correlation ──────────────────────────\n")

pw_salt_cor <- function(mat, meta, label, end_tp = "6M post-FMT") {

  # 計算 pathway 變化量（end_tp - Baseline）
  meta_wide <- meta %>%
    filter(Timepoint %in% c("Baseline", end_tp)) %>%
    dplyr::select(SampleID, PatientID, Timepoint, Treatment, SALT_score)

  common <- intersect(colnames(mat), meta_wide$SampleID)
  mat_sub <- mat[, common]

  delta_pw <- meta_wide %>%
    dplyr::select(SampleID, PatientID, Timepoint) %>%
    filter(SampleID %in% common) %>%
    pivot_wider(names_from = Timepoint, values_from = SampleID) %>%
    rename(sid_bl = Baseline, sid_et = !!sym(end_tp)) %>%
    filter(!is.na(sid_bl) & !is.na(sid_et)) %>%
    mutate(delta = purrr::map2(sid_et, sid_bl, function(et, bl) {
      mat_sub[, et] - mat_sub[, bl]
    })) %>%
    dplyr::select(PatientID, delta) %>%
    mutate(delta = purrr::map(delta, ~ setNames(.x, rownames(mat_sub))))

  delta_mat <- do.call(cbind, delta_pw$delta)
  colnames(delta_mat) <- delta_pw$PatientID

  # SALT 變化量
  salt_delta <- meta_wide %>%
    dplyr::select(PatientID, Timepoint, SALT_score) %>%
    distinct() %>%
    pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
    rename(salt_bl = Baseline, salt_et = !!sym(end_tp)) %>%
    mutate(dSALT = salt_et - salt_bl) %>%
    filter(PatientID %in% colnames(delta_mat))

  # Spearman correlation：每個 pathway vs dSALT
  cor_res <- apply(delta_mat[, salt_delta$PatientID], 1, function(pw_delta) {
    ct <- cor.test(pw_delta, salt_delta$dSALT,
                   method = "spearman", exact = FALSE)
    c(rho = ct$estimate, pval = ct$p.value)
  }) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Pathway") %>%
    mutate(
      padj = p.adjust(pval, method = "BH"),
      sig  = case_when(padj < 0.05 ~ "*", padj < 0.1 ~ "†", TRUE ~ "")
    ) %>%
    rename_with(~"rho", contains("rho")) %>% 
    arrange(padj)

  cat(sprintf("[%s] Pathways correlated with ΔSALT (BH-adj p<0.05): %d\n",
              label, sum(cor_res$padj < 0.05)))
  cat(sprintf("[%s] Pathways correlated with ΔSALT (BH-adj p<0.1):  %d\n",
              label, sum(cor_res$padj < 0.1)))
  
  # Top 20 plot
  top_cor <- cor_res %>% slice_max(abs(rho), n = 20)

  p_cor <- ggplot(top_cor, aes(x = reorder(Pathway, rho), y = rho,
                                fill = rho > 0)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig),
              hjust = if_else(top_cor$rho > 0, -0.2, 1.2),
              size  = 3.5) +
    geom_hline(yintercept = 0, colour = "grey40") +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "#2166AC")) +
    scale_y_continuous(limits = c(-1.1, 1.1)) +
    labs(title    = sprintf("[%s] Pathway Δabundance vs ΔSALT (Spearman ρ)", label),
         subtitle = "Negative ΔSALT = clinical improvement  |  * BH-adj p<0.05",
         x = NULL, y = "Spearman ρ") +
    theme_fmt +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8))

  list(cor_res = cor_res, p_cor = p_cor, delta_mat = delta_mat)
}

salt_metacyc <- pw_salt_cor(metacyc, meta_raw, "MetaCyc")
salt_ec      <- pw_salt_cor(ec,      meta_raw, "EC")

write.csv(salt_metacyc$cor_res, "outputs/pathway/pathway_salt_correlation_metacyc.csv", row.names = FALSE)
write.csv(salt_ec$cor_res,      "outputs/pathway/pathway_salt_correlation_ec.csv",      row.names = FALSE)


# =============================================================================
#  SECTION 5 — Taxa–Pathway 連結（與 DESeq2 結果對照）
# =============================================================================

cat("\n── Section 5: Taxa–Pathway overlap with DESeq2 ─────────────────────────\n")

# 找出 DESeq2 FMT-specific taxa 和 MaAsLin2 FMT-specific pathway
# 都往同一方向（Up）且都與 SALT 改善相關的

if (exists("fmt_specific") && is.data.frame(fmt_specific) && nrow(fmt_specific) > 0) {
  cat("DESeq2 FMT-specific taxa:\n")
  print(fmt_specific$Genus)
}

cat("\nMaAsLin2 FMT-specific pathways:\n")
if (exists("fmt_specific_pw") && is.data.frame(fmt_specific_pw)) print(fmt_specific_pw)

# Pathway 和 SALT correlation 的交集
pw_salt_sig <- salt_metacyc$cor_res %>%
  filter(padj < 0.1, rho < 0)   # 負 rho = pathway 增加對應 SALT 下降（改善）

pw_fmt_and_salt <- intersect(fmt_specific_pw, pw_salt_sig$Pathway)

cat(sprintf("\nPathways that are FMT-specific AND negatively correlated with ΔSALT: %d\n",
            length(pw_fmt_and_salt)))
print(pw_fmt_and_salt)

write.csv(
  data.frame(Pathway = pw_fmt_and_salt),
  "outputs/pathway/pathway_FMT_specific_and_SALT_correlated.csv",
  row.names = FALSE
)


# =============================================================================
#  Save
# =============================================================================

save_fig <- function(plot, filename, w = 10, h = 6) {
  if (is.null(plot)) { message("Skipping NULL: ", filename); return(invisible()) }
  ggsave(filename, plot = plot, width = w, height = h, dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", filename))
}

dir.create("outputs/pathway", recursive = T)
# Beta diversity
save_fig(beta_metacyc$p_bray,  "outputs/pathway/pw_fig1a_metacyc_bray_pcoa.pdf")
save_fig(beta_metacyc$p_ait,   "outputs/pathway/pw_fig1b_metacyc_ait_pcoa.pdf")
save_fig(beta_metacyc$p_arrow, "outputs/pathway/pw_fig1c_metacyc_shift_envfit.pdf", w = 9, h = 7)
save_fig(beta_ec$p_bray,       "outputs/pathway/pw_fig1d_ec_bray_pcoa.pdf")
save_fig(beta_ec$p_arrow,      "outputs/pathway/pw_fig1e_ec_shift_envfit.pdf",      w = 9, h = 7)

# Comparison panel
p_beta_compare <- (beta_metacyc$p_bray | beta_ec$p_bray) +
  plot_annotation(title = "Pathway beta diversity: MetaCyc vs EC",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
save_fig(p_beta_compare, "outputs/pathway/pw_fig1_beta_comparison.pdf", w = 14, h = 6)

p_beta_compare1 <- (beta_metacyc$p_ait | beta_ec$p_ait) +
  plot_annotation(title = "Pathway beta diversity: MetaCyc vs EC",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
save_fig(p_beta_compare1, "outputs/pathway/pw_fig1_beta_comparison_ait.pdf", w = 14, h = 6)

# Differential
save_fig(p_vol_metacyc,     "outputs/pathway/pw_fig2a_volcano_metacyc_treatment.pdf")
save_fig(p_vol_metacyc_fmt, "outputs/pathway/pw_fig2b_volcano_metacyc_FMT_temporal.pdf")
if (exists("p_pw_traj"))
  save_fig(p_pw_traj, "outputs/pathway/pw_fig3_FMT_specific_trajectories.pdf", w = 14, h = 10)

# SALT correlation
save_fig(salt_metacyc$p_cor, "outputs/pathway/pw_fig4a_metacyc_salt_correlation.pdf", w = 10, h = 8)
save_fig(salt_ec$p_cor,      "outputs/pathway/pw_fig4b_ec_salt_correlation.pdf",      w = 10, h = 8)

cat("\n✓ Pathway analysis complete.\n")
cat("  Key results:\n")
cat("    fmt_specific_pw         — FMT-specific pathways\n")
cat("    pw_fmt_and_salt         — FMT-specific + SALT-correlated\n")
cat("    salt_metacyc$cor_res    — All MetaCyc × SALT correlations\n")
cat("    beta_metacyc$fit_scores — envfit significant pathways\n")
