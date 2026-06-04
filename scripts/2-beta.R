# =============================================================================
#  beta_diversity.R
#  Beta diversity analysis — Bray-Curtis vs Aitchison
#
#  主要 function：run_beta_diversity(dist_type = "bray" | "aitchison")
#  最後呼叫兩次並產出並排比較圖
#
#  前置條件：
#    - ps_filt, meta_raw, col_treatment, col_timepoint, tp_labels, theme_fmt
#      已在主 script 定義
#    - 需要 vegan, phyloseq, ggplot2, ggrepel, patchwork, tidyverse
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")


# =============================================================================
#  核心 function
# =============================================================================

run_beta_diversity <- function(
    dist_type = c("bray", "aitchison"),  # 距離算法
    ps_obj    = ps_filt,                 # phyloseq object（未 rarefy 原始 counts）
    meta      = meta_raw,
    end_tp    = "6M post-FMT",           # 位移分析的終點時間
    n_envfit  = 10,                      # envfit top N genera
    seed      = 42
) {
  
  dist_type <- match.arg(dist_type)
  cat(sprintf("\n%s\n── Running beta diversity: %s ──\n%s\n",
              strrep("=", 60), toupper(dist_type), strrep("=", 60)))
  
  # ── Step 1：建立距離矩陣 ───────────────────────────────────────────────────
  
  if (dist_type == "bray") {
    
    dist_mat    <- phyloseq::distance(ps_obj, method = "bray")
    ps_for_pcoa <- ps_obj
    dist_label  <- "Bray-Curtis"
    
  } else {
    
    # Aitchison：agglomerate to genus → pseudo-count → CLR → euclidean
    ps_for_pcoa <- tax_glom(ps_obj, taxrank = "Genus") %>%
      transform_sample_counts(function(x) x + 0.5) %>%
      transform_sample_counts(function(x) log(x) - mean(log(x)))
    
    clr_mat  <- t(as.matrix(otu_table(ps_for_pcoa)))
    dist_mat <- dist(clr_mat, method = "euclidean")
    dist_label <- "Aitchison (CLR)"
    
  }
  
  # ── Step 2：PCoA ──────────────────────────────────────────────────────────
  
  ord     <- ordinate(ps_for_pcoa, method = "PCoA", distance = dist_mat)
  var_exp <- ord$values$Relative_eig * 100
  
  cat(sprintf("[%s] PC1: %.1f%%  PC2: %.1f%%  PC3: %.1f%%\n",
              dist_label, var_exp[1], var_exp[2], var_exp[3]))
  
  scores_raw <- data.frame(
    ord$vectors[, 1:3],
    as.data.frame(sample_data(ps_for_pcoa))
  ) %>%
    rename(PC1 = Axis.1, PC2 = Axis.2, PC3 = Axis.3)
  
  # SampleID 可能已存在於 sample_data，先檢查再決定是否從 rownames 補
  if (!"SampleID" %in% names(scores_raw)) {
    scores_raw <- tibble::rownames_to_column(scores_raw, "SampleID")
  }
  
  scores <- scores_raw %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_labels)) %>%
    arrange(PatientID, Timepoint)
  
  xlab <- sprintf("PC1 (%.1f%%)", var_exp[1])
  ylab <- sprintf("PC2 (%.1f%%)", var_exp[2])
  
  # ── Step 3：PERMANOVA（with blocks）──────────────────────────────────────
  
  meta_ord  <- meta %>%
    filter(SampleID %in% sample_names(ps_for_pcoa)) %>%
    arrange(match(SampleID, sample_names(ps_for_pcoa)))
  perm_ctrl <- permute::how(blocks = meta_ord$PatientID, nperm = 999)
  
  set.seed(seed)
  perm_res <- adonis2(
    dist_mat ~ Treatment * Timepoint + Age,
    data         = meta_ord,
    permutations = perm_ctrl,
    by           = "terms"
  )
  cat(sprintf("\n[%s] PERMANOVA (blocks = PatientID):\n", dist_label))
  print(perm_res)
  
  # betadisper
  bd      <- betadisper(dist_mat, meta_ord$Treatment)
  bd_test <- permutest(bd, permutations = 999)
  cat(sprintf("\n[%s] betadisper (Treatment):\n", dist_label))
  print(bd_test)
  
  # ── Step 4：Within-subject distance change ────────────────────────────────
  
  bc_mat    <- as.matrix(dist_mat)
  change_df <- meta_ord %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_labels)) %>%
    group_by(PatientID) %>%
    group_modify(~ {
      sids <- .x$SampleID
      tp   <- .x$Timepoint
      expand.grid(from = sids, to = sids, stringsAsFactors = FALSE) %>%
        filter(from != to) %>%
        mutate(
          bc      = bc_mat[cbind(from, to)],
          from_tp = tp[match(from, sids)],
          to_tp   = tp[match(to,   sids)]
        )
    }) %>%
    filter(as.integer(from_tp) < as.integer(to_tp)) %>%
    left_join(meta_ord %>% select(SampleID, Treatment),
              by = c("from" = "SampleID")) %>%
    ungroup()
  
  # ── Step 5：Directional shift (Baseline → end_tp) ────────────────────────
  incomplete_patients <- scores %>%
    filter(Timepoint %in% c("Baseline", end_tp)) %>%
    group_by(PatientID) %>%
    filter(n() < 2) %>%
    pull(PatientID) %>%
    unique()
  
  if (length(incomplete_patients) > 0)
    message(sprintf("[%s] 以下病人缺少 Baseline 或 %s，排除於 shift 分析：%s",
                    dist_label, end_tp,
                    paste(incomplete_patients, collapse = ", ")))
  
  shift_wide <- scores %>%
    filter(Timepoint %in% c("Baseline", end_tp)) %>%
    select(PatientID, Treatment, Timepoint, PC1, PC2) %>%
    pivot_wider(names_from  = Timepoint,
                values_from = c(PC1, PC2))
  
  # 動態 rename：欄位名稱含空格，用 backtick
  shift_df <- shift_wide %>%
    rename(
      PC1_bl = `PC1_Baseline`,
      PC1_et = !!sym(paste0("PC1_", end_tp)),
      PC2_bl = `PC2_Baseline`,
      PC2_et = !!sym(paste0("PC2_", end_tp))
    ) %>%
    mutate(dPC1 = PC1_et - PC1_bl,
           dPC2 = PC2_et - PC2_bl) %>%
    filter(!is.na(dPC1) & !is.na(dPC2))   # ← 加這行，排除缺 baseline 的 GM03
  
  # 方向一致性
  cat(sprintf("\n[%s] PC1 shift (Baseline → %s):\n", dist_label, end_tp))
  shift_df %>%
    group_by(Treatment) %>%
    summarise(
      mean_dPC1          = round(mean(dPC1), 4),
      sd_dPC1            = round(sd(dPC1),   4),
      n_positive         = sum(dPC1 > 0),
      n_negative         = sum(dPC1 < 0),
      all_same_direction = all(dPC1 > 0) | all(dPC1 < 0),
      .groups = "drop"
    ) %>% print()
  
  # Wilcoxon dPC1 FMT vs placebo
  wx_dpc1 <- wilcox.test(dPC1 ~ Treatment, data = shift_df, exact = FALSE) #動的方向（Direction）
  cat(sprintf("[%s] Wilcoxon dPC1: W=%.1f, p=%.4f\n",
              dist_label, wx_dpc1$statistic, wx_dpc1$p.value))
  
  # Wilcoxon on distance (Baseline → end_tp)
  bc_end  <- change_df %>% filter(from_tp == "Baseline", to_tp == end_tp)
  wx_fmt  <- wilcox.test(bc_end %>% filter(Treatment == "FMT")     %>% pull(bc),
                         mu = 0, alternative = "greater", exact = FALSE)
  wx_plc  <- wilcox.test(bc_end %>% filter(Treatment == "placebo") %>% pull(bc),
                         mu = 0, alternative = "greater", exact = FALSE)
  wx_grp  <- wilcox.test(bc ~ Treatment, data = bc_end, exact = FALSE)# important:在不同treatment下移動的距離有沒有差異?
  
  cat(sprintf("[%s] dist>0 FMT:     p=%.4f\n", dist_label, wx_fmt$p.value)) # 這群人的位移（Distance）平均而言是不是 0 / 動了多遠（Magnitude）
  cat(sprintf("[%s] dist>0 placebo: p=%.4f\n", dist_label, wx_plc$p.value))
  cat(sprintf("[%s] FMT vs placebo: p=%.4f\n", dist_label, wx_grp$p.value))
  
  # ── Step 6：SALT change + Spearman ───────────────────────────────────────
  
  salt_change <- meta %>%
    filter(Timepoint %in% c("Baseline", end_tp)) %>%
    select(PatientID, Timepoint, SALT_score) %>%
    pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
    rename(SALT_bl = Baseline, SALT_et = !!sym(end_tp)) %>%
    mutate(dSALT = SALT_et - SALT_bl)
  
  shift_salt <- left_join(shift_df, salt_change, by = "PatientID")
  
  cor_shift <- cor.test(shift_salt$dPC1, shift_salt$dSALT,
                        method = "spearman", exact = FALSE)
  cat(sprintf("[%s] Spearman dPC1~dSALT: rho=%.3f, p=%.4f\n",
              dist_label, cor_shift$estimate, cor_shift$p.value)) # 菌相沿著 PC1 軸移動的距離（dPC1），是否能預測脫髮改善的程度（dSALT）？
  
  # ── Step 7：Plots ─────────────────────────────────────────────────────────
  
  # 7a. PCoA by treatment（trajectory）
  p_pcoa_trt <- ggplot(scores, aes(x = PC1, y = PC2, colour = Treatment)) +
    geom_path(aes(group = PatientID),
              colour = "gray40", linewidth = 0.4, linetype = "dashed",
              arrow  = arrow(length = unit(0.2, "cm"), type = "closed"),
              alpha  = 0.6) +
    geom_point(aes(shape = Timepoint), size = 3.5, alpha = 0.85) +
    stat_ellipse(aes(group = Treatment), level = 0.75, linewidth = 0.8) +
    geom_text(data    = filter(scores, Timepoint == "Baseline"),
              aes(label = PatientID),
              nudge_x = 0.01, nudge_y = -0.01, size = 2.8) +
    scale_colour_manual(values = col_treatment) +
    scale_shape_manual(values  = c(16, 17, 15)) +
    labs(title    = sprintf("PCoA (%s) — Treatment trajectory", dist_label),
         subtitle = sprintf("PERMANOVA Treatment*Timepoint p = %.3f",
                            perm_res["Treatment:Timepoint", "Pr(>F)"]),
         x = xlab, y = ylab,
         colour = "Treatment", shape = "Time point") +
    theme(aspect.ratio = 1)+
    theme_fmt
  
  # 7b. PCoA by timepoint
  p_pcoa_tp <- ggplot(scores, aes(x = PC1, y = PC2, colour = Timepoint)) +
    geom_path(aes(group = PatientID),
              colour = "gray40", linewidth = 0.4, linetype = "dashed",
              arrow  = arrow(length = unit(0.2, "cm"), type = "closed"),
              alpha  = 0.6) +
    geom_point(aes(shape = Treatment), size = 3.5, alpha = 0.85) +
    stat_ellipse(aes(group = Timepoint), level = 0.75, linewidth = 0.8) +
    geom_text(data    = filter(scores, Timepoint == "Baseline"),
              aes(label = PatientID),
              nudge_x = 0.01, nudge_y = -0.01, size = 2.8) +
    scale_colour_manual(values = col_timepoint) +
    scale_shape_manual(values  = c(16, 17)) +
    labs(title  = sprintf("PCoA (%s) — Time point", dist_label),
         x = xlab, y = ylab,
         colour = "Time point", shape = "Treatment") +
    theme(aspect.ratio = 1)+
    theme_fmt
  
  # 7c. Within-subject distance change boxplot
  p_bc_change <- ggplot(change_df,
                        aes(x = paste(from_tp, "→", to_tp),
                            y = bc, fill = Treatment)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_point(
      aes(colour = Treatment, shape = Treatment),
      position = position_jitterdodge(jitter.width = 0.3,
                                      dodge.width  = 0.75, seed = seed),
      size = 3, alpha = 0.6
    ) +
    geom_text_repel(
      aes(label = PatientID, colour = Treatment),
      position       = position_jitterdodge(jitter.width = 0.3,
                                            dodge.width  = 0.75, seed = seed),
      segment.colour = NA, box.padding   = 0.2,
      point.padding  = 0.3, max.overlaps = Inf,
      size = 3, show.legend = FALSE, fontface = "bold"
    ) +
    scale_fill_manual(values   = col_treatment) +
    scale_colour_manual(values = col_treatment) +
    scale_shape_manual(values  = c(FMT = 16, placebo = 17)) +
    labs(title = sprintf("Within-subject %s change", dist_label),
         x = "Transition", y = "Distance") +
    theme_fmt
  
  # 7d. Arrow plot（Baseline → end_tp）
  p_arrow <- ggplot() +
    geom_point(
      data = filter(scores, Timepoint == "Baseline"),
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
    labs(title    = sprintf("PCoA (%s) displacement", dist_label),
         subtitle = sprintf("Circle = Baseline | Arrowhead = %s", end_tp),
         x = xlab, y = ylab) +
    theme(aspect.ratio = 1)+
    theme_fmt
  
  # 7e. dPC1 bar per patient
  p_dpc1 <- ggplot(shift_df,
                   aes(x = reorder(PatientID, dPC1), y = dPC1,
                       fill = Treatment)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = col_treatment) +
    labs(title = sprintf("[%s] ΔPC1 per patient", dist_label),
         x = NULL, y = "ΔPC1") +
    theme_fmt +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  
  # 7f. dPC1 vs dSALT scatter
  p_salt <- ggplot(shift_salt,
                   aes(x = dPC1, y = dSALT, colour = Treatment)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "grey40", fill = "grey85", linewidth = 0.8) +
    geom_text_repel(aes(label = PatientID), size = 3, segment.colour = NA) +
    scale_colour_manual(values = col_treatment) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("ρ = %.3f\np = %.3f",
                             cor_shift$estimate, cor_shift$p.value),
             hjust = 1.1, vjust = 1.5, size = 3.5, colour = "grey30") +
    labs(title    = sprintf("[%s] ΔPC1 vs ΔSALT", dist_label),
         subtitle = "Negative dSALT = clinical improvement",
         x = "ΔPC1", y = "ΔSALT") +
    theme_fmt
  
  # ── Step 8a：envfit ────────────────────────────────────────────────────────
  
  # Aitchison 的 envfit 要用 CLR genus matrix，不是 rel abund
  if (dist_type == "aitchison") {
    # 直接用已經建好的 CLR matrix（ps_for_pcoa 裡就是 CLR）
    genus_mat_env <- as.data.frame(t(otu_table(ps_for_pcoa)))
    colnames(genus_mat_env) <- as.data.frame(tax_table(ps_for_pcoa))$Genus
  } else {
    # Bray-Curtis 用 relative abundance
    ps_genus_rel_env <- tax_glom(ps_obj, taxrank = "Genus") %>%
      transform_sample_counts(function(x) x / sum(x))
    genus_mat_env <- as.data.frame(t(otu_table(ps_genus_rel_env)))
    colnames(genus_mat_env) <- as.data.frame(tax_table(ps_genus_rel_env))$Genus
  }
  
  samples_sub   <- scores %>%
    filter(Timepoint %in% c("Baseline", end_tp)) %>%
    pull(SampleID)
  pcoa_mat      <- scores %>%
    filter(SampleID %in% samples_sub) %>%
    select(PC1, PC2) %>%
    as.matrix()
  
  genus_mat_sub <- genus_mat_env[rownames(pcoa_mat), ]
  
  set.seed(seed)
  fit <- envfit(pcoa_mat, genus_mat_sub, permutations = 999, na.rm = TRUE)
  
  fit_scores <- as.data.frame(scores(fit, display = "vectors")) %>%
    tibble::rownames_to_column("Genus") %>%
    mutate(r2   = fit$vectors$r,
           pval = fit$vectors$pvals) %>%
    filter(pval < 0.05) %>%
    arrange(desc(r2))
  
  cat(sprintf("[%s] envfit significant genera (p<0.05): %d\n",
              dist_label, nrow(fit_scores)))
  
  top_taxa     <- slice_max(fit_scores, r2, n = n_envfit)
  pcoa_range  <- max(abs(c(range(scores$PC1), range(scores$PC2))))
  vector_range <- max(sqrt(top_taxa$PC1^2 + top_taxa$PC2^2))
  scale_factor <- (pcoa_range * 0.35) / vector_range
  
  cat(sprintf("[%s] auto scale_factor: %.4f\n", dist_label, scale_factor))
  
  p_envfit <- p_arrow +
    geom_segment(
      data = top_taxa,
      aes(x = 0, y = 0,
          xend = PC1 * scale_factor,
          yend = PC2 * scale_factor),
      arrow       = arrow(length = unit(0.012, "npc"), type = "closed"),
      colour      = "grey20", linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    geom_text_repel(
      data = top_taxa,
      aes(x = PC1 * scale_factor,
          y = PC2 * scale_factor,
          label = Genus),
      colour = "grey10", size = 3, fontface = "italic",
      segment.colour = NA, inherit.aes = FALSE
    ) +
    theme(aspect.ratio = 1) +
    labs(subtitle = sprintf("envfit p<0.05 top %d genera | BL=circle, %s=head",
                            n_envfit, end_tp))
  
  # ── Step 8b：Genus × SALT correlation heatmap ────────────────────────────
  
  cat(sprintf("\n[%s] Genus–SALT correlation heatmap\n", dist_label))
  
  top_genera_names <- c(top_taxa$Genus, "Bifidobacterium")
  
  genus_change <- genus_mat_env %>%
    tibble::rownames_to_column("SampleID") %>%
    left_join(
      meta %>% select(SampleID, PatientID, Timepoint),
      by = "SampleID"
    ) %>%
    filter(Timepoint %in% c("Baseline", end_tp)) %>%
    select(PatientID, Timepoint, all_of(top_genera_names)) %>%
    pivot_longer(cols = all_of(top_genera_names),
                 names_to = "Genus", values_to = "abund") %>%
    pivot_wider(names_from = Timepoint, values_from = abund) %>%
    rename(bl = Baseline, m_end = !!sym(end_tp)) %>%
    mutate(delta_abund = m_end - bl) %>%
    select(PatientID, Genus, delta_abund) %>%
    pivot_wider(names_from = Genus, values_from = delta_abund)
  
  cor_input <- genus_change %>%
    left_join(salt_change %>% select(PatientID, dSALT), by = "PatientID")
  
  # Spearman per treatment group（分組，exploratory）
  cor_results <- c("FMT", "placebo") %>%
    purrr::map_dfr(function(trt) {
      patients_in_group <- meta %>%
        filter(Timepoint == "Baseline", Treatment == trt) %>%
        pull(PatientID)
      
      sub_input <- cor_input %>%
        filter(PatientID %in% patients_in_group)
      
      top_genera_names %>%
        purrr::map_dfr(function(g) {
          test <- cor.test(sub_input[[g]], sub_input$dSALT,
                           method = "spearman", exact = FALSE)
          tibble(
            Treatment = trt,
            Genus     = g,
            rho       = round(test$estimate, 3),
            pval      = round(test$p.value,  3),
            sig       = case_when(
              test$p.value < 0.05 ~ "*",
              test$p.value < 0.1  ~ "†",
              TRUE                ~ ""
            )
          )
        })
    })
  
  cat(sprintf("[%s] Spearman correlation by treatment group (exploratory):\n", dist_label))
  print(cor_results)
  
  # Heatmap matrix（z-score）
  heatmap_mat <- cor_input %>%
    select(PatientID, all_of(top_genera_names), dSALT) %>%
    mutate(
      across(all_of(top_genera_names), ~ scale(.)[, 1]),
      dSALT = scale(dSALT)[, 1]
    ) %>%
    left_join(
      meta %>% filter(Timepoint == "Baseline") %>%
        select(PatientID, Treatment) %>% distinct(),
      by = "PatientID"
    ) %>%
    filter(!is.na(Treatment))   # 排除無 baseline 的樣本（如 GM03）
  
  genus_order <- cor_results %>%
    filter(Treatment == "FMT") %>%
    arrange(rho) %>%
    pull(Genus)
  
  heatmap_long <- heatmap_mat %>%
    pivot_longer(cols = c(all_of(top_genera_names), "dSALT"),
                 names_to = "Variable", values_to = "z") %>%
    left_join(
      cor_results %>% select(Treatment, Genus, sig),
      by = c("Treatment", "Variable" = "Genus")
    ) %>%
    mutate(
      sig      = replace_na(sig, ""),
      var_type = if_else(Variable == "dSALT", "Clinical", "Genus"),
      Variable = factor(Variable, levels = c("dSALT", genus_order))
    )
  
  p_cor_heatmap <- ggplot(heatmap_long,
                          aes(x = Variable, y = PatientID, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig), size = 4, colour = "gray10", fontface = "bold") +
    # geom_vline(xintercept = 1.5, colour = "grey30", linewidth = 1) +
    facet_grid(Treatment ~ var_type, scales = "free", space = "free") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = 0, name = "z-score"
    ) +
    labs(
      title = sprintf("[%s] Δ %s vs Δ SALT (%s − Baseline)",
                      dist_label,
                      if_else(dist_type == "aitchison", "CLR abundance", "Rel. abundance"),
                      end_tp),
      subtitle = "* p<0.05  † p<0.1  |  Columns sorted by Spearman ρ (FMT group)",
      x = NULL, y = NULL
    ) +
    theme_fmt +
    theme(
      axis.text.x     = element_text(angle = 40, hjust = 1, face = "italic"),
      strip.text.x    = element_text(face = "bold"),
      legend.position = "right"
    )
  
  p_rho_bar <- ggplot(cor_results,
                      aes(x = reorder(Genus, rho), y = rho, fill = rho > 0)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(sprintf("%.2f", rho), sig)),
              hjust = if_else(cor_results$rho > 0, -0.1, 1.1),
              size  = 2.8) +
    geom_hline(yintercept = 0, colour = "grey40") +
    coord_flip() +
    facet_wrap(~ Treatment, ncol = 1) +
    scale_fill_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "#2166AC")) +
    scale_y_continuous(limits = c(-1.2, 1.2)) +
    labs(title = sprintf("[%s] Spearman ρ\n(Δabundance vs ΔSALT)", dist_label),
         x = NULL, y = "ρ") +
    theme_fmt +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "italic", size = 9))
  
  p_cor_panel <- p_cor_heatmap + p_rho_bar +
    plot_layout(widths = c(3, 1)) +
    plot_annotation(
      title = sprintf("Genus–SALT correlation panel [%s]", dist_label),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  # 分開畫每個bac變動下salt score的改變
  tmp <- meta_raw %>%
    select(PatientID, Treatment) %>%
    distinct() %>%
    left_join(cor_input, by = "PatientID") %>%
    filter(!is.na(dSALT)) %>%
    filter(if_all(all_of(top_genera_names), ~ !is.na(.x)))
  
  plot_list <- list() # 建立一個空的 list 來存圖
  
  for (i in top_genera_names) {
    y_lab <- paste0("Δ", i)
    
    ct <- cor.test(tmp[[i]], tmp$dSALT, method = "spearman", exact = FALSE)
    
    # p <- tmp %>% 
    #   ggplot(aes(x = dSALT, y = !!sym(i), color = Treatment)) +
    #   geom_point(size = 3, alpha = 0.6) +
    #   geom_smooth(method = "lm", se = FALSE) + 
    #   theme_bw() +
    #   theme(legend.position = "none") + # 先關掉 legend，最後再統一顯示
    #   labs(y = y_lab, x = "ΔSALT")+
    #   scale_color_manual(values = col_treatment)
    
    p <- ggplot(tmp, aes(x = dSALT, y = .data[[i]], colour = Treatment)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_point(size = 3, alpha = 0.75) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("ρ=%.2f\np=%.3f", ct$estimate, ct$p.value),
               hjust = 1.1, vjust = 1.4, size = 2.8, colour = "grey30") +
      scale_colour_manual(values = col_treatment) +
      labs(y = bquote(Delta ~ italic(.(i))),
           x = "ΔSALT (lower = better)") +
      theme_bw() +
      theme(legend.position = "none",
            axis.title      = element_text(size = 9),
            plot.title      = element_blank())

    
    plot_list[[i]] <- p # 將圖存入 list
  }
  
  
  combined_plot <- wrap_plots(plot_list, ncol = 5) + 
    plot_layout(guides = "collect") + # 統一收集 legend
    plot_annotation(
      title    = "Δ Genus abundance vs Δ SALT score",
      subtitle = "Negative dSALT = clinical improvement",
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))
    ) &
    theme(legend.position = "bottom")
  
  # ── Step 9：回傳所有物件 ──────────────────────────────────────────────────
  
  list(
    dist_type    = dist_type,
    dist_label   = dist_label,
    dist_mat     = dist_mat,
    scores       = scores,
    var_exp      = var_exp,
    perm         = perm_res,
    change_df    = change_df,
    shift_df     = shift_df,
    shift_salt   = shift_salt,
    salt_change  = salt_change,
    fit_scores   = fit_scores,
    top_taxa     = top_taxa,
    genus_mat    = genus_mat_env,
    cor_results  = cor_results,
    cor_input    = cor_input,
    p_pcoa_trt   = p_pcoa_trt,
    p_pcoa_tp    = p_pcoa_tp,
    p_bc_change  = p_bc_change,
    p_arrow      = p_arrow,
    p_dpc1       = p_dpc1,
    p_salt       = p_salt,
    p_envfit     = p_envfit,
    p_cor_heatmap = p_cor_heatmap,
    p_rho_bar    = p_rho_bar,
    p_top_10     = combined_plot,
    p_cor_panel  = p_cor_panel
  )
}


# =============================================================================
#  執行兩種距離
# =============================================================================

cat("\n\n===== BRAY-CURTIS =====\n")
res_bray <- run_beta_diversity(dist_type = "bray")

cat("\n\n===== AITCHISON =====\n")
res_ait  <- run_beta_diversity(dist_type = "aitchison")


# =============================================================================
#  並排比較 panels
# =============================================================================

make_compare <- function(p_left, p_right, title) {
  (p_left | p_right) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )&
    theme(legend.position = "bottom")
}

p_compare_pcoa   <- make_compare(res_bray$p_pcoa_trt,  res_ait$p_pcoa_trt,
                                 "PCoA: Bray-Curtis vs Aitchison")
p_compare_arrow  <- make_compare(res_bray$p_arrow,     res_ait$p_arrow,
                                 "Directional shift: Bray-Curtis vs Aitchison")
p_compare_change <- make_compare(res_bray$p_bc_change, res_ait$p_bc_change,
                                 "Within-subject distance change: Bray-Curtis vs Aitchison")
p_compare_envfit <- make_compare(res_bray$p_envfit,    res_ait$p_envfit,
                                 "envfit vectors: Bray-Curtis vs Aitchison")
p_compare_salt   <- make_compare(res_bray$p_salt,      res_ait$p_salt,
                                 "ΔPC1 vs ΔSALT: Bray-Curtis vs Aitchison")
p_compare_panel   <- make_compare(res_bray$p_cor_panel,      res_ait$p_cor_panel,
                                 "Correlation between relative abundance and ΔSALT: Bray-Curtis vs Aitchison")

p_full_comparison <- (
  res_bray$p_pcoa_trt  | res_ait$p_pcoa_trt
) / (
  res_bray$p_arrow     | res_ait$p_arrow
) / (
  res_bray$p_bc_change | res_ait$p_bc_change
) +
  plot_annotation(
    title    = "Beta diversity full comparison",
    subtitle = "Left = Bray-Curtis  |  Right = Aitchison (CLR, genus-level)",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  )


# =============================================================================
#  Save
# =============================================================================.

library(Cairo)
save_fig <- function(plot, filename, w = 12, h = 6) {
  #dev.off()
  ggsave(filename, plot = plot, width = w, height = h, dpi = 300, bg = "white", device = cairo_pdf)
  cat(sprintf("  Saved: %s\n", filename))
}
dir.create(path = "outputs/tax_beta/", recursive = T)
# 比較圖
save_fig(p_compare_pcoa,    "outputs/tax_beta/beta_fig1_pcoa_comparison.pdf",    w = 14, h = 6)
save_fig(p_compare_arrow,   "outputs/tax_beta/beta_fig2_arrow_comparison.pdf",   w = 14, h = 7)
save_fig(p_compare_change,  "outputs/tax_beta/beta_fig3_change_comparison.pdf",  w = 14, h = 6)
save_fig(p_compare_envfit,  "outputs/tax_beta/beta_fig4_envfit_comparison.pdf",  w = 14, h = 7)
save_fig(p_compare_salt,    "outputs/tax_beta/beta_fig5_salt_comparison.pdf",    w = 12, h = 6)
save_fig(p_full_comparison, "outputs/tax_beta/beta_fig0_full_comparison.pdf",    w = 14, h = 18)

# 個別圖
for (res in list(res_bray, res_ait)) {
  tag <- res$dist_type
  save_fig(res$p_pcoa_trt,  sprintf("outputs/tax_beta/beta_%s_pcoa_treatment.pdf",  tag))
  save_fig(res$p_pcoa_tp,   sprintf("outputs/tax_beta/beta_%s_pcoa_timepoint.pdf",  tag))
  save_fig(res$p_bc_change, sprintf("outputs/tax_beta/beta_%s_distance_change.pdf", tag), w = 10, h = 6)
  save_fig(res$p_arrow,     sprintf("outputs/tax_beta/beta_%s_shift_arrow.pdf",     tag), w = 9,  h = 7)
  save_fig(res$p_dpc1,      sprintf("outputs/tax_beta/beta_%s_dPC1_bar.pdf",        tag), w = 8,  h = 5)
  save_fig(res$p_salt,      sprintf("outputs/tax_beta/beta_%s_dPC1_vs_dSALT.pdf",   tag), w = 7,  h = 6)
  save_fig(res$p_envfit,    sprintf("outputs/tax_beta/beta_%s_envfit.pdf",          tag), w = 9,  h = 7)
  save_fig(res$p_cor_panel, sprintf("outputs/tax_beta/beta_%s_cor_heatmap.pdf", tag), w = 14, h = 7)
  save_fig(res$p_top_10, sprintf("outputs/tax_beta/beta_%s_top10_genus.pdf", tag), w = 14, h = 7)
}


cat("\n✓ Beta diversity complete.\n")
cat("  Results: res_bray, res_ait\n")
cat("  Plots:   res_bray$p_pcoa_trt, res_ait$p_envfit, ...\n")
