


# =============================================================================
#  SECTION 2 — Beta Diversity PCoA
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

cat("\n── Section 2: Beta diversity ───────────────────────────────────────────\n")

# Bray-Curtis distance on rarefied table
bray_dist <- phyloseq::distance(ps_filt, method = "bray")

# PCoA ordination
ord_pcoa <- ordinate(ps_filt, method = "PCoA", distance = bray_dist)

# Variance explained
var_exp <- ord_pcoa$values$Relative_eig * 100
cat(sprintf("PC1: %.1f%%  PC2: %.1f%%  PC3: %.1f%%\n",
            var_exp[1], var_exp[2], var_exp[3]))

# Scores data frame
scores_df <- data.frame(
  ord_pcoa$vectors[, 1:3],
  sample_data(ps_filt)
) %>% rename(PC1 = Axis.1, PC2 = Axis.2, PC3 = Axis.3)

# ── 4a. PCoA coloured by treatment ───────────────────────────────────────────

p_pcoa_treatment <- ggplot(scores_df,
                            aes(x = PC1, y = PC2,
                                colour = Treatment, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_line(aes(group = PatientID), colour = "grey60",
            linewidth = 0.4, linetype = "dashed") +
  stat_ellipse(aes(group = Treatment), level = 0.75, linewidth = 0.8) +
  scale_colour_manual(values = col_treatment) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(
    title  = "PCoA (Bray-Curtis) — coloured by treatment",
    x      = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y      = sprintf("PC2 (%.1f%%)", var_exp[2]),
    colour = "Treatment",
    shape  = "Time point"
  ) +
  theme_fmt

# ── 4b. PCoA coloured by time point ──────────────────────────────────────────

p_pcoa_time <- ggplot(scores_df,
                       aes(x = PC1, y = PC2,
                           colour = Timepoint, shape = Treatment)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_line(aes(group = PatientID), colour = "grey60",
            linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(values = col_timepoint) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    title  = "PCoA (Bray-Curtis) — coloured by time point",
    x      = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y      = sprintf("PC2 (%.1f%%)", var_exp[2]),
    colour = "Time point",
    shape  = "Treatment"
  ) +
  theme_fmt


# ── 4c. PERMANOVA（含 blocks 控制 repeated measures）────────────────────────

meta_ordered <- meta_raw[sample_names(ps_filt), ]

# how() 限制置換只在 subject 間進行，正確處理 repeated measures
perm_ctrl <- how(
  blocks = meta_ordered$PatientID,
  nperm  = 999
)

set.seed(42)
perm_full <- adonis2(
  bray_dist ~ Treatment * Timepoint + Age,
  data         = meta_ordered,
  permutations = perm_ctrl
)
cat("\nPERMANOVA with blocks (full model):\n")
print(perm_full)

# Homogeneity of dispersions
bd <- betadisper(bray_dist, meta_ordered$Treatment)
cat("\nbetadisper (Treatment):\n")
print(permutest(bd, permutations = 999))

# Within-subject Bray-Curtis change
bc_mat <- as.matrix(bray_dist)
change_df <- meta_ordered %>%
  group_by(PatientID) %>%
  group_modify(~{
    sids <- .x$SampleID
    tp   <- .x$Timepoint
    expand.grid(from = sids, to = sids, stringsAsFactors = FALSE) %>%
      filter(from != to) %>%
      mutate(
        bc   = bc_mat[cbind(from, to)],
        from_tp = tp[match(from, sids)],
        to_tp   = tp[match(to,   sids)]
      )
  }) %>%
  filter(as.integer(from_tp) < as.integer(to_tp)) %>%
  left_join(meta_ordered %>% select(SampleID, Treatment), by = c("from" = "SampleID")) %>%
  ungroup()

p_bc_change <- ggplot(change_df,
                       aes(x = paste(from_tp, "→", to_tp),
                           y = bc, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  # geom_point + position_jitter(seed) 確保點與標籤位置完全一致
  geom_point(
    aes(colour = Treatment, shape = Treatment),
    position = position_jitter(width = 0.15, seed = 42),
    size = 3, alpha = 0.6
  ) +
  geom_text_repel(
    aes(label = PatientID, colour = Treatment),
    position       = position_jitter(width = 0.15, seed = 42),  # seed 必須相同
    segment.colour = NA,      # 隱藏連線
    box.padding    = 0.2,     # 文字貼近點
    point.padding  = 0.3,
    max.overlaps   = Inf,     # n=10 強制全部顯示
    size           = 3,
    show.legend    = FALSE
  ) +
  scale_fill_manual(values   = col_treatment) +
  scale_colour_manual(values = col_treatment) +
  scale_shape_manual(values  = c(FMT = 16, placebo = 17)) +
  labs(title = "Within-subject Bray-Curtis change",
       x = "Transition", y = "Bray-Curtis dissimilarity") +
  theme_fmt



# =============================================================================
#  SECTION 2b — Aitchison Distance PCoA
#  (Compositional-aware alternative to Bray-Curtis)
# =============================================================================

cat("\n── Section 2b: Aitchison PCoA ──────────────────────────────────────────\n")

# CLR 轉換（先加 pseudo-count 0.5 處理零值，再 agglomerate 到 genus 降低稀疏性）
ps_genus_clr <- tax_glom(ps_filt, taxrank = "Genus") %>%
  transform_sample_counts(function(x) x + 0.5) %>%
  transform_sample_counts(function(x) log(x) - mean(log(x)))   # CLR

# Aitchison distance = CLR 後的歐氏距離
clr_mat        <- t(as.matrix(otu_table(ps_genus_clr)))         # rows = samples
aitchison_dist <- dist(clr_mat, method = "euclidean")

# PCoA
ord_ait    <- ordinate(ps_genus_clr, method = "PCoA", distance = aitchison_dist)
var_exp_a  <- ord_ait$values$Relative_eig * 100
cat(sprintf("Aitchison PC1: %.1f%%  PC2: %.1f%%\n", var_exp_a[1], var_exp_a[2]))

scores_ait <- data.frame(
  ord_ait$vectors[, 1:2],
  sample_data(ps_genus_clr)
) %>% rename(PC1 = Axis.1, PC2 = Axis.2)

p_pcoa_ait <- ggplot(scores_ait,
                      aes(x = PC1, y = PC2,
                          colour = Treatment, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_line(aes(group = PatientID), colour = "grey60",
            linewidth = 0.4, linetype = "dashed") +
  stat_ellipse(aes(group = Treatment), level = 0.75, linewidth = 0.8) +
  scale_colour_manual(values = col_treatment) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(
    title    = "PCoA (Aitchison / CLR) — coloured by treatment",
    subtitle = "Genus-level CLR, pseudo-count = 0.5",
    x        = sprintf("PC1 (%.1f%%)", var_exp_a[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp_a[2]),
    colour   = "Treatment", shape = "Time point"
  ) +
  theme_fmt

# PERMANOVA on Aitchison distance
meta_ait <- meta_raw[rownames(clr_mat), ]
perm_ctrl_ait <- how(blocks = meta_ait$PatientID, nperm = 999)

set.seed(42)
perm_ait <- adonis2(
  aitchison_dist ~ Treatment * Timepoint + Age,
  data         = meta_ait,
  permutations = perm_ctrl_ait
)
cat("\nPERMANOVA (Aitchison):\n")
print(perm_ait)


# =============================================================================
#  SECTION 2c — PCoA Directional Shift (Baseline → 6M)
# =============================================================================

cat("\n── Section 2c: PCoA directional shift (Baseline→6M) ───────────────────\n")

# 使用 Bray-Curtis PCoA scores（scores_df 已在 Section 2 建立）
pcoa_scores <- scores_df   # rename for clarity in this section

shift_df <- pcoa_scores %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  select(PatientID, Treatment, Timepoint, PC1, PC2) %>%
  pivot_wider(names_from  = Timepoint,
              values_from = c(PC1, PC2)) %>%
  rename(
    PC1_bl = `PC1_Baseline`,
    PC1_6m = `PC1_6M post-FMT`,
    PC2_bl = `PC2_Baseline`,
    PC2_6m = `PC2_6M post-FMT`
  ) %>%
  mutate(
    dPC1 = PC1_6m - PC1_bl,
    dPC2 = PC2_6m - PC2_bl
  )

# 方向一致性摘要
cat("\nPC1 shift direction summary:\n")
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

# Wilcoxon：兩組 dPC1 是否有差異
wilcox_dpc1 <- wilcox.test(dPC1 ~ Treatment, data = shift_df, exact = FALSE)
cat(sprintf("\nWilcoxon (dPC1 FMT vs placebo): W = %.1f, p = %.4f\n",
            wilcox_dpc1$statistic, wilcox_dpc1$p.value))

# Wilcoxon Signed-Rank：Baseline→6M BC distance per group
bc_bl_6m <- change_df %>%
  filter(from_tp == "Baseline", to_tp == "6M post-FMT")

wilcox_bc_fmt <- wilcox.test(
  bc_bl_6m %>% filter(Treatment == "FMT")     %>% pull(bc),
  mu = 0, alternative = "greater", exact = FALSE
)
wilcox_bc_plc <- wilcox.test(
  bc_bl_6m %>% filter(Treatment == "placebo") %>% pull(bc),
  mu = 0, alternative = "greater", exact = FALSE
)
wilcox_bc_grp <- wilcox.test(bc ~ Treatment, data = bc_bl_6m, exact = FALSE)

cat(sprintf("Wilcoxon BC>0 (FMT):     p = %.4f\n", wilcox_bc_fmt$p.value))
cat(sprintf("Wilcoxon BC>0 (placebo): p = %.4f\n", wilcox_bc_plc$p.value))
cat(sprintf("Wilcoxon BC FMT vs placebo (Baseline→6M): p = %.4f\n", wilcox_bc_grp$p.value))

# SALT score 變化量
salt_change <- meta_raw %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  select(PatientID, Timepoint, SALT_score) %>%
  pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
  rename(SALT_bl = Baseline, SALT_6m = `6M post-FMT`) %>%
  mutate(dSALT = SALT_6m - SALT_bl)   # 負值 = 臨床改善

shift_salt <- left_join(shift_df, salt_change, by = "PatientID") %>%
  left_join(meta_raw %>% filter(Timepoint == "Baseline") %>%
              select(PatientID, Treatment) %>% distinct(),
            by = "PatientID", suffix = c("", ".meta")) %>%
  mutate(Treatment = coalesce(Treatment, Treatment.meta)) %>%
  select(-Treatment.meta)

# Spearman: dPC1 vs dSALT
cor_shift <- cor.test(shift_salt$dPC1, shift_salt$dSALT, method = "spearman", exact = FALSE)
cat(sprintf("\nSpearman dPC1 ~ dSALT: rho = %.3f, p = %.4f\n",
            cor_shift$estimate, cor_shift$p.value))

# ── Plots ─────────────────────────────────────────────────────────────────────

# Arrow plot：Baseline → 6M 位移向量
p_shift_arrow <- ggplot() +
  geom_point(
    data    = pcoa_scores %>% filter(Timepoint == "Baseline"),
    aes(x = PC1, y = PC2, colour = Treatment),
    size = 3.5, shape = 1, stroke = 1.2
  ) +
  geom_segment(
    data = shift_df,
    aes(x = PC1_bl, y = PC2_bl,
        xend = PC1_6m, yend = PC2_6m,
        colour = Treatment),
    arrow     = arrow(length = unit(0.015, "npc"), type = "closed"),
    linewidth = 0.9, alpha = 0.85
  ) +
  geom_text_repel(
    data = shift_df,
    aes(x = PC1_6m, y = PC2_6m, label = PatientID, colour = Treatment),
    size = 3, segment.colour = NA
  ) +
  scale_colour_manual(values = col_treatment) +
  labs(
    title    = "PCoA displacement: Baseline → 6M post-FMT",
    subtitle = "Circle = baseline, arrowhead = 6M",
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2])
  ) +
  theme_fmt

# dPC1 bar per patient
p_dPC1 <- ggplot(shift_df,
                  aes(x = reorder(PatientID, dPC1), y = dPC1, fill = Treatment)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = col_treatment) +
  labs(title = "PC1 shift (6M − Baseline) per patient",
       x = NULL, y = "ΔPC1") +
  theme_fmt +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# dPC1 vs dSALT scatter
p_shift_salt <- ggplot(shift_salt,
                        aes(x = dPC1, y = dSALT, colour = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_text_repel(aes(label = PatientID), size = 3, segment.colour = NA) +
  scale_colour_manual(values = col_treatment) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("ρ = %.3f\np = %.3f", cor_shift$estimate, cor_shift$p.value),
           hjust = 1.1, vjust = 1.5, size = 3.5, colour = "grey30") +
  labs(
    title    = "PC1 shift vs SALT score change",
    subtitle = "Negative dSALT = clinical improvement",
    x        = "ΔPC1 (6M − Baseline)",
    y        = "ΔSALT (6M − Baseline)"
  ) +
  theme_fmt


# =============================================================================
#  SECTION 2d — envfit（找出驅動 PCoA 位移的 genus）
# =============================================================================

cat("\n── Section 2d: envfit ──────────────────────────────────────────────────\n")

# Genus relative abundance matrix（僅 Baseline & 6M）
ps_genus_rel_env <- tax_glom(ps_filt, taxrank = "Genus") %>%
  transform_sample_counts(function(x) x / sum(x))

genus_mat <- as.data.frame(t(otu_table(ps_genus_rel_env)))
tax_df    <- as.data.frame(tax_table(ps_genus_rel_env))
colnames(genus_mat) <- tax_df$Genus

# 只保留 Baseline 和 6M 樣本
samples_bl_6m <- pcoa_scores %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  pull(SampleID)

pcoa_scores_mat <- pcoa_scores %>%
  filter(SampleID %in% samples_bl_6m) %>%
  column_to_rownames("SampleID") %>%
  select(PC1, PC2) %>%
  as.matrix()

genus_mat_sub <- genus_mat[rownames(pcoa_scores_mat), ]

set.seed(42)
fit       <- envfit(pcoa_scores_mat, genus_mat_sub, permutations = 999, na.rm = TRUE)

fit_scores <- as.data.frame(scores(fit, display = "vectors")) %>%
  tibble::rownames_to_column("Genus") %>%
  mutate(
    r2   = fit$vectors$r,
    pval = fit$vectors$pvals
  ) %>%
  filter(pval < 0.05) %>%
  arrange(desc(r2))

cat(sprintf("Significant genera (envfit p<0.05): %d\n", nrow(fit_scores)))
print(fit_scores)

top_taxa    <- fit_scores %>% slice_max(r2, n = 10)
scale_factor <- 0.3

p_envfit <- p_shift_arrow +
  geom_segment(
    data = top_taxa,
    aes(x = 0, y = 0,
        xend = PC1 * scale_factor,
        yend = PC2 * scale_factor),
    arrow       = arrow(length = unit(0.012, "npc"), type = "closed"),
    colour      = "grey20",
    linewidth   = 0.7,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = top_taxa,
    aes(x = PC1 * scale_factor,
        y = PC2 * scale_factor,
        label = Genus),
    colour         = "grey10",
    size           = 3,
    fontface       = "italic",
    segment.colour = NA,
    inherit.aes    = FALSE
  ) +
  labs(subtitle = "Arrows: envfit vectors (p<0.05, top 10 r²) | Circle=Baseline, head=6M")


# =============================================================================
#  SECTION 2e — Genus × SALT Correlation Heatmap（by treatment group）
# =============================================================================

cat("\n── Section 2e: Genus–SALT correlation heatmap ──────────────────────────\n")

top_genera_names <- top_taxa$Genus

# Genus 豐度變化量（6M − Baseline）
genus_change <- genus_mat %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(
    meta_raw %>% select(SampleID, PatientID, Timepoint),
    by = "SampleID"
  ) %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  select(PatientID, Timepoint, all_of(top_genera_names)) %>%
  pivot_longer(cols = all_of(top_genera_names),
               names_to = "Genus", values_to = "abund") %>%
  pivot_wider(names_from = Timepoint, values_from = abund) %>%
  rename(bl = Baseline, m3 = `6M post-FMT`) %>%
  mutate(delta_abund = m3 - bl) %>%
  select(PatientID, Genus, delta_abund) %>%
  pivot_wider(names_from = Genus, values_from = delta_abund)

cor_input <- genus_change %>%
  left_join(salt_change %>% select(PatientID, dSALT), by = "PatientID") %>%
  column_to_rownames("PatientID")

# Spearman per treatment group（n=5 per group，exploratory）
cor_results <- c("FMT", "placebo") %>%
  purrr::map_dfr(function(trt) {
    patients_in_group <- meta_raw %>%
      filter(Timepoint == "Baseline", Treatment == trt) %>%
      pull(PatientID)

    sub_input <- cor_input %>%
      tibble::rownames_to_column("PatientID") %>%
      filter(PatientID %in% patients_in_group) %>%
      column_to_rownames("PatientID")

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

cat("\nSpearman correlation by treatment group (n=5, exploratory):\n")
print(cor_results)

# Heatmap matrix（z-score）
heatmap_mat <- cor_input %>%
  select(all_of(top_genera_names)) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(dSALT = scale(cor_input$dSALT)[, 1]) %>%
  tibble::rownames_to_column("PatientID") %>%
  left_join(
    meta_raw %>% filter(Timepoint == "Baseline") %>%
      select(PatientID, Treatment) %>% distinct(),
    by = "PatientID"
  )

# Genus 欄位按照 FMT 組的 rho 排序
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
  geom_text(aes(label = sig), size = 4, colour = "white", fontface = "bold") +
  geom_vline(xintercept = 1.5, colour = "grey30", linewidth = 1) +
  facet_grid(Treatment ~ var_type, scales = "free", space = "free") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#D6604D",
    midpoint = 0,
    name     = "z-score"
  ) +
  labs(
    title    = "Δ Genus abundance vs Δ SALT score (6M − Baseline)",
    subtitle = "* p<0.05  † p<0.1  |  Columns sorted by Spearman ρ (FMT group)",
    x        = NULL, y = NULL
  ) +
  theme_fmt +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, face = "italic"),
    strip.text.x = element_text(face = "bold"),
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
  labs(title = "Spearman ρ\n(Δabundance vs ΔSALT)",
       x = NULL, y = "ρ") +
  theme_fmt +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic", size = 9))

p_cor_panel <- p_cor_heatmap + p_rho_bar +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "Genus–SALT correlation panel (by treatment group)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )



cat("\n── Section 3: Taxa composition ─────────────────────────────────────────\n")

# Aggregate to genus level
ps_genus <- tax_glom(ps_filt, taxrank = "Genus", NArm = FALSE)

# Relative abundance
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x) * 100)

# Top 10 genera across all samples
top10_genera <- names(sort(taxa_sums(ps_genus_rel), decreasing = TRUE))[1:10]

# Melt to long format
genus_melt <- psmelt(ps_genus_rel) %>%
  mutate(Genus2 = if_else(OTU %in% top10_genera, as.character(Genus), "Other")) %>%
  group_by(SampleID, Genus2, Timepoint, Treatment, PatientID) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

genus_summary <- genus_melt %>%
  group_by(Treatment, Timepoint, Genus2) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

genus_cols <- c(
  setNames(
    colorRampPalette(brewer.pal(10, "Paired"))(10),
    unique(genus_melt$Genus2[genus_melt$Genus2 != "Other"])
  ),
  Other = "grey80"
)

p_barplot <- ggplot(genus_summary,
                     aes(x = Timepoint, y = mean_abund, fill = Genus2)) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~ Treatment) +
  scale_fill_manual(values = genus_cols) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "Genus-level relative abundance",
       x = NULL, y = "Mean relative abundance (%)", fill = "Genus") +
  theme_fmt +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Firmicutes / Bacteroidetes ratio
ps_phylum <- tax_glom(ps_filt, taxrank = "Phylum", NArm = FALSE)
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x) * 100)

fb_df <- psmelt(ps_phylum_rel) %>%
  filter(Phylum %in% c("Firmicutes", "Bacteroidetes")) %>%
  select(SampleID, Phylum, Abundance, Timepoint, Treatment, PatientID) %>%
  pivot_wider(names_from = Phylum, values_from = Abundance) %>%
  mutate(FB_ratio = Firmicutes / (Bacteroidetes + 0.01))

p_fb <- ggplot(fb_df, aes(x = Timepoint, y = FB_ratio,
                            colour = Treatment, group = PatientID)) +
  geom_line(alpha = 0.4, linewidth = 0.6) +
  geom_point(size = 2.5) +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "line", linewidth = 1.5) +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "point", size = 4, shape = 18) +
  scale_colour_manual(values = col_treatment) +
  labs(title = "Firmicutes / Bacteroidetes ratio",
       x = NULL, y = "F/B ratio") +
  theme_fmt


# =============================================================================
#  SECTION 4 — SALT Score Correlation
# =============================================================================

cat("\n── Section 4: SALT score correlation ───────────────────────────────────\n")

# Join alpha diversity with SALT scores
salt_alpha <- alpha_df %>% select(SampleID, Shannon, Chao1, Timepoint, Treatment,
                                   PatientID, SALT_score)

# Scatter: Shannon vs SALT coloured by time point
p_salt_shannon <- ggplot(salt_alpha,
                          aes(x = Shannon, y = SALT_score,
                              colour = Timepoint, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              colour = "grey40", fill = "grey85") +
  geom_text_repel(aes(label = PatientID), size = 2.5, max.overlaps = 6) +
  scale_colour_manual(values = col_timepoint) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Shannon diversity vs SALT score",
       x = "Shannon index", y = "SALT score (lower = better)",
       colour = "Time point", shape = "Treatment") +
  theme_fmt

# Spearman correlation per time point
cat("\nSpearman correlation (Shannon ~ SALT) per time point:\n")
salt_alpha %>%
  group_by(Timepoint) %>%
  summarise(
    rho   = cor(Shannon, SALT_score, method = "spearman"),
    pval  = cor.test(Shannon, SALT_score, method = "spearman")$p.value,
    n     = n(),
    .groups = "drop"
  ) %>%
  print()

# Mixed model: SALT ~ Shannon + Timepoint + Treatment + (1|PatientID)
lmm_salt <- lmer(
  SALT_score ~ Shannon + Timepoint + Treatment + Age + (1 | PatientID),
  data = salt_alpha, REML = TRUE
)
cat("\nLMM SALT score — fixed effects:\n")
print(summary(lmm_salt)$coefficients)

# Trajectories per patient coloured by treatment
p_salt_traj <- ggplot(salt_alpha,
                       aes(x = Timepoint, y = SALT_score,
                           colour = Treatment, group = PatientID)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 2.5) +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "line", linewidth = 1.5) +
  stat_summary(aes(group = Treatment), fun.data = mean_se,
               geom = "ribbon", alpha = 0.15, colour = NA,
               aes(fill = Treatment)) +
  scale_colour_manual(values = col_treatment) +
  scale_fill_manual(values   = col_treatment) +
  labs(title = "SALT score trajectory",
       x = NULL, y = "SALT score", colour = "Treatment") +
  theme_fmt


# =============================================================================
#  SECTION 5 — Differential Taxa (DESeq2)
# =============================================================================

cat("\n── Section 5: Differential taxa (DESeq2) ───────────────────────────────\n")

run_deseq2 <- function(ps_obj, contrast_var, ref_level,
                        subset_expr = NULL, tag = "") {
  if (!is.null(subset_expr)) {
    ps_obj <- subset_samples(ps_obj, eval(parse(text = subset_expr)))
  }

  sample_data(ps_obj)[[contrast_var]] <- relevel(
    factor(sample_data(ps_obj)[[contrast_var]]), ref = ref_level
  )

  # DESeq2 needs un-rarefied counts
  dds <- phyloseq_to_deseq2(ps_obj, as.formula(paste("~", contrast_var)))
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)

  res <- results(dds, contrast = c(contrast_var, setdiff(levels(dds[[contrast_var]]), ref_level), ref_level),
                 alpha = 0.05) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ASV_ID") %>%
    left_join(
      tax_table(ps_obj) %>% as.data.frame() %>% tibble::rownames_to_column("ASV_ID"),
      by = "ASV_ID"
    ) %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  cat(sprintf("\n[DESeq2 %s] Significant ASVs (padj<0.05): %d\n", tag,
              sum(res$padj < 0.05, na.rm = TRUE)))
  return(res)
}

# 5a. FMT vs Placebo at month 3 (primary comparison)
ps_m3   <- subset_samples(ps_filt, Timepoint == "6M post-FMT")
res_m3  <- run_deseq2(ps_m3, "Treatment", "placebo", tag = "FMT vs Placebo @ month6")

# 5b. Within FMT: baseline vs month 3
ps_fmt  <- subset_samples(ps_filt, Treatment == "FMT")
res_fmt <- run_deseq2(ps_fmt, "Timepoint", "Baseline", tag = "FMT baseline→month6")

# ── Volcano plot ──────────────────────────────────────────────────────────────

make_volcano <- function(res_df, title) {
  res_df <- res_df %>%
    mutate(
      sig   = case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "Up in FMT",
        padj < 0.05 & log2FoldChange < -1 ~ "Down in FMT",
        TRUE ~ "NS"
      ),
      label = if_else(padj < 0.01, as.character(Genus), NA_character_)
    )

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                      colour = sig, label = label)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(size = 2.8, max.overlaps = 12, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c(
      "Up in FMT"   = col_treatment["FMT"],
      "Down in FMT" = col_treatment["placebo"],
      "NS"          = "grey70"
    )) +
    labs(title = title,
         x = "Log₂ fold change", y = "-log₁₀(adj. p-value)",
         colour = NULL) +
    theme_fmt
}

p_volcano_m3  <- make_volcano(res_m3,  "FMT vs placebo @ month 3")
p_volcano_fmt <- make_volcano(res_fmt, "FMT group: baseline → month 3")

# ── LFC bar chart (top 15 significant genera) ─────────────────────────────────

make_lfc_bar <- function(res_df, n = 15, title) {
  top <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(n) %>%
    mutate(
      Genus2    = coalesce(Genus, paste0("ASV_", str_sub(ASV_ID, 1, 6))),
      direction = if_else(log2FoldChange > 0, "Up in FMT", "Down in FMT")
    )

  ggplot(top, aes(x = reorder(Genus2, log2FoldChange),
                   y = log2FoldChange, fill = direction)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                      ymax = log2FoldChange + lfcSE),
                  width = 0.25, colour = "grey40") +
    coord_flip() +
    scale_fill_manual(values = c(
      "Up in FMT"   = col_treatment["FMT"],
      "Down in FMT" = col_treatment["placebo"]
    )) +
    labs(title = title, x = NULL, y = "Log₂ fold change (±SE)",
         fill = NULL) +
    theme_fmt
}

p_lfc_m3  <- make_lfc_bar(res_m3,  title = "Top diff. genera: FMT vs placebo @ month 3")
p_lfc_fmt <- make_lfc_bar(res_fmt, title = "Top diff. genera: FMT baseline → month 3")


# =============================================================================
#  SECTION 6 — Demographic summary plots
# =============================================================================

cat("\n── Section 6: Demographics ─────────────────────────────────────────────\n")

# 每位病人只取一列（baseline），避免重複計算
demo <- meta_raw %>%
  filter(Timepoint == "Baseline") %>%          # factor label after formatting
  select(PatientID, Treatment, Age, Race, Sex) %>%
  distinct()

# ── 6a. Treatment group counts ────────────────────────────────────────────────
p_demo_treatment <- ggplot(demo, aes(x = Treatment, fill = Treatment)) +
  geom_bar(width = 0.55, colour = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = col_treatment) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Treatment allocation", x = NULL, y = "n patients") +
  theme_fmt + theme(legend.position = "none")

# ── 6b. Race / ethnicity composition ─────────────────────────────────────────
race_counts <- demo %>%
  count(Treatment, Race) %>%
  group_by(Treatment) %>%
  mutate(pct = n / sum(n) * 100)

p_demo_race <- ggplot(race_counts,
                       aes(x = Treatment, y = pct, fill = Race)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Race / ethnicity", x = NULL, y = "% of group", fill = "Race") +
  theme_fmt

# ── 6c. Sex composition ───────────────────────────────────────────────────────
sex_counts <- demo %>%
  count(Treatment, Sex) %>%
  group_by(Treatment) %>%
  mutate(pct = n / sum(n) * 100)

p_demo_sex <- ggplot(sex_counts,
                      aes(x = Treatment, y = pct, fill = Sex)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c(Male = "#4A90D9", Female = "#E87D72")) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Sex", x = NULL, y = "% of group", fill = "Sex") +
  theme_fmt

# ── 6d. Age distribution ──────────────────────────────────────────────────────
p_demo_age <- ggplot(demo, aes(x = Treatment, y = Age, fill = Treatment)) +
  geom_boxplot(width = 0.45, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(colour = Treatment), width = 0.1, size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = PatientID), size = 2.8,
                  box.padding = 0.4, segment.colour = "grey60") +
  scale_fill_manual(values   = col_treatment) +
  scale_colour_manual(values = col_treatment) +
  labs(title = "Age distribution", x = NULL, y = "Age (years)") +
  theme_fmt + theme(legend.position = "none")

# ── 6e. Combined demographic panel ───────────────────────────────────────────
p_demo_panel <- (p_demo_treatment | p_demo_sex | p_demo_race | p_demo_age) +
  plot_annotation(
    title = "Patient demographics",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

cat("Demographic summary (baseline only):\n")
demo %>%
  group_by(Treatment) %>%
  summarise(
    n          = n(),
    Age_mean   = round(mean(Age, na.rm = TRUE), 1),
    Age_sd     = round(sd(Age,  na.rm = TRUE), 1),
    .groups    = "drop"
  ) %>% print()


# =============================================================================
#  SECTION 7 — Save all plots
# =============================================================================

cat("\n── Saving figures ──────────────────────────────────────────────────────\n")

save_fig <- function(plot, filename, w = 10, h = 6) {
  ggsave(filename, plot = plot, width = w, height = h,
         dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", filename))
}

# Fig 0 — Demographics
save_fig(p_demo_panel,     "fig0_demographics.pdf",          w = 14, h = 5)
save_fig(p_demo_age,       "fig0b_age_distribution.pdf",     w = 6,  h = 5)

# Fig 1 — Alpha diversity
save_fig(p_alpha_combined, "fig1_alpha_diversity.pdf", w = 14, h = 8)

# Fig 2 — Beta diversity
save_fig(p_pcoa_treatment, "fig2a_pcoa_bray_treatment.pdf")
save_fig(p_pcoa_time,      "fig2b_pcoa_bray_timepoint.pdf")
save_fig(p_bc_change,      "fig2c_bc_within_change.pdf",        w = 10, h = 6)
save_fig(p_pcoa_ait,       "fig2d_pcoa_aitchison.pdf")
save_fig(p_shift_arrow,    "fig2e_pcoa_shift_arrow.pdf",        w = 9,  h = 7)
save_fig(p_dPC1,           "fig2f_dPC1_per_patient.pdf",        w = 8,  h = 5)
save_fig(p_shift_salt,     "fig2g_dPC1_vs_dSALT.pdf",          w = 7,  h = 6)
save_fig(p_envfit,         "fig2h_envfit_overlay.pdf",          w = 9,  h = 7)
save_fig(p_cor_panel,      "fig6_genus_salt_correlation.pdf",   w = 14, h = 7)

# Fig 3 — Taxa composition
save_fig(p_barplot,        "fig3_genus_barplot.pdf",          w = 12, h = 6)
save_fig(p_fb,             "fig3b_firmicutes_bacteroidetes.pdf", w = 10, h = 5)

# Fig 4 — SALT score
save_fig(p_salt_shannon,   "fig4a_salt_shannon_scatter.pdf")
save_fig(p_salt_traj,      "fig4b_salt_trajectory.pdf")

# Fig 5 — Differential taxa
save_fig(p_volcano_m3,     "fig5a_volcano_fmt_vs_placebo_m3.pdf")
save_fig(p_volcano_fmt,    "fig5b_volcano_fmt_time.pdf")
save_fig(p_lfc_m3,         "fig5c_lfc_bar_m3.pdf",           w = 10, h = 6)
save_fig(p_lfc_fmt,        "fig5d_lfc_bar_fmt.pdf",           w = 10, h = 6)

# Combined poster panel
poster <- (p_demo_panel) /
  (p_alpha_combined) /
  (p_pcoa_treatment | p_salt_traj) /
  (p_lfc_m3 | p_barplot) +
  plot_annotation(
    title    = "Alopecia oral FMT — Gut microbiome dynamics",
    subtitle = "n=10 | 3 time points | FMT vs placebo",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )
save_fig(poster, "fig_poster_panel.pdf", w = 16, h = 22)

# Save DESeq2 results tables
write.csv(res_m3,  "deseq2_fmt_vs_placebo_month6.csv", row.names = FALSE)
write.csv(res_fmt, "deseq2_fmt_baseline_vs_month6.csv", row.names = FALSE)
write.csv(alpha_df, "alpha_diversity_table.csv", row.names = FALSE)

cat("\nAll done! Output files written to working directory.\n")


# =============================================================================
#  QUICK-START: simulate data to verify pipeline runs
# =============================================================================
# Uncomment the block below to generate synthetic data and test the script
# without your real files.
#
# source("simulate_data.R")   # creates asv_table.csv, taxonomy.csv, metadata.csv
