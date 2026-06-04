# =============================================================================
#  16S ASV Analysis — Alopecia Oral FMT Study
#  n=10 patients | 3 time points | FMT vs Placebo
#
#  Priority order:
#    1. Alpha diversity (Shannon, Chao1)
#    2. Beta diversity PCoA (Bray-Curtis)
#    3. FMT vs Placebo group comparison
#    4. SALT score ~ microbiome correlation
#    5. Differential taxa (DESeq2)
#
#  Input files expected:
#    asv_table.csv    — rows = ASVs, cols = sample IDs (first col = ASV_ID)
#    taxonomy.csv     — cols: ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species
#    metadata.csv     — cols: SampleID, PatientID, Timepoint, Treatment, SALT_score, Age, Race
#                       Timepoint values: "baseline", "month2", "month3"
#                       Treatment values: "FMT", "placebo"
# =============================================================================


# ── 0. Install & load packages ───────────────────────────────────────────────

packages <- c(
  "phyloseq", "vegan", "DESeq2", "ggplot2", "dplyr", "tidyr",
  "patchwork", "ggrepel", "RColorBrewer", "scales", "lme4",
  "lmerTest", "emmeans", "broom", "stringr"
)

# Bioconductor packages
bioc_packages <- c("phyloseq", "DESeq2")

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(missing, ask = FALSE)
    } else {
      install.packages(missing)
    }
  }
}

cran_pkgs <- setdiff(packages, bioc_packages)
install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_packages, bioc = TRUE)

suppressPackageStartupMessages({
  lapply(packages, library, character.only = TRUE)
})


# ── 1. Load data ─────────────────────────────────────────────────────────────

# Set working directory to project root
setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

asv_raw  <- read.csv("inputs/asv_table.csv",  row.names = 1, check.names = FALSE)
tax_raw  <- read.csv("inputs/taxonomy.csv",   stringsAsFactors = FALSE)
meta_raw <- read.csv("inputs/metadata.csv",   stringsAsFactors = FALSE)

# ── 1a. Validate and format input ───────────────────────────────────────────

# ASV table: ensure integer counts
asv_mat <- as.matrix(asv_raw)
mode(asv_mat) <- "integer"

# Taxonomy: rownames must match ASV table rownames
tax_mat <- as.matrix(tax_raw[, -1])
rownames(tax_mat) <- tax_raw$ASV_ID
colnames(tax_mat) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_mat <- tax_mat[rownames(asv_mat), ]          # align order

# Metadata: rownames = SampleID
fast_responder <- c("GM01","GM05","GM06")
meta_raw <- meta_raw %>% mutate(fast_responder = if_else(PatientID %in% fast_responder, TRUE, FALSE))
rownames(meta_raw) <- meta_raw$SampleID
meta_raw <- meta_raw[colnames(asv_mat), ]        # align order

# Ordered factor for time points
meta_raw$Timepoint <- factor(
  meta_raw$Timepoint,
  levels = c("baseline", "month2", "month6"),
  labels = c("Baseline", "2M post-FMT", "6M post-FMT")
)
meta_raw$Treatment <- factor(meta_raw$Treatment, levels = c("FMT", "placebo"))


# ── 2. Build phyloseq object ─────────────────────────────────────────────────

ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_raw)
)

cat("Phyloseq object created:\n")
print(ps)
cat(sprintf("  Total reads: %s\n", formatC(sum(otu_table(ps)), format = "d", big.mark = ",")))
cat(sprintf("  Mean reads/sample: %s\n", formatC(mean(sample_sums(ps)), format = "d", big.mark = ",")))


# ── 2a. Filtering ────────────────────────────────────────────────────────────
rowSums(asv_mat) %>% quantile(probs = c(0.25, 0.5, 0.75)) # c(0,4,81) use 4 
# Remove ASVs present in < 2 samples or with < 4 total counts
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= 2 & sum(x) >= 4, TRUE)
cat(sprintf("After filtering: %d ASVs (from %d)\n", ntaxa(ps_filt), ntaxa(ps)))

# Rarefy to even depth for alpha/beta diversity
# set.seed(42)
# rarefaction_depth <- floor(min(sample_sums(ps_filt)) * 0.9)
# cat(sprintf("Rarefying to: %d reads\n", rarefaction_depth))
# ps_rare <- rarefy_even_depth(ps_filt, sample.size = rarefaction_depth, verbose = FALSE)


# ── 3. Colour palette ────────────────────────────────────────────────────────

col_treatment <- c(FMT = "#1D9E75", placebo = "#D85A30")
col_timepoint <- c(Baseline = "#378ADD", `2M post-FMT` = "#9F3FAF", `6M post-FMT` = "#BA7517")
tp_labels     <- c("Baseline", "2M post-FMT", "6M post-FMT")

theme_fmt <- theme_bw(base_size = 11) + 
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black") # 強制指定 axis.line 為 element_line 預防報錯
  )

# =============================================================================
#  SECTION 1 — Alpha Diversity
# =============================================================================

cat("\n── Section 1: Alpha diversity ──────────────────────────────────────────\n")

alpha_df <- estimate_richness(ps_filt, measures = c("Shannon", "Chao1", "Observed")) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(meta_raw %>% select(SampleID, PatientID, Timepoint, Treatment, SALT_score, Age, Race, fast_responder),
            by = "SampleID")

# Summary statistics
alpha_summary <- alpha_df %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    Shannon_mean = mean(Shannon),
    Shannon_sd   = sd(Shannon),
    Chao1_mean   = mean(Chao1),
    Chao1_sd     = sd(Chao1),
    .groups = "drop"
  )
cat("Alpha diversity summary:\n")
print(alpha_summary)

# ── 3a. Trajectory plots ─────────────────────────────────────────────────────

p_shannon <- ggplot(alpha_df, aes(x = Timepoint, y = Shannon,
                                  colour = Treatment, group = PatientID)) +
  geom_line(alpha = 0.35, linewidth = 0.6) +
  geom_point(alpha = 0.6, size = 2) +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "line", linewidth = 1.4, linetype = "solid") +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "point", size = 4, shape = 18) +
  scale_colour_manual(values = col_treatment) +
  labs(title = "Shannon diversity across time points",
       x = NULL, y = "Shannon index", colour = "Treatment") +
  theme_fmt + facet_wrap(~Treatment) + 
  geom_text(data = filter(alpha_df, Timepoint == "6M post-FMT"),aes(label =PatientID ), position = position_dodge(width =0.1))

p_chao1 <- ggplot(alpha_df, aes(x = Timepoint, y = Chao1,
                                  colour = Treatment, group = PatientID)) +
  geom_line(alpha = 0.35, linewidth = 0.6) +
  geom_point(alpha = 0.6, size = 2) +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "line", linewidth = 1.4, linetype = "solid") +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "point", size = 4, shape = 18) +
  scale_colour_manual(values = col_treatment) +
  labs(title = "Chao1 richness across time points",
       x = NULL, y = "Chao1", colour = "Treatment") +
  theme_fmt+ facet_wrap(~Treatment) + 
  geom_text(data = filter(alpha_df, Timepoint == "6M post-FMT"),aes(label =PatientID ), position = position_dodge(width =0.1))


p_alpha_combined <- p_shannon + p_chao1 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ── 3b. Linear mixed model: Shannon ~ Timepoint * Treatment + (1|PatientID) ──
# Linear Mixed-Effects Regression
lmm_shannon <- lmer(Shannon ~ Timepoint * Treatment + Age + (1 | PatientID),
                    data = alpha_df, REML = TRUE)
cat("\nLMM Shannon — fixed effects:\n")
print(summary(lmm_shannon)$coefficients)

# Pairwise comparisons within each treatment
emm_shannon <- emmeans(lmm_shannon, ~ Timepoint | Treatment)
cat("\nEMM pairwise (Shannon within treatment):\n")
print(pairs(emm_shannon))


# =============================================================================
#  SECTION 2 — Beta Diversity PCoA
# =============================================================================

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
) %>% dplyr::rename(PC1 = Axis.1, PC2 = Axis.2, PC3 = Axis.3)

# ── 4a. PCoA coloured by treatment ───────────────────────────────────────────

p_pcoa_treatment <- ggplot(scores_df,
                            aes(x = PC1, y = PC2,
                                colour = Treatment)) +
  geom_point(aes(shape = factor(Timepoint)), size = 3.5, alpha = 0.85) +
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
  theme_fmt + theme(aspect.ratio = 1)+
  geom_text(data = filter(scores_df, Timepoint == "Baseline"),aes(label =PatientID ), nudge_x = 0.01, nudge_y = -0.01)

# ── 4b. PCoA coloured by time point ──────────────────────────────────────────

p_pcoa_time <- ggplot(scores_df,
                       aes(x = PC1, y = PC2,
                           colour = Timepoint)) +
  geom_point(aes(shape = Treatment), size = 3.5, alpha = 0.85) +
  geom_line(aes(group = PatientID), colour = "grey60",
            linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(values = col_timepoint) +
  stat_ellipse(aes(group = Timepoint), level = 0.75, linewidth = 0.8) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    title  = "PCoA (Bray-Curtis) — coloured by time point",
    x      = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y      = sprintf("PC2 (%.1f%%)", var_exp[2]),
    colour = "Time point",
    shape  = "Treatment"
  ) +
  theme_fmt+ theme(aspect.ratio = 1)+
  geom_text(data = filter(scores_df, Timepoint == "Baseline"),aes(label =PatientID ), nudge_x = 0.01, nudge_y = -0.01)


# ── 4c. PERMANOVA ─────────────────────────────────────────────────────────────

meta_ordered <- meta_raw[sample_names(ps_filt), ]

# Full model

perm_full <- adonis2(
  bray_dist ~ Treatment * Timepoint + Age,
  data  = meta_ordered,
  permutations = 999,  by  = "margin")

# 考慮重複受試者
set.seed(42)
perm_ctrl <- how(
  blocks = meta_ordered$PatientID,   # 置換只在 block 內進行
  nperm  = 999
)

perm_full <- adonis2(
  bray_dist ~ Treatment * Timepoint + Age,
  data         = meta_ordered,
  permutations = perm_ctrl,
  by = "margin"
)
cat("\nPERMANOVA (full model):\n")
print(perm_full)

# Homogeneity of dispersions : FMT vs placebo
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

p_bc_change <- ggplot(change_df, aes(x = paste(from_tp, "→", to_tp), y = bc, fill = Treatment)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + # 隱藏 boxplot 的離群點避免重複
  geom_point(
    aes(colour = Treatment, shape = Treatment),
    # 使用 position_jitterdodge 讓點對齊分組的 box
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75, seed = 42), 
    size = 3, 
    alpha = 0.6
  ) + 
  geom_text_repel(
    aes(label = PatientID, colour = Treatment),
    # 這裡的 position 參數必須與 geom_point 完全一致
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75, seed = 42), 
    segment.colour = NA, 
    box.padding = 0.2, 
    point.padding = 0.3, 
    max.overlaps = Inf, 
    size = 3, 
    show.legend = FALSE,
    fontface = "bold"
  ) + 
  scale_fill_manual(values = col_treatment) + 
  scale_colour_manual(values = col_treatment) + 
  scale_shape_manual(values = c(FMT = 16, placebo = 17)) + 
  labs(title = "Within-subject Bray-Curtis change", x = "Transition", y = "Bray-Curtis dissimilarity") + 
  theme_fmt


# 用 wilcoxon signed-rank 看是否在 樣本數小/paired/b-c dist有0-1的界線不適合做t-test

## 第一層 — FMT vs placebo 的 BC change（baseline→6M）組間比較
bc_bl_2m <- change_df %>%
  filter(from_tp == "Baseline", to_tp == "2M post-FMT")

wilcox.test(bc ~ Treatment, data = bc_bl_2m,
            exact = FALSE) 

bc_bl_6m <- change_df %>%
  filter(from_tp == "Baseline", to_tp == "6M post-FMT")

wilcox.test(bc ~ Treatment, data = bc_bl_6m,
            exact = FALSE)   # unpaired，比較兩組的 shift 幅度

## 第二層 — 每組內部 baseline vs 6M 的 within-subject change
# FMT 組：每個病人的 baseline→6M BC 值 or baseline→2M

fmt_bc <- bc_bl_2m %>% filter(Treatment == "FMT") %>% pull(bc)
plc_bc <- bc_bl_2m %>% filter(Treatment == "placebo") %>% pull(bc)

wilcox.test(fmt_bc, mu = 0, alternative = "greater")   # FMT 組 BC change 顯著 > 0？ 0.0625
wilcox.test(plc_bc, mu = 0, alternative = "greater")   # placebo 組呢？ 0.03125


fmt_bc <- bc_bl_6m %>% filter(Treatment == "FMT") %>% pull(bc)
plc_bc <- bc_bl_6m %>% filter(Treatment == "placebo") %>% pull(bc)

wilcox.test(fmt_bc, mu = 0, alternative = "greater")   # FMT 組 BC change 顯著 > 0？0.0625
wilcox.test(plc_bc, mu = 0, alternative = "greater")   # placebo 組呢？0.03125


# 計算每個病人 6M 與 Baseline 的「重心位移」。 在 PCoA 座標軸上計算： Shift = (PC1_6M - PC1_Baseline)
# 然後比較：
# FMT 組的位移是否都往同一個方向（例如 PC1 增加）？
# Placebo 組的位移是否雜亂無章（有的增加有的減少）？ 
# 如果 FMT 組的位移方向一致，而 Placebo 亂跳，這就能完美解釋為什麼「變動量」差不多，但只有 FMT 組症狀改善！
# 從 PCoA scores 取 PC1 / PC2
pcoa_scores <- data.frame(
  ord_pcoa$vectors[, 1:2],
  sample_data(ps_filt)
) %>%
  dplyr::rename(PC1 = Axis.1, PC2 = Axis.2)

# 計算每個病人 6M - Baseline 的位移
shift_df <- pcoa_scores %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  select(PatientID, Treatment, Timepoint, PC1, PC2) %>%
  pivot_wider(names_from = Timepoint,
              values_from = c(PC1, PC2)) %>%
  dplyr::rename(
    PC1_bl = `PC1_Baseline`,
    PC1_6m = `PC1_6M post-FMT`,
    PC2_bl = `PC2_Baseline`,
    PC2_6m = `PC2_6M post-FMT`
  ) %>%
  mutate(
    dPC1 = PC1_6m - PC1_bl,
    dPC2 = PC2_6m - PC2_bl
  ) %>% 
  filter(!is.na(PC1_bl))

# 方向一致性：FMT 組 dPC1 是否都同號？
shift_df %>%
  group_by(Treatment) %>%
  summarise(
    mean_dPC1  = mean(dPC1),
    sd_dPC1    = sd(dPC1),
    all_same_direction = all(dPC1 > 0) | all(dPC1 < 0),   # 方向一致性
    n_positive = sum(dPC1 > 0),
    n_negative = sum(dPC1 < 0)
  )

# Wilcoxon：兩組的 dPC1 是否有差異
wilcox.test(dPC1 ~ Treatment, data = shift_df, exact = FALSE)
wilcox.test(dPC2 ~ Treatment, data = shift_df, exact = FALSE)

# 與 SALT score 改善量的相關性
salt_change2 <- meta_raw %>%
  filter(Timepoint %in% c("Baseline", "2M post-FMT")) %>%
  select(PatientID, Timepoint, SALT_score) %>%
  pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
  dplyr::rename(SALT_bl = Baseline, SALT_2m = `2M post-FMT`) %>%
  mutate(dSALT = SALT_2m - SALT_bl) %>%    # 負值 = 改善
  filter(!is.na(dSALT))

shift_salt <- left_join(shift_df, salt_change, by = "PatientID")

cor.test(shift_salt$dPC1, shift_salt$dSALT,
         method = "spearman")
# M2 failed / M6 worked but sample數量應該要再增加 不然pval不顯著
salt_change6 <- meta_raw %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  select(PatientID, Timepoint, SALT_score) %>%
  pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
  dplyr::rename(SALT_bl = Baseline, SALT_6m = `6M post-FMT`) %>%
  mutate(dSALT = SALT_6m - SALT_bl) %>%    # 負值 = 改善
  filter(!is.na(dSALT))

wilcox.test(salt_change2$dSALT, salt_change6$dSALT, paired = TRUE)

# viz arror pc
# 箭頭圖：每個病人在 PCoA 空間的位移向量
p_shift_arrow <- ggplot() +
  # baseline 點
  geom_point(data = pcoa_scores %>% filter(Timepoint == "Baseline"),
             aes(x = PC1, y = PC2, colour = Treatment),
             size = 3, shape = 1, stroke = 1.2) +
  # 箭頭：baseline → 6M
  geom_segment(data = shift_df,
               aes(x = PC1_bl, y = PC2_bl,
                   xend = PC1_6m, yend = PC2_6m,
                   colour = Treatment),
               arrow = arrow(length = unit(0.015, "npc"),
                             type = "closed"),
               linewidth = 0.9, alpha = 0.8) +
  geom_text_repel(data = shift_df,
                  aes(x = PC1_6m, y = PC2_6m,
                      label = PatientID, colour = Treatment),
                  size = 3, segment.colour = NA) +
  scale_colour_manual(values = col_treatment) +
  labs(title    = "PCoA displacement: Baseline → 6M",
       subtitle = "Circle = baseline, arrowhead = 6M post-FMT",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme_fmt

# dPC1 分佈：bar + 參考線
p_dPC1 <- ggplot(shift_df,
                 aes(x = reorder(PatientID, dPC1),
                     y = dPC1, fill = Treatment)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = col_treatment) +
  labs(title = "PC1 shift (6M − Baseline) per patient",
       x = NULL, y = "ΔPC1") +
  theme_fmt +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# dPC1 vs dSALT scatter
p_shift_salt <- ggplot(shift_salt,
                       aes(x = dPC1, y = -dSALT, colour = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              fill = "grey85", linewidth = 0.8) +
  geom_text_repel(aes(label = PatientID), size = 3,
                  segment.colour = NA) +
  scale_colour_manual(values = col_treatment) +
  labs(title    = "PC1 shift vs SALT score change",
       subtitle = "Negative dSALT = improvement",
       x = "ΔPC1 (6M − Baseline)",
       y = "ΔSALT (6M − Baseline) x (-1)") +
  theme_fmt


library(vegan)

# Step 1：genus level relative abundance matrix
ps_genus_rel <- tax_glom(ps_filt, taxrank = "Genus") %>%
  transform_sample_counts(function(x) x / sum(x))

genus_mat <- as.data.frame(t(otu_table(ps_genus_rel)))   # rows = samples
tax_df    <- as.data.frame(tax_table(ps_genus_rel))

# 欄位名稱改成 Genus 名稱（方便看圖）
colnames(genus_mat) <- tax_df$Genus

# Step 2：確認樣本順序對齊
pcoa_scores_mat <- pcoa_scores %>% filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>% select(PC1, PC2) %>% as.matrix()
genus_mat       <- genus_mat[rownames(pcoa_scores_mat), ]   # 強制對齊

# Step 3：envfit
set.seed(42)
fit <- envfit(pcoa_scores_mat, genus_mat, permutations = 999, na.rm = TRUE)

# Step 4：抽出顯著結果
fit_scores <- as.data.frame(scores(fit, display = "vectors")) %>%
  tibble::rownames_to_column("Genus") %>%
  mutate(
    r2   = fit$vectors$r,
    pval = fit$vectors$pvals
  ) %>%
  filter(pval < 0.05) %>%
  arrange(desc(r2))

cat("顯著 taxa (p < 0.05):", nrow(fit_scores), "\n")
print(fit_scores)

#  視覺化：疊加在 PCoA 箭頭圖上
# 只取 r2 前 10 的 genus 畫箭頭
top_taxa <- fit_scores %>% slice_max(r2, n = 10)

# 縮放因子（讓箭頭長度跟 PCoA 座標尺度相符）
scale_factor <- 0.3

p_envfit <- p_shift_arrow +   # 你原本的 PCoA 圖
  geom_segment(
    data = top_taxa,
    aes(x = 0, y = 0,
        xend = PC1 * scale_factor,
        yend = PC2 * scale_factor),
    arrow     = arrow(length = unit(0.012, "npc"), type = "closed"),
    colour    = "grey20",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = top_taxa,
    aes(x = PC1 * scale_factor,
        y = PC2 * scale_factor,
        label = Genus),
    colour   = "grey10",
    size     = 3,
    fontface = "italic",
    segment.colour = NA,
    inherit.aes = FALSE
  ) +
  labs(subtitle = "Arrows: envfit vectors (p < 0.05, top 10 by r²) Genus")



# =============================================================================
#  Correlation Heatmap：top envfit genera × SALT score change
# =============================================================================

# Step 1：取出 top 10 envfit genus 名稱
top_genera_names <- top_taxa$Genus   # 來自pcoa結果 top_taxa dataframe 直接取

# Step 2：建立 genus 豐度變化量 (6M - Baseline)
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
  rename(bl = Baseline, m6 = `6M post-FMT`) %>%
  mutate(delta_abund = m6 - bl) %>%           # 豐度變化量
  select(PatientID, Genus, delta_abund) %>%
  pivot_wider(names_from = Genus, values_from = delta_abund) %>% 
  filter(!is.na(Faecalibacterium))

# Step 3：SALT score 變化量 (6M - Baseline)
salt_change <- meta_raw %>%
  filter(Timepoint %in% c("Baseline", "6M post-FMT")) %>%
  select(PatientID, Timepoint, SALT_score) %>%
  pivot_wider(names_from = Timepoint, values_from = SALT_score) %>%
  rename(SALT_bl = Baseline, SALT_6m = `6M post-FMT`) %>%
  mutate(dSALT = SALT_6m - SALT_bl) %>%   # 負值 = 改善
  filter(!is.na(dSALT))

# Step 4：合併，只保留有完整資料的病人
cor_input <- genus_change %>%
  left_join(salt_change %>% select(PatientID, dSALT), by = "PatientID") %>%
  tibble::column_to_rownames("PatientID")

# Step 5：Spearman 相關性計算（每個 genus vs dSALT）
cor_results <- top_genera_names %>%
  purrr::map_dfr(function(g) {
    test <- cor.test(cor_input[[g]], -cor_input$dSALT,
                     method = "spearman", exact = FALSE)
    tibble(
      Genus = g,
      rho   = round(test$estimate, 3),
      pval  = round(test$p.value,  3),
      sig   = case_when(
        test$p.value < 0.05  ~ "*",
        test$p.value < 0.1   ~ "†",
        TRUE                 ~ ""
      )
    )
  }) %>%
  arrange(rho)

cat("Spearman correlation (Δgenus abundance vs ΔSALT):\n")
print(cor_results)

tmp1 <- meta_raw %>% 
  select(PatientID, Treatment) %>% 
  distinct() %>% 
  left_join(.,tibble::rownames_to_column(cor_input, "PatientID"), by = "PatientID") %>% 
  filter(!is.na(Faecalibacterium))

plot_list <- list() # 建立一個空的 list 來存圖

for (i in top_genera_names) {
  y_lab <- paste0("Δ", i)
  
  p <- tmp1 %>% 
    ggplot(aes(x = dSALT, y = !!sym(i), color = Treatment)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) + 
    theme_bw() +
    theme(legend.position = "none") + # 先關掉 legend，最後再統一顯示
    labs(y = y_lab, x = "ΔSALT")+
    scale_color_manual(values = col_treatment)
  
  plot_list[[i]] <- p # 將圖存入 list
}


combined_plot <- wrap_plots(plot_list, ncol = 5) + 
  plot_layout(guides = "collect") & # 統一收集 legend
  theme(legend.position = "bottom")

# 顯示圖片
print(combined_plot)

# Step 6：準備 heatmap 矩陣
# 行 = 病人，欄 = genus，值 = delta_abund（z-score 標準化方便視覺比較）
heatmap_mat <- cor_input %>%
  select(all_of(top_genera_names)) %>%
  scale() %>%                          # z-score per genus
  as.data.frame()

# 加入 dSALT 欄（也標準化）
heatmap_mat$dSALT <- scale(cor_input$dSALT)[, 1]

# 加入 Treatment 標注
heatmap_mat <- heatmap_mat %>%
  tibble::rownames_to_column("PatientID") %>%
  left_join(
    meta_raw %>% filter(Timepoint == "Baseline") %>%
      select(PatientID, Treatment) %>% distinct(),
    by = "PatientID"
  ) %>% 
  filter(!is.na(Treatment))

# Step 7：長格式 for ggplot
heatmap_long <- heatmap_mat %>%
  pivot_longer(cols = c(all_of(top_genera_names), "dSALT"),
               names_to = "Variable", values_to = "z") %>%
  mutate(
    # 欄位分組：genus vs clinical
    var_type = if_else(Variable == "dSALT", "Clinical", "Genus"),
    # 加入 rho 和顯著性標注（只對 genus 欄）
    label = if_else(
      Variable %in% cor_results$Genus,
      cor_results$sig[match(Variable, cor_results$Genus)],
      ""
    ),
    # Genus 欄位按照 rho 排序
    Variable = factor(Variable,
                      levels = c("dSALT",
                                 cor_results %>% arrange(rho) %>% pull(Genus)))
  )

# Step 8：畫圖
p_cor_heatmap <- ggplot(heatmap_long,
                        aes(x = Variable, y = PatientID, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  # 顯著性標記
  geom_text(aes(label = label), size = 4, colour = "grey30", fontface = "bold") +
  # 分隔線區分 genus vs SALT
  # geom_vline(xintercept = 1.5, colour = "grey30", linewidth = 1) +
  # 病人按 Treatment 分組
  facet_grid(Treatment ~ var_type, scales = "free", space = "free") +
  scale_fill_gradient2(
    low      = "#2166AC",   # 負值（藍）= 豐度下降 / SALT 改善
    mid      = "white",
    high     = "#D6604D",   # 正值（紅）= 豐度上升 / SALT 惡化
    midpoint = 0,
    name     = "z-score"
  ) +
  labs(
    title    = "Δ Genus abundance vs Δ SALT score (3M − Baseline)",
    subtitle = "* p<0.05  † p<0.1  |  Columns sorted by Spearman ρ with ΔSALT",
    x        = NULL,
    y        = NULL
  ) +
  theme_fmt +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, face = "italic"),
    strip.text.x = element_text(face = "bold"),
    legend.position = "right"
  )

# Step 9：加上 rho 數值 bar（側邊小圖）
p_rho_bar <- ggplot(cor_results,
                    aes(x = reorder(Genus, rho), y = rho,
                        fill = rho > 0)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0(rho, sig)),
            hjust = if_else(cor_results$rho > 0, -0.1, 1.1),
            size  = 3) +
  geom_hline(yintercept = 0, colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "#2166AC")) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  labs(title = "Spearman ρ (Δabundance vs ΔSALT)",
       x = NULL, y = "ρ") +
  theme_fmt +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))

# 合併兩張圖
p_cor_panel <- p_cor_heatmap + p_rho_bar +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "Genus–SALT correlation panel",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

save_fig(p_cor_panel, "fig6_genus_salt_correlation.pdf", w = 14, h = 7)




# =============================================================================
#  SECTION 3 — Taxa Composition
# =============================================================================

cat("\n── Section 3: Taxa composition ─────────────────────────────────────────\n")

plot_constitution_stack_bar <- function(ps_obj = ps_filt, level = "Genus", top_N = 10 ){
  
  # Aggregate to genus level
  ps_glom <- tax_glom(ps_obj, taxrank = level, NArm = TRUE)
  
  # Relative abundance
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)

  # Top N genera across all samples
  topN <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:top_N]
  
  # Melt to long format
  name2 = paste0(level, "2")
  ps_milt <- psmelt(ps_rel) %>%
    mutate(!!sym(name2) := if_else(OTU %in% topN, as.character(!!sym(level)), "Other")) %>%
    group_by(SampleID, !!sym(name2), Timepoint, Treatment, PatientID) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  lev_summary <- ps_milt %>%
    group_by(Treatment, Timepoint, !!sym(name2)) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop")
  
  # 顏色設定 (安全機制：確保 top_N 在 3~12 之間，或改用 colorRampPalette 彈性生成)
  pal_n <- max(3, min(top_N, 12)) 
  unique_groups <- unique(ps_milt[[name2]][ps_milt[[name2]] != "Other"])
  
  genus_cols <- c(
    setNames(
      colorRampPalette(brewer.pal(pal_n, "Paired"))(length(unique_groups)), 
      unique_groups
    ),
    Other = "grey80"
  )
  
  gn <- lev_summary[[name2]] %>% unique()
  gn1 <- ! gn %in% c("Other", NA)
  lev_summary[[name2]] <- factor(lev_summary[[name2]], levels = c(gn[gn1], "Other", NA))
  
  
  p_barplot <- ggplot(lev_summary,
                      aes(x = Timepoint, y = mean_abund, fill = !!sym(name2))) +
    geom_col(position = "stack", width = 0.75) +
    facet_wrap(~ Treatment) +
    scale_fill_manual(values = genus_cols) +
    scale_y_continuous(labels = label_percent(scale = 1)) +
    labs(title = sprintf("%s-level relative abundance", level),
         x = NULL, y = "Mean relative abundance (%)", fill = level) +
    theme_fmt +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  list(
    bar_plt = p_barplot,
    ps_milt = ps_milt,
    lev_summary = lev_summary
  )
}

p_barplot_genus <- plot_constitution_stack_bar(ps_obj = ps_filt, level = "Genus", top_N = 10)
p_barplot_phylum <- plot_constitution_stack_bar(ps_obj = ps_filt, level = "Phylum", top_N = 10)

# remove top 3
tmp =p_barplot_phylum$lev_summary %>% filter(! Phylum2 %in% c("Actinobacteria","Bacteroidetes","Firmicutes"))
ggplot(tmp, aes(x = Timepoint, y = mean_abund, fill = Phylum2)) + 
  geom_col(position = "stack", width = 0.75)+
  facet_wrap(~ Treatment)+
  theme_fmt


# Firmicutes / Bacteroidetes ratio
ps_phylum <- tax_glom(ps_filt, taxrank = "Phylum", NArm = T)
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x) * 100)

fb_df <- psmelt(ps_phylum_rel) %>%
  filter(Phylum %in% c("Firmicutes", "Bacteroidetes")) %>%
  select(SampleID, Phylum, Abundance, Timepoint, Treatment, PatientID) %>%
  pivot_wider(names_from = Phylum, values_from = Abundance) %>%
  mutate(FB_ratio = Firmicutes / (Bacteroidetes + 0.01))
retained_id <- filter(fb_df, FB_ratio<=100) %>% group_by(PatientID) %>% count() %>% filter(n==3) %>% pull(PatientID)
retain_normal_fb <- filter(fb_df, FB_ratio<=100, PatientID %in% retained_id)
p_fb <- ggplot(retain_normal_fb, aes(x = Timepoint, y = FB_ratio,
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
  geom_text(data = filter(retain_normal_fb, Timepoint == "6M post-FMT"),aes(label = PatientID), 
            position = position_jitterdodge())+
  theme_fmt



tp <- ggplot(retain_normal_fb, aes(x = Timepoint, y = FB_ratio, fill = Treatment ))+
  geom_boxplot(outliers = F)+
  scale_fill_manual(values = col_treatment)+
  theme_fmt

trt <- ggplot(retain_normal_fb, aes(x = Treatment, y = FB_ratio, fill = Treatment ))+
  geom_boxplot(outliers = F)+
  scale_fill_manual(values = col_treatment)+
  theme_fmt

tp + tp +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = sprintf("Firmicutes / Bacteroidetes ratio for different features"),
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# =============================================================================
#  SECTION 4 — SALT Score Correlation
# =============================================================================

cat("\n── Section 4: SALT score correlation ───────────────────────────────────\n")

# Join alpha diversity with SALT scores
salt_alpha <- alpha_df %>% select(SampleID, Shannon, Chao1, Timepoint, Treatment,
                                   PatientID, SALT_score, Age)

# Scatter: Shannon vs SALT coloured by time point $ smooth line is based ont Treatment group????
p_salt_shannon <- ggplot(salt_alpha,
                          aes(x = Shannon, y = SALT_score,
                              colour = Timepoint, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              colour = "grey40", fill = "grey85") +
  #geom_text_repel(aes(label = PatientID), size = 2.5, max.overlaps = 6) +
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
  stat_summary(data = salt_alpha, aes(group = Treatment), fun = mean,
               geom = "line", linewidth = 1.5) +
  stat_summary(data = salt_alpha, aes(group = Treatment, fill = Treatment), fun.data = mean_se,
               geom = "ribbon", alpha = 0.15, colour = NA) +
  scale_colour_manual(values = col_treatment) +
  scale_fill_manual(values   = col_treatment) +
  labs(title = "SALT score trajectory",
       x = "", y = "SALT score", colour = "Treatment") +
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
ps_m3   <- subset_samples(ps_filt, Timepoint == "3M post-FMT")
res_m3  <- run_deseq2(ps_m3, "Treatment", "placebo", tag = "FMT vs Placebo @ month3")

# 5b. Within FMT: baseline vs month 3
ps_fmt  <- subset_samples(ps_filt, Treatment == "FMT")
res_fmt <- run_deseq2(ps_fmt, "Timepoint", "Baseline", tag = "FMT baseline→month3")

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
#  SECTION 6 — Save all plots
# =============================================================================

cat("\n── Saving figures ──────────────────────────────────────────────────────\n")

save_fig <- function(plot, filename, w = 10, h = 6) {
  ggsave(filename, plot = plot, width = w, height = h,
         dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", filename))
}

save_fig(p_alpha_combined, "fig1_alpha_diversity.pdf", w = 12, h = 5)
save_fig(p_pcoa_treatment, "fig2a_pcoa_treatment.pdf")
save_fig(p_pcoa_time,      "fig2b_pcoa_timepoint.pdf")
save_fig(p_bc_change,      "fig2c_bc_within_change.pdf", w = 8, h = 5)
save_fig(p_barplot,        "fig3_genus_barplot.pdf", w = 12, h = 6)
save_fig(p_fb,             "fig3b_firmicutes_bacteroidetes.pdf", w = 10, h = 5)
save_fig(p_salt_shannon,   "fig4a_salt_shannon_scatter.pdf")
save_fig(p_salt_traj,      "fig4b_salt_trajectory.pdf")
save_fig(p_volcano_m3,     "fig5a_volcano_fmt_vs_placebo_m3.pdf")
save_fig(p_volcano_fmt,    "fig5b_volcano_fmt_time.pdf")
save_fig(p_lfc_m3,         "fig5c_lfc_bar_m3.pdf", w = 10, h = 6)
save_fig(p_lfc_fmt,        "fig5d_lfc_bar_fmt.pdf", w = 10, h = 6)

# Combined poster panel
poster <- (p_alpha_combined) /
  (p_pcoa_treatment | p_salt_traj) /
  (p_lfc_m3 | p_barplot) +
  plot_annotation(
    title    = "Alopecia oral FMT — Gut microbiome dynamics",
    subtitle = "n=10 | 3 time points | FMT vs placebo",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )
save_fig(poster, "fig_poster_panel.pdf", w = 16, h = 18)

# Save DESeq2 results tables
write.csv(res_m3,  "deseq2_fmt_vs_placebo_month3.csv", row.names = FALSE)
write.csv(res_fmt, "deseq2_fmt_baseline_vs_month3.csv", row.names = FALSE)
write.csv(alpha_df, "alpha_diversity_table.csv", row.names = FALSE)

cat("\nAll done! Output files written to working directory.\n")


# =============================================================================
#  QUICK-START: simulate data to verify pipeline runs
# =============================================================================
# Uncomment the block below to generate synthetic data and test the script
# without your real files.
#
# source("simulate_data.R")   # creates asv_table.csv, taxonomy.csv, metadata.csv
