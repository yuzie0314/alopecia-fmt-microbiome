# =============================================================================
#  SECTION 1 — Alpha Diversity
#  Richness  : Observed ASVs, Chao1, ACE
#  Diversity : Shannon
#  Evenness  : Pielou's J (Shannon / log(Observed))
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

cat("\n── Section 1: Alpha diversity ──────────────────────────────────────────\n")

alpha_df <- estimate_richness(ps_filt,
                              measures = c("Observed", "Chao1", "ACE", "Shannon")) %>%
  tibble::rownames_to_column("SampleID") %>%
  mutate(Pielou_J = Shannon / log(Observed)) %>%
  left_join(meta_raw %>% select(SampleID, PatientID, Timepoint, Treatment,
                                SALT_score, Age, Race, fast_responder),
            by = "SampleID") %>% 
  mutate(Treatment = relevel(factor(Treatment), ref = "placebo"))

# Summary statistics
alpha_summary <- alpha_df %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    Observed_mean = round(mean(Observed),  1),
    Observed_sd   = round(sd(Observed),    1),
    Chao1_mean    = round(mean(Chao1),     1),
    Chao1_sd      = round(sd(Chao1),       1),
    ACE_mean      = round(mean(ACE),       1),
    ACE_sd        = round(sd(ACE),         1),
    Shannon_mean  = round(mean(Shannon),   3),
    Shannon_sd    = round(sd(Shannon),     3),
    PielouJ_mean  = round(mean(Pielou_J),  3),
    PielouJ_sd    = round(sd(Pielou_J),    3),
    .groups = "drop"
  )
cat("Alpha diversity summary:\n")
print(alpha_summary)

# ── Helper：共用的 trajectory plot function ────────────────────────────────────

plot_alpha_traj <- function(df, yvar, ytitle) {
  ggplot(df, aes(x = Timepoint, y = .data[[yvar]],
                 colour = Treatment, group = PatientID)) +
    geom_line(alpha = 0.35, linewidth = 0.6) +
    geom_point(alpha = 0.6, size = 2) +
    stat_summary(aes(group = Treatment), fun = mean,
                 geom = "line",  linewidth = 1.4) +
    stat_summary(aes(group = Treatment), fun = mean,
                 geom = "point", size = 4, shape = 18) +
    scale_colour_manual(values = col_treatment) +
    labs(title = ytitle, x = NULL, y = ytitle, colour = "Treatment") +
    theme_fmt+ 
    geom_text(data = filter(df, Timepoint == "6M post-FMT", fast_responder == FALSE),aes(label =PatientID ), position = position_dodge(width =0.2))+
    geom_text(data = filter(df, Timepoint == "Baseline", fast_responder == TRUE),aes(label =PatientID ), position = position_dodge(width =0.3))
}

# ── 1a. Individual trajectory plots ──────────────────────────────────────────

p_shannon  <- plot_alpha_traj(alpha_df, "Shannon",   "Shannon index")
p_chao1    <- plot_alpha_traj(alpha_df, "Chao1",     "Chao1 (richness)")
p_ace      <- plot_alpha_traj(alpha_df, "ACE",       "ACE (richness)")
p_observed <- plot_alpha_traj(alpha_df, "Observed",  "Observed ASVs (richness)")
p_pielou   <- plot_alpha_traj(alpha_df, "Pielou_J",  "Pielou's J (evenness)")

# ── 1b. Combined panel：2 × 3 layout ─────────────────────────────────────────
#  Row 1: Richness (Observed | Chao1 | ACE)
#  Row 2: Diversity + Evenness (Shannon | Pielou_J | blank)

p_alpha_combined <- (p_observed | p_chao1 | p_ace) /
  (p_shannon | p_pielou | plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Alpha diversity — richness, diversity & evenness",
    subtitle = "Bold diamond = group mean",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  ) &
  theme(legend.position = "bottom")

# ── 1c. Linear mixed models ───────────────────────────────────────────────────

run_lmm_alpha <- function(metric) {
  lmer(as.formula(
    paste(metric, "~ Timepoint * Treatment + Age + (1 | PatientID)")),
    data = alpha_df, REML = TRUE)
}

lmm_shannon  <- run_lmm_alpha("Shannon")
lmm_pielou   <- run_lmm_alpha("Pielou_J")
lmm_observed <- run_lmm_alpha("Observed")
lmm_chao1    <- run_lmm_alpha("Chao1")

cat("\nLMM Shannon — fixed effects:\n");  print(summary(lmm_shannon)$coefficients)
cat("\nLMM Pielou J — fixed effects:\n"); print(summary(lmm_pielou)$coefficients)
cat("\nLMM Observed — fixed effects:\n"); print(summary(lmm_observed)$coefficients)
cat("\nLMM Chao1 — fixed effects:\n");    print(summary(lmm_chao1)$coefficients)

# Pairwise EMM for each metric
for (lmm_obj in list(Shannon  = lmm_shannon,
                     Pielou_J = lmm_pielou,
                     Observed = lmm_observed,
                     Chao1    = lmm_chao1)) {
  idx <- deparse(formula(lmm_obj)[[2]])
  emm <- emmeans(lmm_obj, ~ Treatment | Timepoint) # 可以改Timepoint | Treatment
  cat(sprintf("\nEMM pairwise (%s within treatment):\n",
              idx))
  print(pairs(emm))
  y_lab <- paste0(idx, "Index (Estimated Mean)")
  p <- ggplot(as.data.frame(emm), aes(x = Timepoint, y = emmean, color = Treatment, group = Treatment)) +
    geom_line(position = position_dodge(0.3)) +
    geom_point(position = position_dodge(0.3), size =3, alpha = 0.6) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(0.3))+
    scale_color_manual(values = col_treatment) +
    labs(title = sprintf("Estimated Marginal Means of Alpha Diversity (%s)",idx),
         y = y_lab,
         x = "Timepoint") +
    theme(legend.position = "bottom")+theme_fmt
  print(p)
}


# 2. 開始繪圖
ggplot(as.data.frame(emm) , aes(x = Timepoint, y = emmean, color = Treatment, group = Treatment)) +
  # 畫出點（平均值）
  geom_point(size = 3, position = position_dodge(0.3)) +
  # 畫出誤差棒 (95% CI)
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.3)) +
  # 畫出折線觀察趨勢
  geom_line(position = position_dodge(0.3)) +
  # 分面（依處理組分開看，或是畫在同一張）
  facet_wrap(~Treatment, scales = "free") +
  # 美化
  theme_fmt +
  labs(title = "Estimated Marginal Means of Alpha Diversity",
       y = "Chao1 Index (Estimated Mean)",
       x = "Timepoint") +
  theme(legend.position = "bottom")+
  coord_flip()




# ── 1d. Wilcoxon per metric（baseline vs 6M，within each treatment）────────────
#  n=5 per group → 非參數更保守但適當

cat("\nWilcoxon: Baseline vs 6M post-FMT within each treatment group\n") # 可以看2M: 在treatment group下都沒有差異

metrics_list <- c("Shannon", "Pielou_J", "Observed", "Chao1", "ACE")

wilcox_alpha <- c("FMT", "placebo") %>%
  purrr::map_dfr(function(trt) {
    sub <- alpha_df %>% filter(Treatment == trt) %>% filter(PatientID != 'GM03')
    bl  <- sub %>% filter(Timepoint == "Baseline") %>% arrange(PatientID)
    m6  <- sub %>% filter(Timepoint == "6M post-FMT") %>% arrange(PatientID)
    
    metrics_list %>%
      purrr::map_dfr(function(m) {
        w <- wilcox.test(bl[[m]], m6[[m]], paired = TRUE, exact = FALSE)
        tibble(Treatment = trt, Metric = m,
               W = w$statistic, pval = round(w$p.value, 4))
      })
  })

cat("\n"); print(wilcox_alpha)


# ── 1e. Non-parametric tests ──────────────────────────────────────────────────

cat("\n── 1e: Non-parametric tests ────────────────────────────────────────────\n")

# 1e-i. Mann-Whitney U：FMT vs Placebo 在每個時間點（組間比較）# 組間沒有差異
cat("\nMann-Whitney U: FMT vs Placebo per timepoint\n")

mwu_results <- expand.grid(
  Timepoint = levels(alpha_df$Timepoint),
  Metric    = metrics_list,
  stringsAsFactors = FALSE
) %>%
  purrr::pmap_dfr(function(Timepoint, Metric) {
    sub  <- alpha_df %>% filter(Timepoint == !!Timepoint)
    fmt  <- sub %>% filter(Treatment == "FMT")     %>% pull(!!Metric)
    plc  <- sub %>% filter(Treatment == "placebo") %>% pull(!!Metric)
    w    <- wilcox.test(fmt, plc, exact = FALSE)
    tibble(Timepoint, Metric,
           W    = w$statistic,
           pval = round(w$p.value, 4))
  })

print(mwu_results)

# 1e-ii. Kruskal-Wallis：三個時間點整體差異（within each treatment）
cat("\nKruskal-Wallis: Timepoint effect within each treatment group\n")

kw_results <- expand.grid(
  Treatment = c("FMT", "placebo"),
  Metric    = metrics_list,
  stringsAsFactors = FALSE
) %>%
  purrr::pmap_dfr(function(Treatment, Metric) {
    sub <- alpha_df %>% filter(Treatment == !!Treatment)
    kw  <- kruskal.test(as.formula(paste(Metric, "~ Timepoint")), data = sub)
    tibble(Treatment, Metric,
           chi2 = round(kw$statistic, 3),
           df   = kw$parameter,
           pval = round(kw$p.value, 4))
  })

print(kw_results)

# 1e-iii. Friedman Test：配對資料（同一病人三個時間點），within each treatment
cat("\nFriedman Test: repeated measures across timepoints\n")

friedman_results <- expand.grid(
  Treatment = c("FMT", "placebo"),
  Metric    = metrics_list,
  stringsAsFactors = FALSE
) %>%
  purrr::pmap_dfr(function(Treatment, Metric) {
    sub <- alpha_df %>%
      filter(Treatment == !!Treatment) %>%
      select(PatientID, Timepoint, value = !!Metric) %>%
      # 必須是完整的 balanced block
      group_by(PatientID) %>%
      filter(n() == 3) %>%
      ungroup() %>%
      arrange(PatientID, Timepoint)
    
    mat <- sub %>%
      pivot_wider(names_from = Timepoint, values_from = value) %>%
      select(-PatientID) %>%
      as.matrix()
    
    ft <- friedman.test(mat)
    tibble(Treatment, Metric,
           chi2 = round(ft$statistic, 3),
           df   = ft$parameter,
           pval = round(ft$p.value, 4))
  })

print(friedman_results)

# Post-hoc Dunn test after significant Friedman（pval < 0.1）
if (requireNamespace("dunn.test", quietly = TRUE)) {
  library(dunn.test)
  sig_friedman <- friedman_results %>% filter(pval < 0.1)
  if (nrow(sig_friedman) > 0) {
    cat("\nDunn post-hoc for significant Friedman results:\n")
    sig_friedman %>%
      purrr::pwalk(function(Treatment, Metric, ...) {
        cat(sprintf("\n  %s | %s\n", Treatment, Metric))
        sub <- alpha_df %>%
          filter(Treatment == !!Treatment) %>%
          arrange(PatientID, Timepoint)
        dunn.test(sub[[Metric]], sub$Timepoint,
                  method = "bh", alpha = 0.05, list = TRUE)
      })
  }
}


# ── 1f. Multiple testing correction (BH / FDR) ────────────────────────────────

cat("\n── 1f: Multiple testing correction (Benjamini-Hochberg) ───────────────\n")

# 收集所有 p 值：Wilcoxon paired (BL→6M) + Mann-Whitney + Kruskal-Wallis
all_pvals <- bind_rows(
  wilcox_alpha  %>% mutate(test = "Wilcoxon_paired_BL_6M"),
  mwu_results   %>% rename(Treatment = Timepoint) %>% mutate(test = "MannWhitney_FMTvsPlacebo"),
  kw_results    %>% mutate(test = "KruskalWallis_timepoint"),
  friedman_results %>% mutate(test = "Friedman_repeated")
) %>%
  select(test, Treatment, Metric, pval) %>%
  mutate(
    pval_BH  = p.adjust(pval, method = "BH"),        # Benjamini-Hochberg (FDR)
    pval_bon = p.adjust(pval, method = "bonferroni"), # Bonferroni（保守）
    sig_BH   = case_when(
      pval_BH < 0.05 ~ "**",
      pval_BH < 0.10 ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  arrange(pval_BH)

cat("All tests after BH correction (sorted by adj. p):\n")
print(all_pvals, n = 40)

write.csv(all_pvals, "outputs/alpha_nonparametric_tests_BH.csv", row.names = FALSE)


# ── 1g. Delta analysis：Δ metric (timepoint - baseline) ──────────────────────

cat("\n── 1g: Delta analysis (change from baseline) ───────────────────────────\n")

alpha_delta <- alpha_df %>%
  select(PatientID, Treatment, Timepoint, all_of(metrics_list)) %>%
  pivot_longer(cols = all_of(metrics_list),
               names_to = "Metric", values_to = "value") %>%
  pivot_wider(names_from = Timepoint, values_from = value) %>%
  rename(bl = Baseline, m2 = `2M post-FMT`, m6 = `6M post-FMT`) %>%
  mutate(
    delta_2M = m2 - bl,
    delta_6M = m6 - bl
  )

# Wilcoxon: FMT vs placebo on delta_6M or delta_2M
delta_wilcox <- alpha_delta %>%
  group_by(Metric) %>%
  summarise(
    W    = wilcox.test(delta_6M[Treatment == "FMT"],
                       delta_6M[Treatment == "placebo"],
                       exact = FALSE)$statistic,
    pval = wilcox.test(delta_6M[Treatment == "FMT"],
                       delta_6M[Treatment == "placebo"],
                       exact = FALSE)$p.value %>% round(4),
    .groups = "drop"
  ) %>%
  mutate(pval_BH = p.adjust(pval, method = "BH") %>% round(4))

cat("Wilcoxon on Δ (6M − Baseline): FMT vs placebo\n")
print(delta_wilcox)

# Delta trajectory plot
alpha_delta_long <- alpha_delta %>%
  select(PatientID, Treatment, Metric, delta_2M, delta_6M) %>%
  pivot_longer(cols = c(delta_2M, delta_6M),
               names_to = "Transition", values_to = "delta") %>%
  mutate(Transition = recode(Transition,
                             delta_2M = "BL→2M",
                             delta_6M = "BL→6M"))

p_delta <- ggplot(alpha_delta_long,
                  aes(x = Transition, y = delta,
                      colour = Treatment, group = PatientID)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(alpha = 0.4, linewidth = 0.6) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_text(data = filter(alpha_delta_long, Transition == "BL→2M"), aes(label =PatientID ), position = position_dodge(width =0.1))+
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "line",  linewidth = 1.4) +
  stat_summary(aes(group = Treatment), fun = mean,
               geom = "point", size = 4, shape = 18) +
  scale_colour_manual(values = col_treatment) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  labs(title    = "Change from baseline (Δ metric)",
       subtitle = "Positive = increase from baseline",
       x = NULL, y = "Δ value") +
  theme_fmt +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


# ── 1h. Correlation: alpha index vs SALT score ────────────────────────────────

cat("\n── 1h: Alpha index vs SALT score correlation ───────────────────────────\n")
# 可以替換Timepoint / Treatment
cor_alpha_salt <- expand.grid(
  Metric    = metrics_list,
  Timepoint = levels(alpha_df$Timepoint),
  stringsAsFactors = FALSE
) %>%
  purrr::pmap_dfr(function(Metric, Timepoint) {
    sub <- alpha_df %>% filter(Timepoint == !!Timepoint)
    ct  <- cor.test(sub[[Metric]], sub$SALT_score,
                    method = "spearman", exact = FALSE)
    tibble(Metric, Timepoint,
           rho  = round(ct$estimate, 3),
           pval = round(ct$p.value,  4))
  }) %>%
  mutate(pval_BH = p.adjust(pval, method = "BH") %>% round(4),
         sig     = case_when(pval_BH < 0.05 ~ "*", pval_BH < 0.1 ~ "†", TRUE ~ "")) %>% 
  arrange(pval)

cat("Spearman: alpha metric vs SALT score\n")
print(cor_alpha_salt)

# Heatmap of rho values
p_cor_alpha_salt <- ggplot(cor_alpha_salt,
                           aes(x = Timepoint, y = Metric, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), sig)),
            size = 3.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1), name = "Spearman ρ"
  ) +
  labs(title    = "Alpha diversity vs SALT score — Spearman ρ",
       subtitle = "* BH-adj p<0.05  † BH-adj p<0.1",
       x = NULL, y = NULL) +
  theme_fmt +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

alpha_df %>% 
  select(Shannon, Pielou_J, SALT_score, Treatment) %>% 
  pivot_longer(cols = c("Shannon", "Pielou_J"),names_to = "index",values_to = "value") %>% 
  ggplot(aes(x = value, y = SALT_score, fill = Treatment,color = Treatment))+
  geom_point(size = 3, alpha = 0.6)+
  geom_smooth(se = F)+
  scale_fill_manual(values = col_treatment)+
  scale_color_manual(values = col_treatment)+
  facet_wrap(~index, scale = "free")+
  labs(x = "Alpha index values")+
  theme_fmt
  
# ── 1i. Rarefaction curves ────────────────────────────────────────────────────

cat("\n── 1i: Rarefaction curves ──────────────────────────────────────────────\n")

# 計算 rarefaction（vegan::rarecurve）
otu_for_rare <- t(as.matrix(otu_table(ps_filt)))   # rows = samples

# 關鍵：用 matrix() 重新建構，完全脫離 S4
otu_for_rare <- matrix(
  as.integer(as.vector(otu_for_rare)),
  nrow     = nrow(otu_for_rare),
  ncol     = ncol(otu_for_rare),
  dimnames = list(rownames(otu_for_rare), colnames(otu_for_rare))
)

# step = 200 reads，最大到每個樣本的實際 depth
rare_list <- vegan::rarecurve(otu_for_rare, step = 200,
                              tidy = TRUE,          # 回傳 tidy dataframe
                              label = FALSE)

# tidy = TRUE 回傳 Sample / Sites / Individuals
rare_df <- rare_list %>%
  as_tibble() %>%
  rename(SampleID = Site, richness = Species, depth = Sample) %>%
  # rarecurve tidy 的欄位：Sample(=SampleID), Sites(=richness), Individuals(=depth)
  left_join(meta_raw %>% select(SampleID, Treatment, Timepoint, PatientID),
            by = "SampleID")

# rarecurve tidy 欄位名稱因版本而異，修正一下
if (!"SampleID" %in% names(rare_df)) {
  rare_df <- rare_list %>%
    as_tibble() %>%
    left_join(meta_raw %>% select(SampleID, Treatment, Timepoint, PatientID),
              by = c("Site" = "SampleID")) %>%
    rename(depth = Subsample, richness = Species, SampleID = Site)
}
slice_rare_df <- slice_max(rare_df, by = SampleID , order_by = depth)
p_rarefaction <- ggplot(rare_df,
                        aes(x = depth, y = richness,
                            colour = Treatment,
                            group  = SampleID,
                            linetype = Timepoint)) +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  geom_vline(xintercept = rarefaction_depth,
             linetype = "dashed", colour = "grey30", linewidth = 0.8) +
  geom_text_repel(data = slice_rare_df, aes(label = SampleID))+
  facet_wrap(~Treatment)+
  annotate("text",
           x     = rarefaction_depth * 1.02,
           y     = max(rare_df$richness, na.rm = TRUE) * 0.95,
           label = paste("Rarefaction\ndepth =", formatC(rarefaction_depth,
                                                         format = "d", big.mark = ",")),
           hjust = 0, size = 3, colour = "grey30") +
  scale_colour_manual(values = col_treatment) +
  scale_linetype_manual(values = c(Baseline = "solid",
                                   `2M post-FMT` = "dashed",
                                   `6M post-FMT` = "dotted")) +
  labs(title    = "Rarefaction curves",
       subtitle = "Vertical line = rarefaction depth used in analysis",
       x        = "Sequencing depth (reads)",
       y        = "Observed ASVs",
       colour   = "Treatment",
       linetype = "Time point") +
  theme_fmt
