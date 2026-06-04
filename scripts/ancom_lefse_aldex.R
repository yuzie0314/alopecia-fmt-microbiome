# =============================================================================
#  ancom_lefse_aldex.R
#  Differential abundance — ANCOM-BC2 + LEfSe (lefser) + ALDEx2
#
#  版本確認：
#    ANCOMBC  2.14.0  (R 4.6.0, Bioc 3.23)
#    lefser   1.22.0
#    ALDEx2   1.44.0
#
#  分析設計：
#    Part 1 — A and B not C（FMT-specific @ 6M）
#    Part 2 — Pairwise timepoint（BL→2M / 2M→6M / BL→6M）
#    Data   — Taxonomy (genus) + MetaCyc + EC
#
#  前置條件：
#    ps_filt, meta_raw, col_treatment, tp_labels, theme_fmt 已定義
#    metacyc_raw, ec_raw（rows=features, cols=SampleID, raw RPK）
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

library(ANCOMBC)
library(lefser)
library(ALDEx2)
library(SummarizedExperiment)
library(phyloseq)
library(tidyverse)
library(ggrepel)
library(patchwork)

NO_BL_PATIENTS <- "GM03"

TRANSITIONS <- list(
  BL_2M = c("Baseline",    "2M post-FMT"),
  M2_6M = c("2M post-FMT", "6M post-FMT"),
  BL_6M = c("Baseline",    "6M post-FMT")
)

# 輸出目錄
OUT_ROOT <- file.path(getwd(), "ancom_lefse_aldex_results")
for (dtype in c("taxonomy", "metacyc", "ec", "combined")) {
  if (dtype == "combined") {
    dir.create(file.path(OUT_ROOT, "combined"), recursive = TRUE,
               showWarnings = FALSE)
  } else {
    for (tp in c("A_not_C", names(TRANSITIONS)))
      dir.create(file.path(OUT_ROOT, dtype, tp),
                 recursive = TRUE, showWarnings = FALSE)
  }
}


# =============================================================================
#  輔助 functions：格式轉換
# =============================================================================

# ── phyloseq → SummarizedExperiment（for lefser）────────────────────────────

ps_to_se <- function(ps_obj) {
  if (is.null(ps_obj) || nsamples(ps_obj) == 0 || ntaxa(ps_obj) == 0) {
    message("ps_to_se: empty phyloseq, returning NULL")
    return(NULL)
  }
  
  count_mat <- as.matrix(otu_table(ps_obj))
  if (!taxa_are_rows(ps_obj)) count_mat <- t(count_mat)

  col_data <- as.data.frame(sample_data(ps_obj))

  SummarizedExperiment(
    assays   = list(counts = count_mat),
    colData  = col_data
  )
}

# ── matrix → SummarizedExperiment（for pathway lefser）─────────────────────

mat_to_se <- function(mat_features_samples, meta_sub) {
  # mat_features_samples: rows = features, cols = SampleID
  common   <- intersect(colnames(mat_features_samples), meta_sub$SampleID)
  mat_sub  <- mat_features_samples[, common, drop = FALSE]
  col_data <- meta_sub %>%
    filter(SampleID %in% common) %>%
    remove_rownames() %>% 
    column_to_rownames("SampleID")

  SummarizedExperiment(
    assays  = list(counts = round(mat_sub)),
    colData = col_data
  )
}

# ── 統一結果格式 ──────────────────────────────────────────────────────────────

tidy_ancombc2 <- function(res, comparison, group, data_type,
                          feature_col = "taxon") {
  if (is.null(res)) return(NULL)
  
  res_df <- res$res
  cat("  ANCOMBC2 columns:", paste(names(res_df), collapse = ", "), "\n")
  
  # 排除 Intercept，找真正的 treatment/timepoint effect 欄位
  lfc_col  <- grep("^lfc_",  names(res_df), value = TRUE) %>%
    .[!grepl("Intercept|PatientID", .)] %>% .[1]
  se_col   <- grep("^se_",   names(res_df), value = TRUE) %>%
    .[!grepl("Intercept|PatientID", .)] %>% .[1]
  pval_col <- grep("^p_",    names(res_df), value = TRUE) %>%
    .[!grepl("Intercept|PatientID", .)] %>% .[1]
  padj_col <- grep("^q_",    names(res_df), value = TRUE) %>%
    .[!grepl("Intercept|PatientID", .)] %>% .[1]
  diff_col <- grep("^diff_", names(res_df), value = TRUE) %>%
    .[!grepl("Intercept|PatientID", .)] %>% .[1]
  
  cat(sprintf("  Using: lfc=%s  pval=%s  diff=%s\n",
              lfc_col, pval_col, diff_col))
  
  # 如果還是找不到（例如只有兩組但欄位名稱不含 group name）
  # fallback：用第一個非 Intercept 的 lfc 欄位
  if (is.na(lfc_col)) {
    lfc_col  <- grep("^lfc_",  names(res_df), value = TRUE) %>%
      .[!grepl("Intercept", .)] %>% .[1]
    se_col   <- grep("^se_",   names(res_df), value = TRUE) %>%
      .[!grepl("Intercept", .)] %>% .[1]
    pval_col <- grep("^p_",    names(res_df), value = TRUE) %>%
      .[!grepl("Intercept", .)] %>% .[1]
    padj_col <- grep("^q_",    names(res_df), value = TRUE) %>%
      .[!grepl("Intercept", .)] %>% .[1]
    diff_col <- grep("^diff_", names(res_df), value = TRUE) %>%
      .[!grepl("Intercept", .)] %>% .[1]
    cat(sprintf("  Fallback: lfc=%s  pval=%s  diff=%s\n",
                lfc_col, pval_col, diff_col))
  }
  
  if (is.na(lfc_col) || is.na(pval_col) || is.na(diff_col)) {
    message("  ANCOMBC2: cannot find effect columns, returning NULL")
    cat("  Available columns:\n")
    print(names(res_df))
    return(NULL)
  }
  
  res_df %>%
    rename(feature_id = !!feature_col) %>%
    dplyr::mutate(
      method     = "ANCOM-BC2",
      data_type  = data_type,
      comparison = comparison,
      group      = group,
      lfc        = .data[[lfc_col]],
      se         = if (!is.na(se_col)) .data[[se_col]] else NA_real_,
      pval       = .data[[pval_col]],
      padj       = if (!is.na(padj_col)) .data[[padj_col]] else NA_real_,
      sig_ancom  = .data[[diff_col]],
      direction  = case_when(
        sig_ancom & lfc > 0 ~ "Up",
        sig_ancom & lfc < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) 
  
}

tidy_lefser <- function(res_df, comparison, group, data_type) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  res_df %>%
    as_tibble() %>%
    rename(feature_id = Names, lda_score = scores) %>%
    dplyr::mutate(
      method     = "LEfSe",
      data_type  = data_type,
      comparison = comparison,
      group      = group,
      direction  = if_else(lda_score > 0, "Up", "Down")
    ) %>%
    dplyr::select(method, data_type, comparison, group,
           feature_id, lda_score, direction)
}

tidy_aldex2 <- function(res_df, comparison, group, data_type) {
  if (is.null(res_df)) return(NULL)
  
  # aldex.ttest 欄位：we.ep, we.eBH（Welch）, wi.ep, wi.eBH（Wilcoxon）
  # aldex.effect 欄位：effect, overlap
  # 用 wi.eBH（Wilcoxon BH-corrected）作為主要顯著性
  
  res_df <- tibble::rownames_to_column(as.data.frame(res_df), "feature_id")
  
  res_df$method     <- "ALDEx2"
  res_df$data_type  <- data_type
  res_df$comparison <- comparison
  res_df$group      <- group
  
  # effect 欄位確認存在
  effect_col <- if ("effect" %in% names(res_df)) "effect" else NULL
  pval_col   <- if ("wi.eBH" %in% names(res_df)) "wi.eBH" else
    if ("we.eBH" %in% names(res_df)) "we.eBH" else NULL
  
  if (is.null(pval_col)) {
    message("ALDEx2: cannot find p-value column")
    return(NULL)
  }
  
  res_df$direction <- ifelse(
    res_df[[pval_col]] < 0.05 & !is.null(effect_col) & res_df[[effect_col]] > 0, "Up",
    ifelse(
      res_df[[pval_col]] < 0.05 & !is.null(effect_col) & res_df[[effect_col]] < 0, "Down",
      "NS"
    )
  )
  
  keep_cols <- c("method", "data_type", "comparison", "group",
                 "feature_id", effect_col, "wi.ep", "wi.eBH",
                 "we.ep", "we.eBH", "direction")
  keep_cols <- keep_cols[keep_cols %in% names(res_df)]
  
  res_df[, keep_cols, drop = FALSE] %>%
    dplyr::arrange(.data[[pval_col]])
}

save_result <- function(df, dtype, subfolder, label) {
  if (is.null(df) || nrow(df) == 0) return(invisible())
  fname <- sprintf("%s_%s.csv", label, dtype)
  fpath <- file.path(OUT_ROOT, dtype, subfolder, fname)
  write.csv(df, fpath, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", fpath))
}


# =============================================================================
#  ANCOM-BC2 runner
# =============================================================================

run_ancombc2 <- function(ps_obj, fix_formula, group_var,
                          reference_levels = NULL,
                          struc_zero = TRUE) {
  # round to integer
  otu_table(ps_obj) <- otu_table(
    round(otu_table(ps_obj)), taxa_are_rows = taxa_are_rows(ps_obj)
  )

  tryCatch(
    ancombc2(
      data             = ps_obj,
      fix_formula      = fix_formula,
      rand_formula     = NULL,       # ANCOMBC 2.x: no rand_formula in ancombc2
      p_adj_method     = "BH",
      alpha            = 0.05,
      prv_cut          = 0.10,
      lib_cut          = 0,
      s0_perc          = 0.05,
      group            = group_var,
      struc_zero       = struc_zero,
      neg_lb           = TRUE,
      global           = FALSE,
      pairwise         = FALSE,
      dunnet           = FALSE,
      trend            = FALSE,
      verbose          = FALSE
    ),
    error = function(e) {
      message("ANCOM-BC2 error: ", e$message)
      NULL
    }
  )
}


# =============================================================================
#  LEfSe runner（lefser 1.22）
# =============================================================================

run_lefser <- function(se_obj, class_col, ref_class = NULL,
                       lda_threshold = 2.0) {
  tryCatch({
    col_data <- as.data.frame(colData(se_obj))
    # droplevels 確保 factor 只剩當前樣本的 levels（避免 lefser "must be dichotomous" 錯誤）
    col_data[[class_col]] <- droplevels(factor(col_data[[class_col]]))
    if (!is.null(ref_class) && ref_class %in% levels(col_data[[class_col]]))
      col_data[[class_col]] <- relevel(col_data[[class_col]], ref = ref_class)
    colData(se_obj) <- DataFrame(col_data)

    # lefser 1.22：先轉 relative abundance，再用 relab 參數
    se_relab <- relativeAb(se_obj)
    
    lefser(
      relab             = se_relab,      # ← 1.22 用 relab 不是 expr
      kruskal.threshold = 0.05,
      wilcox.threshold  = 0.05,
      lda.threshold     = lda_threshold,
      groupCol          = class_col,
      blockCol          = NULL
    )
  }, error = function(e) {
    message("LEfSe error: ", e$message)
    NULL
  })
}


# =============================================================================
#  ALDEx2 runner
# =============================================================================

run_aldex2 <- function(count_mat, conditions, paired = FALSE,
                       mc_samples = 128) {
  tryCatch({
    clr <- aldex.clr(
      reads      = round(count_mat),
      conds      = conditions,
      mc.samples = mc_samples, # 64
      denom      = "all",
      verbose    = FALSE
    )
    
    # ttest 回傳 p 值
    tt <- aldex.ttest(clr, paired.test = paired, verbose = FALSE)
    
    # effect size 要另外算
    eff <- aldex.effect(clr, verbose = FALSE)
    
    # 合併
    cbind(tt, eff)
    
  }, error = function(e) {
    message("ALDEx2 error: ", e$message)
    NULL
  })
}


# =============================================================================
#  主分析 function：一個 comparison 跑三個方法
# =============================================================================

run_all_methods <- function(ps_obj,              # phyloseq（taxonomy）or NULL
                             mat_features,         # rows=features cols=samples
                             meta_sub,             # metadata dataframe
                             group_col,            # "Treatment" or "Timepoint"
                             ref_level,            # reference level
                             alt_level,            # comparison level
                             comparison,           # label e.g. "BL_6M"
                             group,                # "FMT" or "placebo"
                             data_type,            # "taxonomy"/"metacyc"/"ec"
                             paired = FALSE,
                             pat_col = "PatientID") {

  cat(sprintf("\n  [%s | %s | %s] ref=%s alt=%s paired=%s\n",
              data_type, group, comparison, ref_level, alt_level, paired))

  # ── 準備共用物件 ──────────────────────────────────────────────────────────

  # phyloseq for taxonomy; build from mat for pathway
  if (is.null(ps_obj) && !is.null(mat_features)) {
    common  <- intersect(colnames(mat_features), meta_sub$SampleID)
    
    # 加這個檢查
    if (length(common) == 0) {
      message("  No common samples between mat_features and meta_sub")
      return(list(ancom = NULL, lefse = NULL, aldex = NULL))
    }
    
    mat_sub <- round(mat_features[, common, drop = FALSE])
    mode(mat_sub) <- "integer"
    mat_sub <- mat_sub[rowSums(mat_sub) > 0, , drop = FALSE]

    meta_ps <- meta_sub %>%
      filter(SampleID %in% common) %>%
      remove_rownames() %>%
      column_to_rownames("SampleID")

    ps_use <- phyloseq(
      otu_table(mat_sub, taxa_are_rows = TRUE),
      sample_data(meta_ps)
    )
    ps_use <- filter_ps_by_prev(ps_obj = ps_use)
  } else {
    # taxonomy: subset phyloseq
    common  <- intersect(sample_names(ps_obj), meta_sub$SampleID)
    ps_use  <- prune_samples(common, ps_obj)
    mat_sub <- as.matrix(otu_table(ps_use))
    if (!taxa_are_rows(ps_use)) mat_sub <- t(mat_sub)
    meta_ps <- meta_sub %>%
      filter(SampleID %in% common) %>%
      remove_rownames() %>%
      column_to_rownames("SampleID")
    sample_data(ps_use) <- sample_data(meta_ps)
  }

  sample_data(ps_use)[[group_col]] <- relevel(
    factor(sample_data(ps_use)[[group_col]]), ref = ref_level
  )
  if (paired) {
    sample_data(ps_use)[[pat_col]] <- factor(
      sample_data(ps_use)[[pat_col]]
    )
  }

  # conditions vector for ALDEx2 (matched to sample order)
  meta_ordered <- meta_sub %>%
    filter(SampleID %in% colnames(mat_sub)) %>%
    arrange(match(SampleID, colnames(mat_sub)))

  conditions <- as.character(meta_ordered[[group_col]])
  # ALDEx2 needs exactly 2 levels
  conditions <- factor(conditions, levels = c(ref_level, alt_level))

  # ── ANCOM-BC2 ────────────────────────────────────────────────────────────

  fix_formula <- if (paired)
    paste0(pat_col, " + ", group_col)
  else
    group_col

  ancom_res  <- run_ancombc2(ps_use, fix_formula = fix_formula,
                               group_var = group_col)
  ancom_tidy <- tidy_ancombc2(ancom_res,
                               comparison = comparison,
                               group = group,
                               data_type = data_type)

  # add taxonomy for taxonomy data type
  if (data_type == "taxonomy" && !is.null(ancom_tidy)) {
    tax_df <- as.data.frame(as(tax_table(ps_use), "matrix")) %>%
      tibble::rownames_to_column("feature_id") %>%
      dplyr::select(feature_id, Genus, Family, Phylum)
    ancom_tidy <- left_join(ancom_tidy, tax_df, by = "feature_id") %>%
      dplyr::mutate(feature_name = coalesce(Genus, feature_id))
  } else if (!is.null(ancom_tidy)) {
    ancom_tidy <- ancom_tidy %>% dplyr::mutate(feature_name = feature_id)
  }

  cat(sprintf("    ANCOM-BC2 sig: %d\n",
              sum(ancom_tidy$direction != "NS", na.rm = TRUE)))

  # ── LEfSe ────────────────────────────────────────────────────────────────
  # 在 run_all_methods 裡，se_obj <- ps_to_se(ps_use) 之前加
  # drop unused factor levels，確保 group_col 只有兩個 level
  sample_data(ps_use)[[group_col]] <- droplevels(
    factor(sample_data(ps_use)[[group_col]])
  )
  
  cat(sprintf("  %s levels: %s\n",
              group_col,
              paste(levels(sample_data(ps_use)[[group_col]]), collapse = ", ")))
  
  se_obj <- ps_to_se(ps_use)
  lefse_tidy <- if (!is.null(se_obj)) {
    lefse_raw  <- run_lefser(se_obj, class_col = group_col,
                             ref_class = ref_level)
    tidy_lefser(lefse_raw, comparison = comparison,
                group = group, data_type = data_type)
  } else {
    message("  LEfSe skipped: empty SummarizedExperiment")
    NULL
  }

  if (data_type == "taxonomy" && !is.null(lefse_tidy)) {
    tax_df <- as.data.frame(as(tax_table(ps_use), "matrix")) %>%
      tibble::rownames_to_column("feature_id") %>%
      dplyr::select(feature_id, Genus)
    lefse_tidy <- left_join(lefse_tidy, tax_df, by = "feature_id") %>%
      dplyr::mutate(feature_name = coalesce(Genus, feature_id))
  } else if (!is.null(lefse_tidy)) {
    lefse_tidy <- lefse_tidy %>% dplyr::mutate(feature_name = feature_id)
  }

  cat(sprintf("    LEfSe sig (LDA>2): %d\n",
              if (!is.null(lefse_tidy)) nrow(lefse_tidy) else 0))

  # ── ALDEx2 ───────────────────────────────────────────────────────────────

  aldex_raw  <- run_aldex2(mat_sub, as.character(conditions),
                             paired = paired)
  aldex_tidy <- tidy_aldex2(aldex_raw,
                              comparison = comparison,
                              group = group,
                              data_type = data_type)

  if (data_type == "taxonomy" && !is.null(aldex_tidy)) {
    tax_df <- as.data.frame(as(tax_table(ps_use), "matrix")) %>%
      tibble::rownames_to_column("feature_id") %>%
      dplyr::select(feature_id, Genus)
    aldex_tidy <- left_join(aldex_tidy, tax_df, by = "feature_id") %>%
      dplyr::mutate(feature_name = coalesce(Genus, feature_id))
  } else if (!is.null(aldex_tidy)) {
    aldex_tidy <- aldex_tidy %>% dplyr::mutate(feature_name = feature_id)
  }

  cat(sprintf("    ALDEx2 sig (wi.eBH<0.05): %d\n",
              sum(aldex_tidy$direction != "NS", na.rm = TRUE)))

  list(ancom = ancom_tidy, lefse = lefse_tidy, aldex = aldex_tidy)
}


# =============================================================================
#  PART 1 — A and B not C（FMT-specific @ 6M）
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("  PART 1: A and B not C\n")
cat(strrep("=", 60), "\n")

# Genus-level phyloseq（genus agglomerate）
ps_genus <- tax_glom(ps_filt, taxrank = "Genus", NArm = FALSE)
ps_genus <- filter_ps_by_prev(ps_obj = ps_genus)
part1_results <- list()

for (dtype in c("taxonomy", "metacyc", "ec")) {

  mat_feat <- switch(dtype,
                      taxonomy = NULL,
                      metacyc  = metacyc_raw,
                      ec       = ec_raw)

  ps_use_tax <- if (dtype == "taxonomy") ps_genus else NULL

  # ── Comparison A：FMT vs placebo @ 6M ────────────────────────────────────

  meta_A <- meta_raw %>%
    filter(Timepoint == "6M post-FMT") %>%
    dplyr::mutate(Treatment = factor(Treatment, levels = c("placebo", "FMT")))

  res_A <- run_all_methods(
    ps_obj      = if (dtype == "taxonomy")
                    subset_samples(ps_genus, Timepoint == "6M post-FMT")
                  else NULL,
    mat_features = mat_feat,
    meta_sub    = meta_A,
    group_col   = "Treatment",
    ref_level   = "placebo",
    alt_level   = "FMT",
    comparison  = "A_FMT_vs_placebo_6M",
    group       = "all",
    data_type   = dtype,
    paired      = FALSE
  )

  # ── Comparison B：FMT Baseline→6M（paired）────────────────────────────────

  meta_B <- meta_raw %>%
    filter(Treatment == "FMT",
           Timepoint %in% c("Baseline", "6M post-FMT"),
           !PatientID %in% NO_BL_PATIENTS) %>%
    dplyr::mutate(Timepoint = factor(Timepoint,
                               levels = c("Baseline", "6M post-FMT")))

  res_B <- run_all_methods(
    ps_obj      = if (dtype == "taxonomy")
                    subset_samples(ps_genus,
                                   Treatment == "FMT" &
                                   Timepoint %in% c("Baseline", "6M post-FMT") &
                                   !PatientID %in% NO_BL_PATIENTS)
                  else NULL,
    mat_features = mat_feat,
    meta_sub    = meta_B,
    group_col   = "Timepoint",
    ref_level   = "Baseline",
    alt_level   = "6M post-FMT",
    comparison  = "B_FMT_BL_vs_6M",
    group       = "FMT",
    data_type   = dtype,
    paired      = TRUE
  )

  # ── Comparison C：Placebo Baseline→6M ─────────────────────────────────────

  meta_C <- meta_raw %>%
    filter(Treatment == "placebo",
           Timepoint %in% c("Baseline", "6M post-FMT")) %>%
    dplyr::mutate(Timepoint = factor(Timepoint,
                               levels = c("Baseline", "6M post-FMT")))

  res_C <- run_all_methods(
    ps_obj      = if (dtype == "taxonomy")
                    subset_samples(ps_genus,
                                   Treatment == "placebo" &
                                   Timepoint %in% c("Baseline", "6M post-FMT"))
                  else NULL,
    mat_features = mat_feat,
    meta_sub    = meta_C,
    group_col   = "Timepoint",
    ref_level   = "Baseline",
    alt_level   = "6M post-FMT",
    comparison  = "C_placebo_BL_vs_6M",
    group       = "placebo",
    data_type   = dtype,
    paired      = TRUE
  )

  # ── A and B not C per method ──────────────────────────────────────────────

  get_sig_features <- function(res_list, comp_letter, method_key,
                               direction_filter = "Up") {
    df <- res_list[[method_key]]
    
    # NULL 或空 dataframe 都回傳 character(0)
    if (is.null(df) || nrow(df) == 0) return(character(0))
    
    feat <- df %>%
      dplyr::filter(direction == direction_filter) %>%
      dplyr::pull(feature_name) %>%
      na.omit() %>%
      unique()
    
    if (length(feat) == 0) character(0) else as.character(feat)
  }

  for (mkey in c("ancom", "lefse", "aldex")) {
    up_A <- get_sig_features(res_A, "A", mkey, "Up")
    up_B <- get_sig_features(res_B, "B", mkey, "Up")
    any_C <- c(
      get_sig_features(res_C, "C", mkey, "Up"),
      get_sig_features(res_C, "C", mkey, "Down")
    )

    fmt_specific <- intersect(up_A, up_B) %>% setdiff(any_C)

    cat(sprintf("  [%s | %s] FMT-specific (A∩B not C): %d features\n",
                dtype, mkey, length(fmt_specific)))

    if (length(fmt_specific) > 0) {
      out_df <- data.frame(
        feature_name = fmt_specific,
        method       = mkey,
        data_type    = dtype
      )
      save_result(out_df, dtype, "A_not_C",
                  sprintf("%s_FMT_specific_%s", mkey, dtype))
    } else {
      cat(sprintf("  [%s | %s] No FMT-specific features found\n", dtype, mkey))
    }
  }

  # save individual comparison results
  for (comp_res in list(A = res_A, B = res_B, C = res_C)) {
    for (mkey in c("ancom", "lefse", "aldex")) {
      df <- comp_res[[mkey]]
      if (!is.null(df))
        save_result(df, dtype, "A_not_C",
                    sprintf("%s_%s_%s", mkey, df$comparison[1], dtype))
    }
  }

  part1_results[[dtype]] <- list(A = res_A, B = res_B, C = res_C)
}


# =============================================================================
#  PART 2 — Pairwise timepoint（BL→2M / 2M→6M / BL→6M）
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("  PART 2: Pairwise timepoint\n")
cat(strrep("=", 60), "\n")

all_pw_results <- list()

for (tp_key in names(TRANSITIONS)) {
  tp_pair <- TRANSITIONS[[tp_key]]
  tp_from <- tp_pair[1]
  tp_to   <- tp_pair[2]

  excl <- if (grepl("Baseline", tp_from)) NO_BL_PATIENTS else character(0)

  for (grp in c("FMT", "placebo")) {
    for (dtype in c("taxonomy", "metacyc", "ec")) {

      meta_sub <- meta_raw %>%
        filter(Treatment == grp,
               Timepoint %in% tp_pair,
               !PatientID %in% excl) %>%
        dplyr::mutate(Timepoint = factor(Timepoint, levels = tp_pair))

      if (nrow(meta_sub) < 4) next

      mat_feat  <- switch(dtype,
                           taxonomy = NULL,
                           metacyc  = metacyc_raw,
                           ec       = ec_raw)

      ps_tax_sub <- if (dtype == "taxonomy")
        subset_samples(ps_genus,
                       Treatment == grp &
                       Timepoint %in% tp_pair &
                       !PatientID %in% excl)
      else NULL

      res_pw <- run_all_methods(
        ps_obj       = ps_tax_sub,
        mat_features = mat_feat,
        meta_sub     = meta_sub,
        group_col    = "Timepoint",
        ref_level    = tp_from,
        alt_level    = tp_to,
        comparison   = tp_key,
        group        = grp,
        data_type    = dtype,
        paired       = TRUE
      )

      # save individual results
      key <- paste(tp_key, grp, dtype, sep = "_")
      all_pw_results[[key]] <- res_pw

      for (mkey in c("ancom", "lefse", "aldex")) {
        df <- res_pw[[mkey]]
        if (!is.null(df))
          save_result(df, dtype, tp_key,
                      sprintf("%s_%s_%s", mkey, grp, tp_key))
      }
    }
  }
}


# =============================================================================
#  PART 3 — 整合 long format + combined CSV
# =============================================================================

cat("\n── Combining all results ────────────────────────────────────────────────\n")

collect_all <- function(results_list, methods = c("ancom", "lefse", "aldex")) {
  purrr::map_dfr(results_list, function(res_grp) {
    purrr::map_dfr(methods, function(mkey) {
      res_grp[[mkey]]
    })
  })
}

# Part 1 combined
part1_combined <- purrr::map_dfr(part1_results, function(x) {
  purrr::map_dfr(list(A = x$A, B = x$B, C = x$C), function(comp) {
    purrr::map_dfr(c("ancom", "lefse", "aldex"), ~ comp[[.x]])
  })
})

# Part 2 combined
part2_combined <- collect_all(all_pw_results)

write.csv(part1_combined,
          file.path(OUT_ROOT, "combined", "part1_A_B_C_all_methods.csv"),
          row.names = FALSE)
write.csv(part2_combined,
          file.path(OUT_ROOT, "combined", "part2_pairwise_all_methods.csv"),
          row.names = FALSE)

cat(sprintf("Part 1 combined: %d rows\n", nrow(part1_combined)))
cat(sprintf("Part 2 combined: %d rows\n", nrow(part2_combined)))


# =============================================================================
#  PART 4 — 視覺化
# =============================================================================

# ── Helper: method concordance dot plot ──────────────────────────────────────

plot_method_concordance <- function(df_list, dtype, comparison_label) {
  # df_list: list(ancom=df, lefse=df, aldex=df)

  ancom_up <- df_list$ancom %>%
    filter(direction == "Up") %>% pull(feature_name) %>% na.omit()
  lefse_up <- df_list$lefse %>%
    filter(direction == "Up") %>% pull(feature_name) %>% na.omit()
  aldex_up <- df_list$aldex %>%
    filter(direction == "Up") %>% pull(feature_name) %>% na.omit()

  all_features <- unique(c(ancom_up, lefse_up, aldex_up))
  if (length(all_features) == 0) return(NULL)

  tbl <- tibble(feature_name = all_features) %>%
    dplyr::mutate(
      ANCOM_BC2 = feature_name %in% ancom_up,
      LEfSe     = feature_name %in% lefse_up,
      ALDEx2    = feature_name %in% aldex_up,
      n_methods = ANCOM_BC2 + LEfSe + ALDEx2
    ) %>%
    arrange(desc(n_methods))

  ggplot(tbl, aes(x = reorder(feature_name, n_methods),
                   y = n_methods, fill = factor(n_methods))) +
    geom_col(width = 0.7) +
    geom_point(data = tbl %>% filter(ANCOM_BC2),
               aes(x = feature_name, y = 0.2),
               colour = "#D6604D", shape = 16, size = 2.5,
               inherit.aes = FALSE) +
    geom_point(data = tbl %>% filter(LEfSe),
               aes(x = feature_name, y = 0.5),
               colour = "#1D9E75", shape = 17, size = 2.5,
               inherit.aes = FALSE) +
    geom_point(data = tbl %>% filter(ALDEx2),
               aes(x = feature_name, y = 0.8),
               colour = "#378ADD", shape = 15, size = 2.5,
               inherit.aes = FALSE) +
    coord_flip() +
    scale_fill_manual(
      values = c("1" = "grey75", "2" = "#BA7517", "3" = "#D79B00"),
      name   = "# methods"
    ) +
    scale_y_continuous(breaks = 1:3) +
    labs(
      title    = sprintf("Method concordance: %s | %s (Up in FMT/post)",
                         dtype, comparison_label),
      subtitle = "● ANCOM-BC2  ▲ LEfSe  ■ ALDEx2",
      x = NULL, y = "# concordant methods"
    ) +
    theme_fmt +
    theme(axis.text.y = element_text(
      face = if (dtype == "taxonomy") "italic" else "plain", size = 8))
}

# ── LDA bar（LEfSe）──────────────────────────────────────────────────────────

plot_lda_bar <- function(df, dtype, label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  top_df <- df %>%
    arrange(desc(abs(lda_score))) %>%
    head(20) %>%
    dplyr::mutate(feature_name = factor(feature_name,
                                  levels = feature_name[order(lda_score)]))

  ggplot(top_df, aes(x = feature_name, y = lda_score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(Up = "#D6604D", Down = "#2166AC")) +
    labs(title    = sprintf("LEfSe LDA score: %s | %s", dtype, label),
         subtitle = "LDA threshold = 2.0",
         x = NULL, y = "LDA score", fill = NULL) +
    theme_fmt +
    theme(axis.text.y = element_text(
      face = if (dtype == "taxonomy") "italic" else "plain", size = 8))
}

# ── ALDEx2 effect plot ────────────────────────────────────────────────────────

plot_aldex_effect <- function(df, dtype, label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = effect, y = -log10(wi.eBH), colour = direction)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_text_repel(
      data = filter(df, wi.eBH < 0.05),
      aes(label = feature_name), size = 2.5,
      fontface = if (dtype == "taxonomy") "italic" else "plain",
      max.overlaps = 12, segment.colour = "grey60"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", colour = "grey50") +
    scale_colour_manual(
      values = c(Up = "#D6604D", Down = "#2166AC", NS = "grey75")
    ) +
    labs(title    = sprintf("ALDEx2 effect: %s | %s", dtype, label),
         subtitle = "wi.eBH < 0.05 labelled",
         x = "Effect size", y = "-log₁₀(BH-adj Wilcoxon p)",
         colour = NULL) +
    theme_fmt
}

# ── Generate and save plots ───────────────────────────────────────────────────

save_pw_fig <- function(plot, filename, w = 10, h = 8) {
  if (is.null(plot)) return(invisible())
  fpath <- file.path(OUT_ROOT, "combined", filename)
  ggsave(fpath, plot = plot, width = w, height = h, dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", fpath))
}

# Part 1 plots
for (dtype in c("taxonomy", "metacyc", "ec")) {
  res_grp <- part1_results[[dtype]]

  p_conc <- plot_method_concordance(res_grp$A, dtype, "FMT vs placebo @ 6M")
  save_pw_fig(p_conc, sprintf("p1_concordance_%s_A.pdf", dtype))

  p_lda  <- plot_lda_bar(res_grp$A$lefse, dtype, "FMT vs placebo @ 6M")
  save_pw_fig(p_lda, sprintf("p1_lefse_lda_%s_A.pdf", dtype))

  p_eff  <- plot_aldex_effect(res_grp$A$aldex, dtype, "FMT vs placebo @ 6M")
  save_pw_fig(p_eff, sprintf("p1_aldex_effect_%s_A.pdf", dtype))
}

# Part 2 plots（LDA heatmap across transitions）
for (dtype in c("taxonomy", "metacyc", "ec")) {
  for (grp in c("FMT", "placebo")) {

    # Collect LEfSe results across transitions
    lefse_timeline <- purrr::map_dfr(names(TRANSITIONS), function(tp) {
      key <- paste(tp, grp, dtype, sep = "_")
      df  <- all_pw_results[[key]]$lefse
      if (!is.null(df)) df else NULL
    })

    if (!is.null(lefse_timeline) && nrow(lefse_timeline) > 0) {
      top_feat <- lefse_timeline %>%
        group_by(feature_name) %>%
        summarise(max_lda = max(abs(lda_score)), .groups = "drop") %>%
        slice_max(max_lda, n = 30) %>% pull(feature_name)

      hm_df <- lefse_timeline %>%
        filter(feature_name %in% top_feat) %>%
        dplyr::mutate(
          lda_signed = if_else(direction == "Up", lda_score, -lda_score),
          transition = factor(comparison, levels = names(TRANSITIONS))
        )

      p_lefse_hm <- ggplot(hm_df,
                            aes(x = transition, y = feature_name,
                                fill = lda_signed)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        scale_fill_gradient2(
          low = "#2166AC", mid = "white", high = "#D6604D",
          midpoint = 0, name = "Signed LDA"
        ) +
        scale_x_discrete(labels = c(BL_2M = "BL→2M",
                                     M2_6M = "2M→6M",
                                     BL_6M = "BL→6M")) +
        labs(
          title = sprintf("LEfSe LDA timeline: %s | %s", dtype, grp),
          x = "Transition", y = NULL
        ) +
        theme_fmt +
        theme(axis.text.y = element_text(
          face = if (dtype == "taxonomy") "italic" else "plain", size = 7))

      save_pw_fig(p_lefse_hm,
                  sprintf("p2_lefse_timeline_%s_%s.pdf", dtype, grp),
                  w = 8, h = 10)
    }
  }
}


# =============================================================================
#  Summary
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("✓ ANCOM-BC2 + LEfSe + ALDEx2 analysis complete\n")
cat(strrep("=", 60), "\n\n")

# FMT-specific summary across all data types
cat("FMT-specific features (A∩B not C) summary:\n")
for (dtype in c("taxonomy", "metacyc", "ec")) {
  for (mkey in c("ancom", "lefse", "aldex")) {
    fpath <- file.path(OUT_ROOT, dtype, "A_not_C",
                       sprintf("%s_FMT_specific_%s.csv", mkey, dtype))
    if (file.exists(fpath)) {
      df <- read.csv(fpath)
      cat(sprintf("  [%s | %s]: %d features\n", dtype, mkey, nrow(df)))
    }
  }
}

cat("\nOutput: ", OUT_ROOT, "\n")
cat("Combined CSVs:\n")
cat("  combined/part1_A_B_C_all_methods.csv\n")
cat("  combined/part2_pairwise_all_methods.csv\n")
