# =============================================================================
#  data_prep.R  (v2 — 對應真實 metadata.xlsx 結構)
#
#  metadata.xlsx 有 4 個 sheet：
#    ids  → analysis_id / lab_id / sequencing_id (sample mapping)
#    v1   → baseline 臨床資料 + SALT score
#    v2   → 2M visit SALT score
#    v3   → 3M visit SALT score + 抗生素記錄
#
#  輸出：
#    asv_table.csv   rows=ASVs, cols=SampleID (lab_id)
#    taxonomy.csv    ASV_ID + 7 rank columns
#    metadata.csv    SampleID / PatientID / Timepoint / Treatment /
#                    SALT_score / Age / Race / antibiotics
#
#  用法：
#    source("data_prep.R")
#    source("alopecia_fmt_16S_analysis.R")
# =============================================================================

library(tidyverse)
library(readxl)

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

# ── 0. File lists ─────────────────────────────────────────────────────────────

asvs     <- list.files("ASV", pattern = "ASV.tsv",         recursive = TRUE, full.names = TRUE)
pathways_ec <- list.files("PATHWAY", pattern = "ec.csv", recursive = TRUE, full.names = TRUE)
pathways_metacyc <- list.files("PATHWAY", pattern = "metacyc.csv", recursive = TRUE, full.names = TRUE)

TAX_COLS <- c("Var.1", "superkingdom", "phylum", "clade", "class",
              "order", "family", "genus", "species", "strain", "taxon")


# ── 1. 讀取 metadata.xlsx ─────────────────────────────────────────────────────

ids_sheet <- read_excel("metadata.xlsx", sheet = "ids") %>%
  select(subject_id, lab_id, analysis_id, status) %>%
  filter(status == "s")          # 只保留 status = "s" 的樣本
# GM03_V1 是 "f"，會被自動排除

v1 <- read_excel("metadata.xlsx", sheet = "v1") %>%
  transmute(
    PatientID   = `Subject ID`,
    SampleID    = lab_id,
    Sex         = Sex,
    Treatment   = Treatment,
    Age         = Age,
    Race        = Race,
    SALT_score  = `Total score`,
    Timepoint   = "baseline",
    antibiotics = NA_character_
  )

v2 <- read_excel("metadata.xlsx", sheet = "v2") %>%
  transmute(
    PatientID   = `Subject ID`,
    SampleID    = lab_id,
    SALT_score  = `Total score`,
    Timepoint   = "month2",
    antibiotics = NA_character_
  )

v3 <- read_excel("metadata.xlsx", sheet = "v3") %>%
  transmute(
    PatientID   = `Subject ID`,
    SampleID    = lab_id,
    SALT_score  = `Total score`,
    Timepoint   = "month6",
    antibiotics = `Did patient take any oral antibiotics medication after the last visit?`
  )

# v2 / v3 缺少 Treatment / Age / Race，從 v1 補入
patient_info <- v1 %>%
  select(PatientID, Sex, Treatment, Age, Race) %>%
  distinct()

metadata <- bind_rows(v1, v2, v3) %>%
  left_join(patient_info, by = "PatientID", suffix = c("", ".fill")) %>%
  mutate(
    Treatment = coalesce(Treatment, Treatment.fill),
    Sex       = coalesce(Sex,       Sex.fill),
    Age       = coalesce(Age,       Age.fill),
    Race      = coalesce(Race,      Race.fill)
  ) %>%
  select(-ends_with(".fill")) %>%
  # 統一 Treatment 大小寫：FMT / placebo
  mutate(
    Treatment = case_when(
      str_to_lower(Treatment) == "fmt"     ~ "FMT",
      str_to_lower(Treatment) == "placebo" ~ "placebo",
      TRUE ~ Treatment
    ),
    Timepoint = factor(Timepoint,
                       levels = c("baseline", "month2", "month6"))
  ) %>%
  # 加入 analysis_id（僅用於對應 ASV table，之後不輸出）
  left_join(
    ids_sheet %>% select(lab_id, analysis_id),
    by = c("SampleID" = "lab_id")
  ) %>%
  # 排除 status="f" 的樣本（ids_sheet 已過濾，left_join 後 analysis_id 為 NA）
  filter(!is.na(analysis_id))

cat("metadata 完成：", nrow(metadata), "筆樣本\n")
print(metadata %>% select(SampleID, PatientID, Timepoint, Treatment, SALT_score))


# ── 2. 合併 ASV tables ────────────────────────────────────────────────────────

combine_asv <- function(file_lst) {
  df_list <- lapply(file_lst, read_csv, show_col_types = FALSE)
  purrr::reduce(df_list, full_join, by = TAX_COLS)
}

tax_raw <- combine_asv(asvs)
names(tax_raw) <- names(tax_raw) %>% str_remove("\\.16S\\.exp\\.")

cat("tax_raw:", nrow(tax_raw), "ASVs ×", ncol(tax_raw), "cols\n")


# ── 3. Mapper：analysis_id → SampleID (lab_id) ────────────────────────────────

mapper <- metadata %>%
  select(analysis_id, SampleID) %>%
  deframe()                        # names = analysis_id, values = SampleID

analysis_cols_in_tax <- intersect(names(mapper), names(tax_raw))

missing_from_tax <- setdiff(names(mapper), names(tax_raw))
missing_from_map <- setdiff(names(tax_raw)[-(seq_along(TAX_COLS))], names(mapper))

if (length(missing_from_tax) > 0)
  warning("Mapper 有但 ASV table 找不到：\n  ",
          paste(missing_from_tax, collapse = "\n  "))
if (length(missing_from_map) > 0)
  message("ASV table 有但 mapper 沒有（略去）：\n  ",
          paste(missing_from_map, collapse = "\n  "))

cat("成功對應樣本數:", length(analysis_cols_in_tax), "\n")


# ── 4. 建立 taxonomy 和 ASV count table ──────────────────────────────────────

taxonomy <- tax_raw %>%
  select(all_of(TAX_COLS)) %>%
  mutate(ASV_ID = paste0("ASV", str_pad(row_number(), 4, pad = "0"))) %>%
  transmute(
    ASV_ID,
    Kingdom = superkingdom,
    Phylum  = phylum,
    Class   = class,
    Order   = order,
    Family  = family,
    Genus   = genus,
    Species = species
  )

asv_table <- tax_raw %>%
  select(all_of(analysis_cols_in_tax)) %>%
  rename_with(~ mapper[.x], .cols = everything()) %>%   # analysis_id → SampleID
  mutate(ASV_ID = taxonomy$ASV_ID) %>%
  select(ASV_ID, everything()) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0L))) %>%
  mutate(across(where(is.numeric), as.integer))

cat("asv_table:", nrow(asv_table), "ASVs ×", ncol(asv_table) - 1, "samples\n")


# ── 5. Pathway table（optional）──────────────────────────────────────────────

process_and_save <- function(file_lst, type_name, mapper) {
  if (length(file_lst) == 0) return(NULL)
  
  # 內部函數：讀取單一無 Header 的 CSV
  read_score <- function(fp) {
    # 提取樣本 ID (從資料夾名稱)
    aid <- basename(dirname(fp)) 
    
    # 讀取檔案，手動給予欄位名稱
    read_csv(fp, col_names = c("pathway", "abund"), show_col_types = FALSE) %>%
      transmute(pathway = pathway, !!aid := abund)
  }
  
  # 1. 批次讀取並合併
  df_list <- lapply(file_lst, read_score)
  pat_raw <- purrr::reduce(df_list, full_join, by = "pathway")
  
  # 2. 根據 mapper 篩選與重命名欄位
  # 找出 pat_raw 中存在於 mapper 內的樣本 ID
  pathway_cols <- intersect(names(mapper), names(pat_raw))
  
  pat <- pat_raw %>%
    select(pathway, all_of(pathway_cols)) %>%
    rename_with(~ mapper[.x], .cols = all_of(pathway_cols))
  
  # 3. 儲存檔案
  out_path <- paste0("inputs/", type_name, "_pathway_table.csv")
  write_csv(pat, out_path)
  
  cat(paste0(out_path, ": "), nrow(pat), "pathways\n")
}

# --- 執行部分 ---

# 處理 EC
process_and_save(pathways_ec, "ec", mapper)

# 處理 MetaCyc
process_and_save(pathways_metacyc, "metacyc", mapper)


# ── 6. 輸出 CSV ───────────────────────────────────────────────────────────────
dir.create(path = "inputs/")
write_csv(asv_table,                          "inputs/asv_table.csv")
write_csv(taxonomy,                           "inputs/taxonomy.csv")
write_csv(metadata %>% select(-analysis_id),  "inputs/metadata.csv")

cat("\n✓ 完成：asv_table.csv / taxonomy.csv / metadata.csv\n")


# ── 7. Sanity check ───────────────────────────────────────────────────────────

cat("\n── Sanity checks ───────────────────────────────────────────────────────\n")

samples_asv  <- names(asv_table)[-1]
samples_meta <- metadata$SampleID

cat(sprintf("ASV table 樣本數  : %d\n", length(samples_asv)))
cat(sprintf("Metadata 樣本數   : %d\n", length(samples_meta)))
cat(sprintf("完全對應          : %d\n", length(intersect(samples_asv, samples_meta))))

only_asv  <- setdiff(samples_asv,  samples_meta)
only_meta <- setdiff(samples_meta, samples_asv)
if (length(only_asv)  > 0) warning("只在 ASV table：",  paste(only_asv,  collapse = ", "))
if (length(only_meta) > 0) warning("只在 metadata：",   paste(only_meta, collapse = ", "))

stopifnot("ASV_ID 不一致" = identical(asv_table$ASV_ID, taxonomy$ASV_ID))

cat(sprintf("Read depth 範圍   : %s – %s\n",
            formatC(min(colSums(asv_table[, -1])), format = "d", big.mark = ","),
            formatC(max(colSums(asv_table[, -1])), format = "d", big.mark = ",")))

cat("\nSALT score（treatment × timepoint）：\n")
metadata %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    mean_SALT = round(mean(SALT_score, na.rm = TRUE), 1),
    n         = n(),
    .groups   = "drop"
  ) %>%
  print()

cat("\n✓ 準備完成，執行 source('alopecia_fmt_16S_analysis.R') 開始分析。\n")