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
#                       Timepoint values: "baseline", "month2", "month6"
#                       Treatment values: "FMT", "placebo"
# =============================================================================


# ── 0. Install & load packages ───────────────────────────────────────────────

packages <- c(
  "phyloseq", "vegan", "DESeq2", "ggplot2", "dplyr", "tidyr",
  "readr", "tibble", "purrr", "forcats",
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

load_pathway <- function(path, type = "metacyc") {
  df <- read_csv(path, show_col_types = FALSE)
  
  # 第一欄是 pathway/EC ID
  id_col <- names(df)[1]
  mat <- df %>%
    column_to_rownames(id_col) %>%
    as.matrix()
  
  cat(sprintf("[%s] %d pathways × %d samples\n", type, nrow(mat), ncol(mat)))
  return(mat)
}

asv_raw  <- read.csv("inputs/asv_table.csv",  row.names = 1, check.names = FALSE)
tax_raw  <- read.csv("inputs/taxonomy.csv",   stringsAsFactors = FALSE)
meta_raw <- read.csv("inputs/metadata.csv",   stringsAsFactors = FALSE)
metacyc_raw <- load_pathway("inputs/metacyc_pathway_table.csv", "MetaCyc")
ec_raw      <- load_pathway("inputs/ec_pathway_table.csv",      "EC")

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

# percentage-based filtering
filter_ps_by_prev <- function(percentage = 0.1, ps_obj = ps){
  prevalence_threshold <- percentage * nsamples(ps_obj)
  ps_filt_prev <- filter_taxa(ps_obj, function(x) sum(x > 0) >= prevalence_threshold, TRUE)
  return(ps_filt_prev)
}


# Rarefy to even depth for alpha/beta diversity
# set.seed(42)
rarefaction_depth <- floor(min(sample_sums(ps_filt)) * 0.9)
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
