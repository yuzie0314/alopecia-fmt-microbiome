# =============================================================================
#  Package requirements — Alopecia FMT Microbiome Pipeline
#  R 4.6.0 / Bioconductor 3.23
#
#  Run once before running the pipeline:
#    Rscript scripts/packages.R
#
#  For containerised alternatives (conda / Docker / Singularity), see INSTALL.md
# =============================================================================

# ── CRAN ──────────────────────────────────────────────────────────────────────
cran_packages <- c(
  # Data wrangling
  "tidyverse",     # dplyr, ggplot2, tidyr, readr, purrr, stringr, forcats
  "tibble",
  "readxl",        # Step 0: parse Excel metadata sheets
  "broom",

  # Visualisation
  "ggrepel",
  "patchwork",
  "RColorBrewer",
  "scales",
  "Cairo",         # PDF device for high-res output (Step 2, optional)

  # Community ecology / statistics
  "vegan",         # PERMANOVA (adonis2), rarefaction curves
  "lme4",          # Linear mixed models (alpha diversity)
  "lmerTest",      # LMM Satterthwaite p-values
  "emmeans",       # Estimated marginal means + post-hoc contrasts
  "rstatix"        # Kruskal-Wallis, Friedman, Wilcoxon helpers
)

# ── Bioconductor ──────────────────────────────────────────────────────────────
bioc_packages <- c(
  # Core microbiome
  "phyloseq",             # ASV / taxonomy / metadata S4 container
  "DESeq2",               # Count-based differential abundance
  "Maaslin2",             # Mixed-effects linear model DA
  "ANCOMBC",              # Compositional DA — use ancombc2() from v2.x

  # Sensitivity / alternative DA (Step 5)
  "lefser",               # LEfSe (v1.22 API: relativeAb + relab= param)
  "ALDEx2",               # ANOVA-like compositional DA (v1.44)

  # Functional enrichment (Steps 6–7)
  "fgsea",                # Fast gene/pathway set enrichment
  "KEGGREST",             # KEGG REST API — EC → KO → module mapping
  "SummarizedExperiment", # Required by lefser

  # Optional
  "ssizeRNA"              # Sample size power analysis (Step 3 only)
)

# ── Install ────────────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.23", ask = FALSE)

new_cran <- cran_packages[!cran_packages %in% installed.packages()[, "Package"]]
new_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[, "Package"]]

if (length(new_cran) > 0) {
  message("Installing CRAN packages: ", paste(new_cran, collapse = ", "))
  install.packages(new_cran, repos = "https://cloud.r-project.org")
}
if (length(new_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(new_bioc, collapse = ", "))
  BiocManager::install(new_bioc, ask = FALSE)
}

cat("\nAll packages installed. R", R.version$major, R.version$minor,
    "/ Bioc", as.character(BiocManager::version()), "\n")
