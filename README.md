# Alopecia Oral FMT — Gut Microbiome Analysis Pipeline

An end-to-end R pipeline for analysing 16S rRNA amplicon sequencing and functional pathway data from a clinical trial of oral faecal microbiota transplantation (FMT) in alopecia areata patients.

> **Note**: This repository contains only analysis code and simulated demonstration data.  
> The original clinical trial data is confidential and not included.

---

## Study Design

| Feature | Detail |
|---------|--------|
| Design | Double-blind RCT pilot |
| n | 10 patients (5 FMT / 5 Placebo) |
| Timepoints | Baseline, 2M post-FMT, 6M post-FMT |
| Primary outcome | SALT score (Severity of Alopecia Tool) |
| Sequencing | 16S rRNA V3-V4 amplicon |
| Functional data | MetaCyc pathways + EC numbers (HUMAnN3) |

---

## Repository Structure

```
.
├── scripts/
│   ├── simulate_demo_data.R   # Generate synthetic demo inputs
│   ├── run_pipeline.R         # Master pipeline runner
│   ├── 0_1-loadData.R         # Data loading & QC
│   ├── 1-alpha.R              # Alpha diversity
│   ├── 2-beta.R               # Beta diversity (Bray-Curtis + Aitchison)
│   ├── 3-da.R                 # Differential abundance (DESeq2 + MaAsLin2)
│   ├── 4-pathway_analysis.R   # Pathway beta diversity + MaAsLin2
│   ├── 5-pairwise_da.R        # Pairwise timepoint DA
│   ├── ancom_lefse_aldex.R    # ANCOM-BC2 + LEfSe + ALDEx2
│   ├── functional_enrichment.R # fgsea + KEGG ORA
│   └── improvement_analyses.R # Publication figures + statistical improvements
├── sim_inputs/                # Synthetic demonstration data (auto-generated)
│   ├── metadata.csv
│   ├── asv_table.csv
│   ├── taxonomy.csv
│   ├── metacyc_pathway_table.csv
│   └── ec_pathway_table.csv
├── CLAUDE.md                  # Full analysis documentation
└── README.md
```

---

## Quick Start (Demo Mode)

### 1. Prerequisites

```r
# Install required packages
install.packages(c("tidyverse","vegan","ggrepel","patchwork","lme4","emmeans"))

if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "phyloseq","DESeq2","ANCOMBC","Maaslin2","lefser","ALDEx2",
  "fgsea","SummarizedExperiment","Biobase"
))
```

**R version**: 4.6.0 | **Bioconductor**: 3.23

### 2. Generate simulated data

```r
source("scripts/simulate_demo_data.R")
# Creates: sim_inputs/metadata.csv, asv_table.csv, taxonomy.csv,
#          metacyc_pathway_table.csv, ec_pathway_table.csv
```

### 3. Update data paths in `0_1-loadData.R`

Change the input file paths from `inputs/` to `sim_inputs/`:

```r
asv_raw     <- read.csv("sim_inputs/asv_table.csv",  row.names=1, check.names=FALSE)
tax_raw     <- read.csv("sim_inputs/taxonomy.csv",   stringsAsFactors=FALSE)
meta_raw    <- read.csv("sim_inputs/metadata.csv",   stringsAsFactors=FALSE)
metacyc_raw <- load_pathway("sim_inputs/metacyc_pathway_table.csv", "MetaCyc")
ec_raw      <- load_pathway("sim_inputs/ec_pathway_table.csv",      "EC")
```

Also update the metadata column format (sim data uses `baseline/month2/month6`):

```r
# In 0_1-loadData.R, the timepoint mapping section:
meta_raw$Timepoint <- recode(meta_raw$Timepoint,
  baseline = "Baseline",
  month2   = "2M post-FMT",
  month6   = "6M post-FMT"
)
```

### 4. Run the pipeline

```r
setwd("/path/to/this/repo")
source("scripts/run_pipeline.R")
```

### 5. (Optional) Functional enrichment extras

Steps 6–7 require MetaCyc and KEGG mapping files:

| File | Source |
|------|--------|
| `map_metacyc-pwy_name.txt.gz` | [MetaCyc downloads](https://metacyc.org/downloads.shtml) |
| `map_metacyc-pwy_lineage.tsv` | [MetaCyc downloads](https://metacyc.org/downloads.shtml) |
| `ko_ec_cache.rds` | Generate via KEGG API (`KEGGREST` package) |
| `ko_module_cache.rds` | Generate via KEGG API |
| `ko_pathway_cache.rds` | Generate via KEGG API |

Update paths in `functional_enrichment.R` accordingly.

---

## Analysis Pipeline

```
Raw 16S ASVs + MetaCyc pathways + EC numbers
        │
Step 0  │  Data loading, QC, metadata merge
Step 1  │  Alpha diversity (Observed, Chao1, ACE, Shannon, Pielou's J)
Step 2  │  Beta diversity (Bray-Curtis PCoA + Aitchison CLR PCoA)
Step 3  │  Differential abundance — taxonomy
        │    DESeq2 (count-based) + MaAsLin2 (mixed model) + ANCOM-BC2
Step 4  │  Pairwise timepoint DA (BL→2M, 2M→6M, BL→6M)
Step 5  │  ANCOM-BC2 + LEfSe + ALDEx2 (taxonomy + pathway)
Step 6  │  Pathway analysis (MetaCyc + EC, MaAsLin2)
Step 7  │  Functional enrichment (fgsea MetaCyc categories, KEGG ORA)
        │
        └→ Clinical integration: ΔPC1 (microbiome shift) ~ ΔSALT (hair regrowth)
```

### Repeated-measures design

All analyses correctly account for the repeated-measures structure (3 timepoints × 10 patients):
- LMM with PatientID random effect (alpha diversity)
- PERMANOVA with `how(blocks = PatientID)` (beta diversity)
- MaAsLin2 with PatientID random effect (DA)

---

## Key Statistical Findings (Simulated Data)

> Results below are from **simulated data** and do not reflect real clinical outcomes.

| Analysis | Finding |
|---------|---------|
| Beta diversity | Aitchison ΔPC1 ~ ΔSALT: Spearman ρ (reported after pipeline run) |
| Top DA taxon | *Klebsiella* ↓ in FMT (planted signal in simulation) |
| Pathway enrichment | Energy-Metabolism fgsea NES (reported after pipeline run) |

---

## Technical Notes

### Known version-specific behaviours

| Package | Version | Note |
|---------|---------|------|
| ANCOMBC | 2.14.0 | `rand_formula = NULL`; 2-group columns use `(Intercept)` naming |
| lefser | 1.22.0 | Requires `relativeAb()` first; use `droplevels()` on grouping factor |
| ALDEx2 | 1.44.0 | `aldex.ttest()` doesn't include effect; call `aldex.effect()` separately |

### Common issues

- **`select` namespace conflict**: Use `select <- dplyr::select` at top of each script
- **phyloseq S4 stripping**: Use `matrix(as.integer(as.vector(otu_table(ps))), ...)` for base R functions
- **MaAsLin2 `reference` parameter**: Must match actual factor levels in the subset

---

## Skills Demonstrated

- Repeated-measures microbiome study design (n=10, 3 timepoints)
- Multi-method differential abundance: DESeq2, MaAsLin2, ANCOM-BC2, LEfSe, ALDEx2
- Compositional data analysis (CLR transformation, Aitchison geometry)
- Functional metagenomics: pathway enrichment (fgsea), KEGG ORA
- Clinical correlation: microbiome community shift ~ clinical outcome (SALT score)
- R pipeline engineering: modular scripts, error handling, >130 min runtime management

---

## License

Code: MIT  
Simulated data: freely reusable (no real patient data)

---

*Analysis pipeline developed for a clinical trial of oral FMT in alopecia areata.*  
*R 4.6.0 / Bioconductor 3.23 / 2026*
