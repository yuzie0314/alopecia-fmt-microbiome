# =============================================================================
#  simulate_demo_data.R
#  Generate synthetic demonstration data for the Alopecia FMT analysis pipeline
#
#  Output: sim_inputs/
#    metadata.csv          - 10 simulated patients (5 FMT / 5 placebo)
#    asv_table.csv         - 200 simulated ASVs x 29 samples
#    taxonomy.csv          - ASV taxonomy (realistic genus names)
#    metacyc_pathway_table.csv - 150 pathways x 29 samples
#    ec_pathway_table.csv  - 300 EC numbers x 29 samples
#
#  NOTE: All data is completely synthetic. Patient IDs, SALT scores,
#        demographics, and microbiome abundances are randomly generated.
#        No real patient data is encoded here.
# =============================================================================

set.seed(2026)
library(tidyverse)

dir.create("sim_inputs", showWarnings = FALSE)

# =============================================================================
#  1. Patient metadata
# =============================================================================

patients <- data.frame(
  PatientID = paste0("SIM", sprintf("%02d", 1:10)),
  Sex       = c("Male","Female","Male","Male","Female","Female","Female","Male","Male","Male"),
  Treatment = c("FMT","placebo","FMT","placebo","placebo","FMT","FMT","placebo","placebo","FMT"),
  Age       = round(runif(10, 25, 65)),
  Race      = sample(c("Chinese","Indian","Malay"), 10, replace = TRUE,
                     prob = c(0.7, 0.15, 0.15)),
  stringsAsFactors = FALSE
)

# Simulate SALT scores: FMT group shows greater improvement on average
# Baseline: 20-90, Month2: slight drop, Month6: larger drop for FMT
gen_salt <- function(bl, trt, seed_offset) {
  set.seed(seed_offset)
  m2_drop <- ifelse(trt == "FMT", runif(1, 0.05, 0.40), runif(1, 0.0, 0.20))
  m6_drop <- ifelse(trt == "FMT", runif(1, 0.20, 0.80), runif(1, 0.05, 0.40))
  m2 <- round(max(0, bl * (1 - m2_drop)), 1)
  m6 <- round(max(0, bl * (1 - m6_drop)), 1)
  c(bl = bl, month2 = m2, month6 = m6)
}

bl_salt <- round(runif(10, 15, 90), 1)
salt_mat <- mapply(gen_salt, bl_salt, patients$Treatment, 1:10 + 100)

# Build long-format metadata
meta_rows <- list()
for (i in seq_len(nrow(patients))) {
  pt  <- patients$PatientID[i]
  trt <- patients$Treatment[i]
  for (tp in c("baseline","month2","month6")) {
    # SIM03 (FMT) has no baseline
    if (pt == "SIM03" && tp == "baseline") next
    v <- switch(tp, baseline="V1", month2="V2", month6="V3")
    salt_val <- switch(tp,
      baseline = salt_mat["bl", i],
      month2   = salt_mat["month2", i],
      month6   = salt_mat["month6", i]
    )
    abx <- if (tp == "month6") sample(c("Yes","No"), 1, prob=c(0.2, 0.8)) else "NA"
    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      PatientID  = pt,
      SampleID   = paste0(pt, "_", v),
      Sex        = patients$Sex[i],
      Treatment  = trt,
      Age        = patients$Age[i],
      Race       = patients$Race[i],
      SALT_score = salt_val,
      Timepoint  = tp,
      antibiotics = abx,
      stringsAsFactors = FALSE
    )
  }
}
meta <- do.call(rbind, meta_rows)
write.csv(meta, "sim_inputs/metadata.csv", row.names = FALSE, quote = FALSE)
cat(sprintf("metadata.csv: %d rows x %d cols\n", nrow(meta), ncol(meta)))

# =============================================================================
#  2. Taxonomy (200 ASVs, realistic genus names)
# =============================================================================

gut_genera <- c(
  "Bacteroides","Prevotella","Faecalibacterium","Blautia","Roseburia",
  "Lachnospira","Ruminococcus","Clostridium","Bifidobacterium","Lactobacillus",
  "Akkermansia","Streptococcus","Veillonella","Dialister","Megamonas",
  "Phocaeicola","Alistipes","Parabacteroides","Klebsiella","Escherichia",
  "Collinsella","Adlercreutzia","Eggerthella","Coriobacterium","Slackia",
  "Coprococcus","Dorea","Hungatella","Lachnoclostridium","Butyrivibrio",
  "Mediterraneibacter","Romboutsia","Erysipelatoclostridium","Faecalibaculum",
  "Holdemanella","Intestinibacter","Mogibacterium","Oribacterium",
  "Peptostreptococcus","Shigella"
)
families <- c(
  "Bacteroidaceae","Prevotellaceae","Ruminococcaceae","Lachnospiraceae",
  "Clostridiaceae","Bifidobacteriaceae","Lactobacillaceae","Akkermansiaceae",
  "Streptococcaceae","Veillonellaceae","Coriobacteriaceae","Erysipelotrichaceae",
  "Peptostreptococcaceae","Enterobacteriaceae"
)
phyla <- c(
  "Bacteroidota","Firmicutes","Actinobacteria","Proteobacteria",
  "Verrucomicrobiota","Fusobacteriota"
)

n_asv <- 200
asv_ids <- paste0("ASV", sprintf("%04d", seq_len(n_asv)))
genus_assign <- sample(gut_genera, n_asv, replace = TRUE)

taxonomy <- data.frame(
  ASV_ID  = asv_ids,
  Kingdom = "Bacteria",
  Phylum  = sapply(genus_assign, function(g) {
    if (g %in% c("Bacteroides","Prevotella","Phocaeicola","Alistipes","Parabacteroides"))
      "Bacteroidota"
    else if (g %in% c("Klebsiella","Escherichia","Shigella")) "Proteobacteria"
    else if (g %in% c("Bifidobacterium","Collinsella","Adlercreutzia","Eggerthella","Coriobacterium","Slackia"))
      "Actinobacteria"
    else if (g == "Akkermansia") "Verrucomicrobiota"
    else "Firmicutes"
  }),
  Class   = "Clostridia",
  Order   = "Clostridiales",
  Family  = sample(families, n_asv, replace = TRUE),
  Genus   = genus_assign,
  Species = NA_character_,
  stringsAsFactors = FALSE
)
write.csv(taxonomy, "sim_inputs/taxonomy.csv", row.names = FALSE, quote = FALSE)
cat(sprintf("taxonomy.csv: %d ASVs\n", n_asv))

# =============================================================================
#  3. ASV count table  (ASV x sample)
# =============================================================================

sample_ids <- meta$SampleID
n_samp <- length(sample_ids)

# Base read depth per sample
read_depth <- round(rnorm(n_samp, mean = 18000, sd = 3000))
read_depth <- pmax(read_depth, 5000)

# Dirichlet-multinomial: prevalence-weighted
prev_weights <- rexp(n_asv, rate = 2) + 0.01
prev_weights <- prev_weights / sum(prev_weights)

# Simulate counts
asv_counts <- matrix(0L, nrow = n_asv, ncol = n_samp,
                     dimnames = list(asv_ids, sample_ids))
for (j in seq_len(n_samp)) {
  alpha <- prev_weights * 50
  theta <- rgamma(n_asv, shape = alpha, rate = 1)
  theta <- theta / sum(theta)
  asv_counts[, j] <- rmultinom(1, read_depth[j], theta)[, 1]
}

# Plant a small FMT signal: Klebsiella lower at 6M in FMT group
kleb_rows <- which(taxonomy$Genus == "Klebsiella")
fmt_6m <- meta %>% filter(Treatment == "FMT", Timepoint == "month6") %>% pull(SampleID)
plc_6m <- meta %>% filter(Treatment == "placebo", Timepoint == "month6") %>% pull(SampleID)
if (length(kleb_rows) > 0) {
  asv_counts[kleb_rows, fmt_6m] <- round(asv_counts[kleb_rows, fmt_6m] * 0.15)
  asv_counts[kleb_rows, plc_6m] <- round(asv_counts[kleb_rows, plc_6m] * 1.20)
}

asv_df <- as.data.frame(asv_counts)
asv_df <- cbind(ASV_ID = rownames(asv_df), asv_df)
write.csv(asv_df, "sim_inputs/asv_table.csv", row.names = FALSE, quote = FALSE)
cat(sprintf("asv_table.csv: %d ASVs x %d samples\n", n_asv, n_samp))

# =============================================================================
#  4. MetaCyc pathway table  (150 pathways x samples)
# =============================================================================

metacyc_pathways <- c(
  "cis-vaccenate biosynthesis","L-methionine biosynthesis I","L-methionine biosynthesis II",
  "thiamine diphosphate biosynthesis II","NAD phosphorylation and transhydrogenation",
  "L-arginine degradation I (arginase pathway)","formate to 5,10-methyleneTHF",
  "(2Fe-2S) iron-sulfur cluster biosynthesis","L-serine degradation",
  "superpathway of glycolysis and the Entner-Doudoroff pathway",
  "TCA cycle (prokaryotic)","pentose phosphate pathway",
  "fatty acid beta-oxidation I","fatty acid biosynthesis initiation I",
  "acetyl-CoA fermentation to butyrate II","L-glutamate degradation V",
  "purine nucleotides de novo biosynthesis II","UMP biosynthesis I",
  "peptidoglycan biosynthesis I","lipid A-core biosynthesis",
  "adenosine nucleotides de novo biosynthesis","CMP-sialic acid biosynthesis",
  "cobalamin biosynthesis","riboflavin biosynthesis I","folate biosynthesis VI",
  paste0("MetaCyc-pathway-", sprintf("%03d", 26:150))
)
metacyc_pathways <- metacyc_pathways[1:150]

rpk_mat <- function(n_features, sids, base_mean = 5000, sd = 3000) {
  mat <- matrix(pmax(0, rnorm(n_features * length(sids),
                              mean = base_mean, sd = sd)),
                nrow = n_features, ncol = length(sids),
                dimnames = list(seq_len(n_features), sids))
  round(mat, 3)
}

mc_mat <- rpk_mat(150, sample_ids)
mc_df <- as.data.frame(mc_mat)
mc_df <- cbind(pathway = metacyc_pathways, mc_df)
write.csv(mc_df, "sim_inputs/metacyc_pathway_table.csv", row.names = FALSE, quote = FALSE)
cat(sprintf("metacyc_pathway_table.csv: 150 pathways x %d samples\n", n_samp))

# =============================================================================
#  5. EC pathway table  (300 ECs x samples)
# =============================================================================

ec_ids <- paste0(
  sample(1:6, 300, replace = TRUE), ".",
  sample(1:99, 300, replace = TRUE), ".",
  sample(1:99, 300, replace = TRUE), ".",
  sample(1:999, 300, replace = TRUE)
)
ec_ids <- paste0("EC:", ec_ids)

ec_mat <- rpk_mat(300, sample_ids, base_mean = 3000, sd = 2000)
ec_df <- as.data.frame(ec_mat)
ec_df <- cbind(pathway = ec_ids, ec_df)
write.csv(ec_df, "sim_inputs/ec_pathway_table.csv", row.names = FALSE, quote = FALSE)
cat(sprintf("ec_pathway_table.csv: 300 ECs x %d samples\n", n_samp))

cat("\n=== Simulation complete ===\n")
cat("Output: sim_inputs/\n")
cat("Run the pipeline with: setwd(); source('scripts/run_pipeline.R')\n")
cat("(update setwd path and input file paths in 0_1-loadData.R first)\n")