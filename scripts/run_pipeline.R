# =============================================================================
#  run_pipeline.R  — Alopecia Oral FMT full analysis pipeline
#  執行順序：Load → Alpha → Beta → DA → Pairwise → ANCOM/LEfSe → Pathway → Enrichment
# =============================================================================

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

t0 <- proc.time()
step_times <- list()

run_step <- function(label, script_path) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("  STEP: %s\n", label))
  cat(sprintf("  %s\n", format(Sys.time(), "%H:%M:%S")))
  cat(sprintf("========================================\n"))
  t_start <- proc.time()
  tryCatch(
    source(script_path, local = FALSE),
    error = function(e) {
      cat(sprintf("\n[ERROR] %s failed: %s\n", label, conditionMessage(e)))
    }
  )
  elapsed <- (proc.time() - t_start)[["elapsed"]]
  step_times[[label]] <<- elapsed
  cat(sprintf("\n  Done in %.1f s\n", elapsed))
}

# ── 0. Load data ──────────────────────────────────────────────────────────────
run_step("0 Load Data",      "scripts/0_1-loadData.R")

# ── 1. Alpha diversity ────────────────────────────────────────────────────────
run_step("1 Alpha",          "scripts/1-alpha.R")

# ── 2. Beta diversity ─────────────────────────────────────────────────────────
run_step("2 Beta",           "scripts/2-beta.R")

# ── 3. Differential abundance (taxonomy) ─────────────────────────────────────
run_step("3 DA taxonomy",    "scripts/3-da.R")

# ── 4. Pairwise timepoint DA ──────────────────────────────────────────────────
run_step("4 Pairwise DA",    "scripts/5-pairwise_da.R")

# ── 5. ANCOM-BC2 + LEfSe + ALDEx2 ────────────────────────────────────────────
run_step("5 ANCOM/LEfSe",    "scripts/ancom_lefse_aldex.R")

# ── 6. Pathway analysis ───────────────────────────────────────────────────────
run_step("6 Pathway",        "scripts/4-pathway_analysis.R")

# ── 7. Functional enrichment ─────────────────────────────────────────────────
run_step("7 Enrichment",     "scripts/functional_enrichment.R")

# ── 8. Improvement analyses (optional) ───────────────────────────────────────
# run_step("8 Improvements",   "scripts/improvement_analyses.R")

# ── Summary ───────────────────────────────────────────────────────────────────
total <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("\n\n====== PIPELINE COMPLETE (%d min %d sec) ======\n",
            floor(total / 60), round(total %% 60)))
for (nm in names(step_times))
  cat(sprintf("  %-25s %5.0f s\n", nm, step_times[[nm]]))
