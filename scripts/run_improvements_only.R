# run_improvements_only.R
# 最小 runner：Load + Beta → improvement_analyses.R
# 省去重跑完整 136 min pipeline

setwd("C:/Users/User/OneDrive - AMILI Pte Ltd/AMILI/project/alopecia_2026May")

t0 <- proc.time()

run_step <- function(label, script) {
  cat(sprintf("\n=== %s (%s) ===\n", label, format(Sys.time(), "%H:%M:%S")))
  t1 <- proc.time()
  tryCatch(source(script, local = FALSE),
           error = function(e) cat(sprintf("[ERROR] %s: %s\n", label, e$message)))
  cat(sprintf("  Done: %.1f s\n", (proc.time() - t1)[["elapsed"]]))
}

run_step("0 Load Data",   "scripts/0_1-loadData.R")
run_step("2 Beta",        "scripts/2-beta.R")
run_step("8 Improvements","scripts/improvement_analyses.R")

total <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("\n=== DONE: %.0f s total ===\n", total))
