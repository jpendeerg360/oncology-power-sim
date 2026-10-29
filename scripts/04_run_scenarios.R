# Stage 4: run every scenario in the pre-specified grid (docs/design_spec.md).
# Outputs:
#   results/reps/scenario_###.rds   per-replicate results, one file per scenario
#   results/scenario_grid.csv       scenario definitions, replicate counts and seeds
#   results/scenario_summary.csv    operating characteristics recomputed from the replicate files
#   results/run_log.txt             timing, seeds and sessionInfo()

suppressPackageStartupMessages({
  library(future)
  library(dplyr)
})
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

control <- readRDS("results/control_model.rds")
dc <- readRDS("results/design_calcs.rds")
D <- dc$D; CAP <- dc$cap; HR <- dc$hr_ph

grid <- scenario_grid(HR)
write.csv(grid, "results/scenario_grid.csv", row.names = FALSE)
dir.create("results/reps", showWarnings = FALSE)

N_WORKERS <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = N_WORKERS)

rep_file <- function(id) sprintf("results/reps/scenario_%03d.rds", id)

t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  scen <- grid[i, ]
  if (file.exists(rep_file(scen$scenario_id))) next   # resume support
  saveRDS(run_scenario(scen, D, CAP, control), rep_file(scen$scenario_id))
}
run_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
plan(sequential)

summary_tbl <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  reps <- readRDS(rep_file(grid$scenario_id[i]))
  stopifnot(nrow(reps) == grid$reps[i])
  cbind(grid[i, ], summarise_scenario(reps, grid$pattern[i], grid$hr[i]))
}))
write.csv(summary_tbl, "results/scenario_summary.csv", row.names = FALSE)

writeLines(c(
  sprintf("Run finished: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Run time (minutes, this invocation): %.2f", run_minutes),
  sprintf("Workers: %d", N_WORKERS),
  sprintf("Master seed: %d; scenario k uses furrr_options(seed = %d + k)", MASTER_SEED, MASTER_SEED),
  sprintf("D = %d, cap = %d months, PH HR = %.2f", D, CAP, HR),
  sprintf("Scenarios: %d; replicates: %d", nrow(grid), sum(grid$reps)),
  "", capture.output(sessionInfo())
), "results/run_log.txt")

cat(sprintf("Done: %d scenarios, %d replicates, %.2f minutes\n",
            nrow(grid), sum(grid$reps), run_minutes))
