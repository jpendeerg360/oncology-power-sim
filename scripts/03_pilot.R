# Stage 3: pilot run on a few scenarios, validation summaries and runtime projection.
# Output: results/pilot_summary.rds, results/runtime_projection.rds

suppressPackageStartupMessages({
  library(future)
  library(dplyr)
})
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

control <- readRDS("results/control_model.rds")
dc <- readRDS("results/design_calcs.rds")
D <- dc$D; CAP <- dc$cap; HR <- dc$hr_ph

N_WORKERS <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = N_WORKERS)

# Pilot scenarios come from the main grid but use a separate seed range, so the pilot is an
# independent check. Null scenarios get more replicates so type I error can be judged.
PILOT_SEED_OFFSET <- 500000
grid <- scenario_grid(HR)
pick <- function(pattern, N, enrol, dropout) {
  grid[grid$pattern == pattern & grid$N == N & grid$enrol == enrol & grid$dropout == dropout, ]
}
pilot <- bind_rows(
  pick("b_ph", 1200, 30, 0),       # PH, cap essentially never binds: Schoenfeld check
  pick("b_ph", 900, 30, 0.05),
  pick("b_ph", 500, 15, 0.15),     # PH, cap almost always binds
  pick("d_delay12", 900, 30, 0.05),  # base case neighbourhood
  pick("a_null", 900, 30, 0.05),
  pick("a_null", 500, 15, 0.15)
)
pilot$reps <- ifelse(pilot$pattern == "a_null", 5000, 500)
pilot$seed <- pilot$seed + PILOT_SEED_OFFSET

t0 <- Sys.time()
pilot_reps <- lapply(seq_len(nrow(pilot)), function(i) run_scenario(pilot[i, ], D, CAP, control))
pilot_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

pilot_summary <- bind_cols(
  pilot[, c("scenario_id", "pattern", "N", "enrol", "dropout", "hr", "delay")],
  bind_rows(Map(function(r, i) summarise_scenario(r, pilot$pattern[i], pilot$hr[i]),
                pilot_reps, seq_len(nrow(pilot))))
)
pilot_summary$schoenfeld_power <- ifelse(
  pilot_summary$pattern == "b_ph", schoenfeld_power(HR, pilot_summary$deaths_mean), NA)

# Control-arm KM in the PH no-cap pilot scenario vs the calibrated curve, using untruncated
# simulated trials (same generator as the pilot).
set.seed(MASTER_SEED)
ctrl_times <- unlist(lapply(1:50, function(i) {
  tr <- simulate_trial(1200, 30, 0, 1, 0, control); tr$death_time[tr$trt == 0]
}))
km <- survival::survfit(survival::Surv(ctrl_times) ~ 1)
t_check <- c(6, 12, 24, 36, 48, 60, 72, 84)
km_check <- data.frame(t = t_check, km = summary(km, times = t_check)$surv,
                       pwe = pwe_surv(t_check, control$cuts, control$rate))

# Runtime projection: time a fixed batch per N (cost scales with N) and extrapolate to the grid.
timing <- lapply(c(400, 800, 1200), function(n) {
  scen <- pick("b_ph", n, 30, 0.05)
  t0 <- Sys.time()
  run_scenario(scen, D, CAP, control, reps = 6000, seed = 1)
  data.frame(N = n, secs_per_rep = as.numeric(difftime(Sys.time(), t0, units = "secs")) / 6000)
})
timing <- bind_rows(timing)
cost_fit <- lm(secs_per_rep ~ N, data = timing)
grid$pred_secs <- grid$reps * predict(cost_fit, newdata = grid)
runtime <- list(n_workers = N_WORKERS, timing = timing,
                projected_minutes = sum(grid$pred_secs) / 60,
                total_reps = sum(grid$reps), pilot_secs = pilot_secs)

saveRDS(list(summary = pilot_summary, km_check = km_check, reps = pilot_reps),
        "results/pilot_summary.rds")
saveRDS(runtime, "results/runtime_projection.rds")
plan(sequential)

# ---- Console summary ----------------------------------------------------------
options(width = 150)
print(pilot_summary |> transmute(
  pattern, N, enrol, dropout, R,
  power_lr = round(power_lr, 4), se_lr = round(mcse_power_lr, 4),
  power_fh = round(power_fh, 4), se_fh = round(mcse_power_fh, 4),
  schoenfeld = round(schoenfeld_power, 3),
  duration = round(duration_mean, 1), pct_cap = round(pct_cap, 3),
  deaths = round(deaths_mean, 1), loghr = round(loghr_mean, 3), cover = round(coverage, 3)))
cat("\nControl-arm KM vs calibrated PWE:\n"); print(round(km_check, 4))
cat(sprintf("\nWorkers: %d. Seconds per replicate (parallel wall time):\n", N_WORKERS))
print(timing, digits = 3)
cat(sprintf("Projected full grid: %s replicates, %.1f minutes\n",
            format(runtime$total_reps, big.mark = ","), runtime$projected_minutes))
