# Stage 2: analytic planning calculations that justify the scenario grid in docs/design_spec.md.
# Output: results/design_calcs.rds

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
cm <- readRDS("results/control_model.rds")
stage1 <- readRDS("results/stage1_summary.rds")

HR_PH  <- round(stage1$cox_hr[["hr"]], 2)
D      <- ceiling(schoenfeld_events(HR_PH, alpha_one_sided = 0.025, power = 0.9))
CAP    <- 84
DUR_MAX <- 72

grid <- expand.grid(N = c(400, 500, 600, 700, 800, 900, 1000, 1200), enrol = c(15, 30),
                    dropout = c(0, 0.05, 0.15), delay = c(0, 6, 12))
grid$accrual_months <- grid$N / grid$enrol
grid$t_D <- mapply(function(N, e, dr, dl)
  time_to_events(D, N, e, cm$cuts, cm$rate, HR_PH, dl, dropout_rate(dr)),
  grid$N, grid$enrol, grid$dropout, grid$delay)
grid$deaths_at_cap <- mapply(function(N, e, dr, dl)
  expected_deaths(CAP, N, e, cm$cuts, cm$rate, HR_PH, dl, dropout_rate(dr)),
  grid$N, grid$enrol, grid$dropout, grid$delay)
grid$analysis_time <- pmin(grid$t_D, CAP)
grid$deaths_at_analysis <- pmin(D, grid$deaths_at_cap)
# Planning approximation only (PH): Schoenfeld power at the expected number of deaths analysed.
grid$approx_power_ph <- schoenfeld_power(HR_PH, grid$deaths_at_analysis)

# Monte Carlo SE targets.
mcse <- function(p, R) sqrt(p * (1 - p) / R)
mcse_table <- data.frame(
  quantity = c("power at 0.5 (worst case)", "power at 0.8", "type I error at 0.025"),
  reps = c(2000, 2000, 10000),
  mcse = c(mcse(0.5, 2000), mcse(0.8, 2000), mcse(0.025, 10000))
)
type1_threshold <- 0.025 + 2 * mcse(0.025, 10000)

design_calcs <- list(hr_ph = HR_PH, D = D, cap = CAP, dur_max = DUR_MAX,
                     d_unrounded = schoenfeld_events(HR_PH), grid = grid,
                     mcse_table = mcse_table, type1_threshold = type1_threshold)
saveRDS(design_calcs, "results/design_calcs.rds")

cat(sprintf("PH HR = %.2f; Schoenfeld D = %.1f -> %d; cap = %d months\n",
            HR_PH, design_calcs$d_unrounded, D, CAP))
print(mcse_table, digits = 3)
cat(sprintf("Type I error eligibility threshold: %.4f\n", type1_threshold))
show <- function(g) print(g[, c("N", "enrol", "dropout", "delay", "accrual_months", "t_D",
                                "deaths_at_cap", "approx_power_ph")], digits = 3, row.names = FALSE)
cat("\nPH (delay 0):\n"); show(subset(grid, delay == 0))
cat("\nBase case neighbourhood (delay 6, dropout 5%):\n"); show(subset(grid, delay == 6 & dropout == 0.05))
cat("\nDelay 12, dropout 5%:\n"); show(subset(grid, delay == 12 & dropout == 0.05))
