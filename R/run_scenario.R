# Scenario grid and replicate runner (spec sections 4, 7, 9).

MASTER_SEED <- 20261029
PATTERNS <- data.frame(pattern = c("a_null", "b_ph", "c_delay6", "d_delay12"),
                       hr_on = c(FALSE, TRUE, TRUE, TRUE), delay = c(0, 0, 6, 12))

# Scenario k is row k; the first factor varies fastest.
scenario_grid <- function(hr = 0.69) {
  g <- expand.grid(pattern = PATTERNS$pattern, N = c(400, 500, 600, 700, 800, 900, 1000, 1200),
                   enrol = c(15, 30), dropout = c(0, 0.05, 0.15), stringsAsFactors = FALSE)
  g <- merge(g, PATTERNS, by = "pattern", sort = FALSE)
  g <- g[order(g$dropout, g$enrol, g$N, match(g$pattern, PATTERNS$pattern)), ]
  g$hr <- ifelse(g$hr_on, hr, 1)
  g$reps <- ifelse(g$pattern == "a_null", 10000, 2000)
  g$scenario_id <- seq_len(nrow(g))
  g$seed <- MASTER_SEED + g$scenario_id
  rownames(g) <- NULL
  g[, c("scenario_id", "pattern", "N", "enrol", "dropout", "hr", "delay", "reps", "seed")]
}

one_replicate <- function(scen, D, cap, control) {
  trial <- simulate_trial(scen$N, scen$enrol, scen$dropout, scen$hr, scen$delay, control)
  cut <- cut_data(trial, D, cap)
  res <- analyze_trial(cut$data)
  data.frame(as.list(res), t_star = cut$t_star, trigger = cut$trigger,
             deaths = cut$deaths, n_enrolled = cut$n_enrolled)
}

# `scen` is one row of scenario_grid(). Replicate streams come from furrr's L'Ecuyer-CMRG seeding,
# so results don't depend on the number of workers.
run_scenario <- function(scen, D, cap, control, reps = scen$reps, seed = scen$seed) {
  out <- furrr::future_map_dfr(
    seq_len(reps),
    function(i) one_replicate(scen, D, cap, control),
    .options = furrr::furrr_options(seed = seed)
  )
  out$rep <- seq_len(reps)
  out$scenario_id <- scen$scenario_id
  out
}
