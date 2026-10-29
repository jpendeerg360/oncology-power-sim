# Single-design calculator used by the exploratory Shiny app (app/app.R).
# Not part of the pre-specified study: it reuses the study's simulator and analyses for one
# user-chosen design, and puts analytic approximations alongside.

# One design as a scenario_grid()-style row, so run_scenario() can be called unchanged.
design_scenario <- function(N, enrol, dropout, hr, delay, reps, seed) {
  pattern <- if (hr == 1) "a_null" else if (delay == 0) "b_ph" else "delay"
  data.frame(scenario_id = 0L, pattern = pattern, N = N, enrol = enrol, dropout = dropout,
             hr = hr, delay = delay, reps = reps, seed = seed, stringsAsFactors = FALSE)
}

# Simulated operating characteristics for one design. With the seed and reps of a grid scenario,
# this reproduces that scenario's row of results/scenario_summary.csv.
simulate_design <- function(N, enrol, dropout, hr, delay, D, cap, control,
                            reps = 1000, seed = MASTER_SEED) {
  scen <- design_scenario(N, enrol, dropout, hr, delay, reps, seed)
  out <- run_scenario(scen, D, cap, control)
  cbind(scen, summarise_scenario(out, scen$pattern, hr))
}

# gsDesign2 inputs for the study's model: PWE control hazard, HR(t) = 1 before `delay` and `hr`
# after, exponential dropout, uniform accrual at `enrol` per month.
gs_rates <- function(N, enrol, dropout, hr, delay, control) {
  cuts <- sort(unique(c(control$cuts, delay)))
  list(
    enroll = gsDesign2::define_enroll_rate(duration = N / enrol, rate = enrol),
    fail = gsDesign2::define_fail_rate(
      duration = c(diff(cuts), 100), fail_rate = control$rate[findInterval(cuts, control$cuts)],
      dropout_rate = dropout_rate(dropout), hr = ifelse(cuts >= delay, hr, 1))
  )
}

# Analytic approximations at the expected analysis time (D deaths or the cap, whichever is first):
# average-HR log-rank and FH(0,1) power from gsDesign2, plus Schoenfeld at D deaths and the input HR.
analytic_design <- function(N, enrol, dropout, hr, delay, D, cap, control, alpha = ALPHA) {
  r <- gs_rates(N, enrol, dropout, hr, delay, control)
  events_at_cap <- gsDesign2::ahr(r$enroll, r$fail, total_duration = cap)$event
  t_an <- if (events_at_cap < D) cap else
    gsDesign2::expected_time(r$enroll, r$fail, target_event = D, interval = c(0.1, cap))$time
  lr <- gsDesign2::fixed_design_ahr(r$enroll, r$fail, alpha = alpha, study_duration = t_an)$analysis
  fh <- gsDesign2::fixed_design_fh(alpha = alpha, study_duration = t_an, enroll_rate = r$enroll,
                                   fail_rate = r$fail, rho = 0, gamma = 1)$analysis
  data.frame(time = t_an, events = lr$event, ahr = lr$ahr,
             power_lr_ahr = lr$power, power_fh_ahr = fh$power,
             power_schoenfeld = if (hr == 1) alpha else schoenfeld_power(hr, D, alpha))
}
