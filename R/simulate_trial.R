# Simulate one two-arm trial before any data cut (spec sections 2.3-2.5).
# Times in months; calendar time 0 = study activation.

# Annual dropout proportion -> monthly exponential rate.
dropout_rate <- function(annual_prop) -log(1 - annual_prop) / 12

# Treatment-arm PWE: control hazard x HR(t), with HR(t) = 1 before `delay` and `hr` after.
treatment_hazard <- function(cuts, rate, hr, delay) {
  new_cuts <- sort(unique(c(cuts, delay)))
  new_rate <- rate[findInterval(new_cuts, cuts)] * ifelse(new_cuts >= delay, hr, 1)
  list(cuts = new_cuts, rate = new_rate)
}

# Inverse of the PWE cumulative hazard: T = H^{-1}(E), so E ~ Exp(1) gives T ~ PWE.
pwe_inverse <- function(E, cuts, rate) {
  H_at_cuts <- c(0, cumsum(rate[-length(rate)] * diff(cuts)))
  k <- findInterval(E, H_at_cuts)
  cuts[k] + (E - H_at_cuts[k]) / rate[k]
}

simulate_trial <- function(N, enrol_rate, dropout_annual, hr, delay, control) {
  stopifnot(N %% 2 == 0)
  entry <- cumsum(rexp(N, enrol_rate))
  trt <- sample(rep(0:1, N / 2))
  trt_haz <- treatment_hazard(control$cuts, control$rate, hr, delay)
  E <- rexp(N)
  death_time <- ifelse(trt == 1,
                       pwe_inverse(E, trt_haz$cuts, trt_haz$rate),
                       pwe_inverse(E, control$cuts, control$rate))
  dropout_time <- if (dropout_annual > 0) rexp(N, dropout_rate(dropout_annual)) else rep(Inf, N)
  data.frame(entry = entry, trt = trt, death_time = death_time, dropout_time = dropout_time)
}
