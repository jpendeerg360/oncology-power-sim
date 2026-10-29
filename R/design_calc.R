# Analytic design calculations: Schoenfeld events and expected deaths over calendar time.
# Used to choose the scenario grid; the simulation itself does not depend on these.

# Required deaths for a 1:1 log-rank test (Schoenfeld).
schoenfeld_events <- function(hr, alpha_one_sided = 0.025, power = 0.9) {
  4 * (qnorm(1 - alpha_one_sided) + qnorm(power))^2 / log(hr)^2
}

# Approximate log-rank power given D deaths under PH (Schoenfeld).
schoenfeld_power <- function(hr, events, alpha_one_sided = 0.025) {
  pnorm(sqrt(events / 4) * abs(log(hr)) - qnorm(1 - alpha_one_sided))
}

# Probability that a patient followed for `u` months dies before dropping out,
# with control PWE hazard multiplied by `hr` after `delay` months.
death_prob <- function(u, cuts, rate, hr = 1, delay = 0, drop = 0, step = 0.25) {
  if (u <= 0) return(0)
  s <- seq(0, u, by = step)
  h <- rate[findInterval(s, cuts)] * ifelse(s >= delay, hr, 1)
  H <- c(0, cumsum((h[-1] + h[-length(h)]) / 2 * diff(s)))
  dens <- h * exp(-H - drop * s)
  sum((dens[-1] + dens[-length(dens)]) / 2 * diff(s))
}

# Expected total deaths at calendar time T (months), uniform accrual of N at `enrol` per month,
# 1:1 randomization.
expected_deaths <- function(T, N, enrol, cuts, rate, hr, delay, drop) {
  acc_end <- N / enrol
  a <- seq(0, min(T, acc_end), length.out = 200)
  if (length(unique(a)) < 2) return(0)
  f <- function(arm_hr) vapply(T - a, death_prob, numeric(1), cuts = cuts, rate = rate,
                               hr = arm_hr, delay = delay, drop = drop)
  per_entry <- 0.5 * (f(1) + f(hr))
  enrol * sum((per_entry[-1] + per_entry[-length(per_entry)]) / 2 * diff(a))
}

# Calendar time at which expected deaths reach D (Inf if not reached by `t_max`).
time_to_events <- function(D, N, enrol, cuts, rate, hr, delay, drop, t_max = 240) {
  g <- function(T) expected_deaths(T, N, enrol, cuts, rate, hr, delay, drop) - D
  if (g(t_max) < 0) return(Inf)
  uniroot(g, c(1, t_max), tol = 0.05)$root
}
