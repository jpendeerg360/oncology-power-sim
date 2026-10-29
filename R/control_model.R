# Piecewise-exponential (PWE) and Weibull helpers for the control-arm model.
# Time is in months throughout. `cuts` are interval start points, beginning at 0;
# the last interval is open-ended.

DAYS_PER_MONTH <- 365.25 / 12

# MLE of a PWE hazard: events / person-time in each interval.
pwe_fit <- function(time, status, cuts) {
  stopifnot(cuts[1] == 0, !is.unsorted(cuts))
  ends <- c(cuts[-1], Inf)
  d  <- vapply(seq_along(cuts), function(k)
    sum(status == 1 & time > cuts[k] & time <= ends[k]), numeric(1))
  pt <- vapply(seq_along(cuts), function(k)
    sum(pmax(0, pmin(time, ends[k]) - cuts[k])), numeric(1))
  rate <- d / pt
  loglik <- sum(ifelse(d > 0, d * log(rate), 0) - rate * pt)
  list(cuts = cuts, rate = rate, events = d, person_time = pt,
       loglik = loglik, n_par = length(cuts))
}

# Cumulative hazard of a PWE at times t.
pwe_cumhaz <- function(t, cuts, rate) {
  ends <- c(cuts[-1], Inf)
  vapply(t, function(x) sum(rate * pmax(0, pmin(x, ends) - cuts)), numeric(1))
}

pwe_surv <- function(t, cuts, rate) exp(-pwe_cumhaz(t, cuts, rate))

pwe_hazard <- function(t, cuts, rate) rate[findInterval(t, cuts)]

# Weibull fit via survreg, returned on the S(t) = exp(-(t/scale)^shape) scale.
weibull_fit <- function(time, status) {
  f <- survival::survreg(survival::Surv(time, status) ~ 1, dist = "weibull")
  shape <- 1 / f$scale
  scale <- exp(unname(coef(f)))
  # Log-likelihood on the time scale, computed directly so it is comparable to pwe_fit().
  loglik <- sum(status * log(weibull_hazard(time, shape, scale))) - sum((time / scale)^shape)
  list(shape = shape, scale = scale, loglik = loglik, n_par = 2)
}

weibull_surv <- function(t, shape, scale) exp(-(t / scale)^shape)

weibull_hazard <- function(t, shape, scale) (shape / scale) * (t / scale)^(shape - 1)
