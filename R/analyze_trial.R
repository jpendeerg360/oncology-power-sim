# Analyses applied to the same cut dataset (spec section 5).

# Log-rank and FH(0,1) Z statistics in one pass. Negative Z favours treatment (trt == 1).
# FH(0,1) weight is 1 - S(t-), with S the pooled Kaplan-Meier estimate.
weighted_logrank <- function(time, status, trt) {
  o <- order(time)
  time <- time[o]; status <- status[o]; trt <- trt[o]

  death_times <- unique(time[status == 1])
  n_risk <- length(time) - findInterval(death_times, time, left.open = TRUE)
  n1_risk <- sum(trt) - findInterval(death_times, time[trt == 1], left.open = TRUE)
  d <- tabulate(match(time[status == 1], death_times), length(death_times))
  d1 <- tabulate(match(time[status == 1 & trt == 1], death_times), length(death_times))

  p1 <- n1_risk / n_risk
  o_minus_e <- d1 - d * p1
  v <- d * p1 * (1 - p1) * (n_risk - d) / pmax(n_risk - 1, 1)

  s_minus <- c(1, cumprod(1 - d / n_risk))[seq_along(death_times)]
  w_fh <- 1 - s_minus

  c(z_lr = sum(o_minus_e) / sqrt(sum(v)),
    z_fh = sum(w_fh * o_minus_e) / sqrt(sum(w_fh^2 * v)))
}

analyze_trial <- function(data) {
  z <- weighted_logrank(data$time, data$status, data$trt)
  cox <- survival::coxph(survival::Surv(time, status) ~ trt, data = data)
  c(z_lr = z[["z_lr"]], p_lr = pnorm(z[["z_lr"]]),
    z_fh = z[["z_fh"]], p_fh = pnorm(z[["z_fh"]]),
    loghr = unname(coef(cox)), se_loghr = sqrt(vcov(cox)[1, 1]))
}
