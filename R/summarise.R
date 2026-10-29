# Operating characteristics for one scenario from its per-replicate results (spec section 6).

ALPHA <- 0.025

mcse_prop <- function(p, R) sqrt(p * (1 - p) / R)

summarise_scenario <- function(reps, pattern, hr) {
  R <- nrow(reps)
  rej_lr <- reps$p_lr < ALPHA
  rej_fh <- reps$p_fh < ALPHA
  p10 <- mean(rej_fh & !rej_lr)
  p01 <- mean(!rej_fh & rej_lr)
  ci_width <- 2 * qnorm(0.975) * reps$se_loghr
  is_ph <- pattern == "b_ph"
  covered <- abs(reps$loghr - log(hr)) <= qnorm(0.975) * reps$se_loghr

  data.frame(
    R = R,
    power_lr = mean(rej_lr), mcse_power_lr = mcse_prop(mean(rej_lr), R),
    power_fh = mean(rej_fh), mcse_power_fh = mcse_prop(mean(rej_fh), R),
    diff_fh_lr = p10 - p01, mcse_diff = sqrt((p10 + p01 - (p10 - p01)^2) / R),
    loghr_mean = mean(reps$loghr), loghr_sd = sd(reps$loghr),
    mcse_loghr_mean = sd(reps$loghr) / sqrt(R),
    ci_width_mean = mean(ci_width), mcse_ci_width = sd(ci_width) / sqrt(R),
    bias = if (is_ph) mean(reps$loghr) - log(hr) else NA_real_,
    coverage = if (is_ph) mean(covered) else NA_real_,
    mcse_coverage = if (is_ph) mcse_prop(mean(covered), R) else NA_real_,
    duration_mean = mean(reps$t_star), mcse_duration = sd(reps$t_star) / sqrt(R),
    pct_cap = mean(reps$trigger == "cap"), mcse_pct_cap = mcse_prop(mean(reps$trigger == "cap"), R),
    deaths_mean = mean(reps$deaths), enrolled_mean = mean(reps$n_enrolled)
  )
}
