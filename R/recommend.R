# Pre-specified recommendation rule and one-factor-at-a-time sensitivity (spec section 8).

POWER_TARGET <- 0.80
METHODS <- c(lr = "Log-rank", fh = "FH(0,1)")

# One row per scenario x method.
to_long <- function(s) {
  rbind(
    transform(s, method = "lr", power = power_lr, mcse_power = mcse_power_lr,
              other_power = power_fh, diff_other = diff_fh_lr),
    transform(s, method = "fh", power = power_fh, mcse_power = mcse_power_fh,
              other_power = power_lr, diff_other = -diff_fh_lr)
  )
}

# Max type I error over the null scenarios sharing each design's N.
design_eligibility <- function(s, threshold) {
  null <- to_long(s[s$pattern == "a_null", ])
  elig <- aggregate(power ~ N + method, data = null, FUN = max)
  names(elig)[names(elig) == "power"] <- "max_type1"
  elig$eligible <- elig$max_type1 <= threshold
  elig
}

# Apply steps 1-5 of the rule in the scenario identified by `base`.
apply_rule <- function(s, base, threshold, dur_max) {
  cell <- to_long(s[s$pattern == base$pattern & s$enrol == base$enrol & s$dropout == base$dropout, ])
  cell <- merge(cell, design_eligibility(s, threshold), by = c("N", "method"))
  cell$meets_power <- cell$power >= POWER_TARGET
  cell$meets_duration <- cell$duration_mean <= dur_max
  cell$qualifies <- cell$eligible & cell$meets_power & cell$meets_duration

  if (any(cell$qualifies)) {
    status <- "qualifies"
    n_star <- min(cell$N[cell$qualifies])
    cand <- cell[cell$N == n_star & cell$qualifies, ]
  } else {
    # Deviation (documented in the report): "closest" = highest power among eligible designs
    # within the duration limit; if none is within it, highest power among all eligible designs.
    # Smaller N breaks ties.
    pool <- cell[cell$eligible & cell$meets_duration, ]
    status <- "closest"
    if (nrow(pool) == 0) {
      pool <- cell[cell$eligible, ]
      status <- "closest_over_duration"
    }
    best <- max(pool$power)
    n_star <- min(pool$N[pool$power == best])
    cand <- pool[pool$N == n_star, ]
  }
  # Step 4: log-rank unless FH(0,1) is better by more than 2 paired MCSEs.
  method <- if (nrow(cand) == 1) cand$method else {
    if (cand$diff_fh_lr[1] > 2 * cand$mcse_diff[1]) "fh" else "lr"
  }
  list(status = status, N = n_star, method = method,
       row = cell[cell$N == n_star & cell$method == method, ], table = cell)
}

# Sensitivity of the recommended design to one factor at a time.
sensitivity_table <- function(s, rec, base, levels, threshold, dur_max) {
  do.call(rbind, lapply(seq_len(nrow(levels)), function(i) {
    lev <- base
    lev[[levels$factor[i]]] <- levels$value[[i]]
    long <- to_long(s[s$pattern == lev$pattern & s$enrol == lev$enrol & s$dropout == lev$dropout, ])
    r <- long[long$N == rec$N & long$method == rec$method, ]
    alt <- apply_rule(s, lev, threshold, dur_max)
    data.frame(
      factor = levels$factor[i], level = as.character(levels$value[[i]]),
      power = r$power, mcse_power = r$mcse_power,
      other_power = r$other_power, diff_other = r$diff_other, mcse_diff = r$mcse_diff,
      duration = r$duration_mean, pct_cap = r$pct_cap,
      below_target = r$power < POWER_TARGET,
      other_wins = r$diff_other > 2 * r$mcse_diff,
      over_duration = r$duration_mean > dur_max,
      rule_status = alt$status, rule_N = alt$N, rule_method = METHODS[[alt$method]],
      rule_power = alt$row$power
    )
  })) |> transform(sensitive = below_target | other_wins | over_duration)
}
