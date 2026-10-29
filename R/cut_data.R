# Event-driven data cut (spec section 3): analyze at the D-th observed death or at the cap,
# whichever comes first. Patients not yet enrolled at the cut are excluded.

cut_data <- function(trial, D, cap) {
  observed_death <- trial$death_time < trial$dropout_time
  calendar_death <- ifelse(observed_death, trial$entry + trial$death_time, Inf)

  if (sum(calendar_death <= cap) >= D) {
    t_star <- sort(calendar_death, partial = D)[D]
    trigger <- "events"
  } else {
    t_star <- cap
    trigger <- "cap"
  }

  keep <- trial$entry <= t_star
  enrolled <- trial[keep, ]
  # Death status is decided on the calendar scale: recomputing t_star - entry in floating point
  # can drop the D-th death, which sits exactly at t_star.
  status <- as.integer(calendar_death[keep] <= t_star)
  censor_time <- pmin(enrolled$dropout_time, t_star - enrolled$entry)
  data <- data.frame(
    time = ifelse(status == 1, enrolled$death_time, censor_time),
    status = status,
    trt = enrolled$trt
  )

  list(data = data, t_star = t_star, trigger = trigger,
       n_enrolled = nrow(data), deaths = sum(data$status))
}
