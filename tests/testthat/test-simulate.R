test_that("simulated control-arm KM matches the calibrated PWE curve", {
  set.seed(101)
  tr <- simulate_trial(N = 40000, enrol_rate = 1000, dropout_annual = 0, hr = 1, delay = 0,
                       control = control)
  ctrl <- tr[tr$trt == 0, ]
  km <- survival::survfit(survival::Surv(death_time) ~ 1, data = ctrl)
  t_check <- c(3, 6, 9, 12, 18, 24, 36, 48, 60, 84)
  km_s <- summary(km, times = t_check)$surv
  # SE of S(t) at n = 20000 is at most ~0.0035; allow ~4 SE.
  expect_lt(max(abs(km_s - pwe_surv(t_check, control$cuts, control$rate))), 0.015)
})

test_that("treatment arm follows the delayed-effect hazard", {
  set.seed(102)
  tr <- simulate_trial(N = 40000, enrol_rate = 1000, dropout_annual = 0, hr = 0.69, delay = 12,
                       control = control)
  treated <- tr[tr$trt == 1, ]
  km <- survival::survfit(survival::Surv(death_time) ~ 1, data = treated)
  th <- treatment_hazard(control$cuts, control$rate, hr = 0.69, delay = 12)
  t_check <- c(6, 12, 18, 24, 36, 48, 60, 84)
  expect_lt(max(abs(summary(km, times = t_check)$surv - pwe_surv(t_check, th$cuts, th$rate))), 0.015)
  # Before the delay, the treatment arm is identical in distribution to control.
  expect_equal(pwe_surv(12, th$cuts, th$rate), pwe_surv(12, control$cuts, control$rate))
})

test_that("pwe_inverse inverts the cumulative hazard", {
  t <- c(0.5, 5.9, 6, 13, 30, 47.9, 50, 120)
  E <- pwe_cumhaz(t, control$cuts, control$rate)
  expect_equal(pwe_inverse(E, control$cuts, control$rate), t)
})

test_that("allocation is exactly 1:1 and dropout follows the annual rate", {
  set.seed(103)
  tr <- simulate_trial(N = 20000, enrol_rate = 30, dropout_annual = 0.15, hr = 1, delay = 0,
                       control = control)
  expect_equal(sum(tr$trt), 10000)
  expect_false(is.unsorted(tr$entry))
  expect_equal(mean(tr$dropout_time <= 12), 0.15, tolerance = 0.01)
  expect_equal(mean(diff(tr$entry)), 1 / 30, tolerance = 0.03)
})
