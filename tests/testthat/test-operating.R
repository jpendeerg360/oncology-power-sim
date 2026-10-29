# Operating-characteristic checks on the pilot (scripts/03_pilot.R must have been run).
pilot_file <- file.path(project_root, "results", "pilot_summary.rds")

test_that("log-rank power under PH without a binding cap is close to Schoenfeld", {
  skip_if_not(file.exists(pilot_file), "pilot not run")
  s <- readRDS(pilot_file)$summary
  no_cap <- s[s$pattern == "b_ph" & s$pct_cap == 0, ]
  expect_gt(nrow(no_cap), 0)
  # Within 3 Monte Carlo SEs of the Schoenfeld prediction at the deaths analysed.
  expect_true(all(abs(no_cap$power_lr - no_cap$schoenfeld_power) < 3 * no_cap$mcse_power_lr))
})

test_that("type I error is about 0.025 for both tests", {
  skip_if_not(file.exists(pilot_file), "pilot not run")
  s <- readRDS(pilot_file)$summary
  null <- s[s$pattern == "a_null", ]
  expect_true(all(abs(null$power_lr - 0.025) < 3 * null$mcse_power_lr))
  expect_true(all(abs(null$power_fh - 0.025) < 3 * null$mcse_power_fh))
})

test_that("Cox CI coverage under PH is about 95%", {
  skip_if_not(file.exists(pilot_file), "pilot not run")
  s <- readRDS(pilot_file)$summary
  ph <- s[s$pattern == "b_ph", ]
  expect_true(all(abs(ph$coverage - 0.95) < 3 * ph$mcse_coverage))
})

test_that("simulated control-arm KM matches the calibrated curve in the pilot", {
  skip_if_not(file.exists(pilot_file), "pilot not run")
  km <- readRDS(pilot_file)$km_check
  expect_lt(max(abs(km$km - km$pwe)), 0.01)
})
