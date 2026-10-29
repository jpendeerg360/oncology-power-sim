# Checks on the exploratory design calculator (R/design_tool.R).
summary_file <- file.path(project_root, "results", "scenario_summary.csv")
D <- 306; CAP <- 84

test_that("simulate_design reproduces the base-case grid scenario exactly", {
  skip_if_not(file.exists(summary_file), "full grid not run")
  grid_row <- read.csv(summary_file)
  grid_row <- grid_row[grid_row$scenario_id == 120, ]   # d_delay12, N = 900, 30/month, 5%/yr
  future::plan(future::sequential)
  s <- simulate_design(900, 30, 0.05, 0.69, 12, D, CAP, control,
                       reps = grid_row$reps, seed = grid_row$seed)
  expect_equal(s$power_lr, grid_row$power_lr)
  expect_equal(s$power_fh, grid_row$power_fh)
  expect_equal(s$duration_mean, grid_row$duration_mean)
})

test_that("simulate_design is reproducible for a fixed seed", {
  future::plan(future::sequential)
  a <- simulate_design(400, 30, 0.05, 0.69, 6, D, CAP, control, reps = 50, seed = 1)
  b <- simulate_design(400, 30, 0.05, 0.69, 6, D, CAP, control, reps = 50, seed = 1)
  expect_identical(a, b)
  expect_equal(a$mcse_power_lr, mcse_prop(a$power_lr, 50))
})

test_that("gsDesign2 inputs reproduce the study's hazards and accrual", {
  r <- gs_rates(900, 30, 0.05, 0.69, 12, control)
  trt <- treatment_hazard(control$cuts, control$rate, 0.69, 12)
  expect_equal(r$fail$fail_rate * r$fail$hr, trt$rate)
  expect_equal(r$fail$dropout_rate[1], dropout_rate(0.05))
  expect_equal(sum(r$enroll$duration * r$enroll$rate), 900)
})

test_that("analytic time to D deaths agrees with design_calc.R", {
  a <- analytic_design(900, 30, 0.05, 0.69, 12, D, CAP, control)
  t_dc <- time_to_events(D, 900, 30, control$cuts, control$rate, 0.69, 12, dropout_rate(0.05))
  expect_lt(abs(a$time - t_dc), 0.5)
  expect_equal(a$events, D, tolerance = 1e-3)
})

test_that("analytic power is alpha under the null and Schoenfeld-like under PH", {
  null <- analytic_design(900, 30, 0.05, 1, 0, D, CAP, control)
  expect_equal(null$power_lr_ahr, ALPHA, tolerance = 1e-3)
  ph <- analytic_design(900, 30, 0.05, 0.69, 0, D, CAP, control)
  expect_lt(abs(ph$power_lr_ahr - ph$power_schoenfeld), 0.01)
})

test_that("analytic power under the base case is within 3 points of the simulated grid", {
  skip_if_not(file.exists(summary_file), "full grid not run")
  s <- read.csv(summary_file)
  s <- s[s$scenario_id %in% c(118, 119, 120), ]
  a <- do.call(rbind, lapply(seq_len(nrow(s)), function(i)
    analytic_design(s$N[i], s$enrol[i], s$dropout[i], s$hr[i], s$delay[i], D, CAP, control)))
  expect_true(all(abs(a$power_lr_ahr - s$power_lr) < 0.03))
  expect_true(all(abs(a$power_fh_ahr - s$power_fh) < 0.03))
})

test_that("analysis time falls back to the cap when D deaths are not expected by then", {
  a <- analytic_design(400, 15, 0.15, 0.69, 12, D, CAP, control)
  expect_equal(a$time, CAP)
  expect_lt(a$events, D)
})
