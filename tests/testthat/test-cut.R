test_that("event-triggered cuts contain exactly D deaths", {
  D <- 306; cap <- 84
  n_events <- 0
  for (s in 1:200) {
    set.seed(s)
    tr <- simulate_trial(900, 30, 0.05, 0.69, 12, control)
    ct <- cut_data(tr, D, cap)
    if (ct$trigger == "events") {
      n_events <- n_events + 1
      expect_equal(ct$deaths, D)
      expect_lte(ct$t_star, cap)
    } else {
      expect_equal(ct$t_star, cap)
      expect_lt(ct$deaths, D)
    }
    expect_true(all(ct$data$time >= 0))
  }
  expect_gt(n_events, 150)
})

test_that("cap-triggered cuts stop at the cap and exclude later enrollees", {
  set.seed(7)
  # 10/month accrues 1200 patients over ~120 months, so enrollment is still open at the cap.
  tr <- simulate_trial(1200, 10, 0.15, 0.69, 0, control)
  ct <- cut_data(tr, 306, 84)
  expect_equal(ct$trigger, "cap")
  expect_equal(ct$t_star, 84)
  expect_equal(ct$n_enrolled, sum(tr$entry <= 84))
  expect_lt(ct$n_enrolled, 1200)
})

test_that("cut data agree with a direct calendar-time calculation", {
  set.seed(11)
  tr <- simulate_trial(600, 30, 0.05, 0.69, 6, control)
  ct <- cut_data(tr, 150, 84)
  cal <- with(tr, ifelse(death_time < dropout_time, entry + death_time, Inf))
  expect_equal(ct$t_star, sort(cal)[150])
  expect_equal(sum(ct$data$status), sum(cal <= ct$t_star))
})
