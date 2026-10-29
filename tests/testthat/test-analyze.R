make_datasets <- function() {
  sets <- list()
  for (s in 1:6) {
    set.seed(200 + s)
    pat <- list(c(1, 0), c(0.69, 0), c(0.69, 12))[[1 + s %% 3]]
    tr <- simulate_trial(c(300, 600)[1 + s %% 2], 30, 0.05, pat[1], pat[2], control)
    sets[[s]] <- cut_data(tr, 150, 84)$data
  }
  # One dataset with tied times (rounded to whole months).
  tied <- sets[[1]]; tied$time <- pmax(1, round(tied$time))
  c(sets, list(tied))
}

test_that("log-rank Z matches survival::survdiff, including sign", {
  for (d in make_datasets()) {
    z <- weighted_logrank(d$time, d$status, d$trt)[["z_lr"]]
    sd <- survival::survdiff(survival::Surv(time, status) ~ trt, data = d)
    expect_equal(z^2, unname(sd$chisq), tolerance = 1e-8)
    # survdiff row 2 is trt == 1; O - E < 0 favours treatment, as does z < 0.
    expect_equal(sign(z), sign(sd$obs[2] - sd$exp[2]))
  }
})

test_that("FH(0,1) Z matches nphRCT::wlrt", {
  for (d in make_datasets()) {
    ours <- weighted_logrank(d$time, d$status, d$trt)[["z_fh"]]
    ref <- nphRCT::wlrt(survival::Surv(time, status) ~ trt, data = d,
                        method = "fh", rho = 0, gamma = 1, timefix = FALSE)
    # nphRCT reports U for the arm named in trt_group; align the sign to trt == 1.
    ref_z <- if (as.character(ref$trt_group) == "1") ref$z else -ref$z
    expect_equal(ours, unname(ref_z), tolerance = 1e-8)
  }
})

test_that("FH(0,1) Z agrees with simtrial::wlr", {
  for (d in make_datasets()) {
    ours <- weighted_logrank(d$time, d$status, d$trt)[["z_fh"]]
    d$arm <- ifelse(d$trt == 1, "experimental", "control")
    ref <- simtrial::wlr(d, weight = simtrial::fh(rho = 0, gamma = 1),
                         formula = survival::Surv(time, status) ~ arm)
    # simtrial's z is oriented so that positive favours experimental.
    expect_equal(abs(ours), abs(ref$z), tolerance = 1e-6)
  }
})

test_that("analyze_trial returns one-sided p-values consistent with Z", {
  d <- make_datasets()[[2]]
  res <- analyze_trial(d)
  expect_equal(res[["p_lr"]], pnorm(res[["z_lr"]]))
  expect_equal(res[["p_fh"]], pnorm(res[["z_fh"]]))
  cox <- survival::coxph(survival::Surv(time, status) ~ trt, data = d)
  expect_equal(res[["loghr"]], unname(coef(cox)))
})
