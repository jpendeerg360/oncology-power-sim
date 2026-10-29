# Reproduce the full study: Rscript run_all.R
# Stage 4 skips scenarios whose replicate file already exists; delete results/reps/ to rerun from scratch.
scripts <- c(
  "scripts/01_explore.R",        # colon data, control-arm calibration
  "scripts/02_design_calcs.R",   # analytic planning numbers used in the spec
  "scripts/03_pilot.R",          # pilot run and runtime projection
  "scripts/04_run_scenarios.R",  # full simulation grid
  "scripts/05_report.R"          # recommendation, figures, report/report.md
)
for (s in scripts) {
  message("== Running ", s)
  source(s, local = new.env())
}
message("== Running tests")
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
