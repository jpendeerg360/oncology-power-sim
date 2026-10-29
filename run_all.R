# Reproduce the full study: Rscript run_all.R
# Stages are added as they are approved.
scripts <- c(
  "scripts/01_explore.R",
  "scripts/02_design_calcs.R",
  "scripts/03_pilot.R"
)
for (s in scripts) {
  message("== Running ", s)
  source(s, local = new.env())
}
