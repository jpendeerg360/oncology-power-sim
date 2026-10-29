# Reproduce the full study: Rscript run_all.R
# Stages are added as they are approved.
scripts <- c(
  "scripts/01_explore.R"
)
for (s in scripts) {
  message("== Running ", s)
  source(s, local = new.env())
}
