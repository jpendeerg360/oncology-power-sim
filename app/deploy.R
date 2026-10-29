# Deploy the design calculator to shinyapps.io: Rscript app/deploy.R
# Requires a one-time rsconnect::setAccountInfo(name, token, secret) on this machine.
# The bundle is assembled in a temporary directory so the app ships with the same R/ files
# and control model that the tests use, without copying them into the repository.

bundle <- file.path(tempdir(), "design-calculator")
unlink(bundle, recursive = TRUE)
dir.create(file.path(bundle, "results"), recursive = TRUE)
file.copy("app/app.R", bundle)
file.copy("R", bundle, recursive = TRUE)
file.copy("results/control_model.rds", file.path(bundle, "results"))

rsconnect::deployApp(bundle, appName = "oncology-design-calculator", forceUpdate = TRUE)
