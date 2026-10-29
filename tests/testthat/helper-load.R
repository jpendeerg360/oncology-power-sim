# Load project functions; tests run with tests/testthat as the working directory.
project_root <- normalizePath(file.path("..", ".."))
for (f in list.files(file.path(project_root, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
control <- readRDS(file.path(project_root, "results", "control_model.rds"))
