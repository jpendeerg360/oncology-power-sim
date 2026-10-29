# Power of a time-to-event oncology trial under enrollment, dropout and delayed effects

**Question:** How do enrollment rate, dropout and a delayed treatment effect change the power of a
two-arm, event-driven overall-survival trial? Does an analysis suited to delayed effects
(Fleming–Harrington FH(0,1)) change which design to recommend, compared with log-rank?

- Report: [`report/report.md`](report/report.md) (source `report/report.Rmd`)
- Pre-specified design: [`docs/design_spec.md`](docs/design_spec.md)

## How to run
```
Rscript -e 'renv::restore()'
Rscript run_all.R
```
`run_all.R` runs scripts `01`–`05` in order and then the `testthat` suite.

Expected runtime is about 3.5 minutes on a 14-core Mac (12 workers). The 768,000-trial simulation
takes about 3 minutes of that. Stage 4 skips scenarios whose file already exists in `results/reps/`;
delete that folder to rerun from scratch. A clean rerun reproduces `results/scenario_summary.csv`
exactly.

No Quarto or pandoc is needed. The report is knitted to Markdown with `knitr::knit()`.

## Layout
```
R/                  functions: control model, simulate_trial, cut_data, analyze_trial,
                    run_scenario, summarise, recommend, design calcs, plot theme
scripts/            01_explore, 02_design_calcs, 03_pilot, 04_run_scenarios, 05_report
docs/design_spec.md pre-specified design (approved before the main run)
results/            control model, planning calcs, pilot, per-replicate results (reps/),
                    scenario summary, recommendation, run log with seeds and sessionInfo
figures/            01_* exploration figures; fig1-fig6 report figures
report/             report.Rmd -> report.md
tests/testthat/     simulator, data cut, test statistics (vs survdiff/nphRCT/simtrial), pilot OCs
```
