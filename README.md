# Power of a time-to-event oncology trial under enrollment, dropout and delayed effects

A pre-specified simulation study of how **enrollment rate**, **dropout** and a **delayed treatment
effect** change the power of a two-arm, event-driven overall-survival (OS) trial, and whether an
analysis built for delayed effects (the Fleming–Harrington FH(0,1) weighted log-rank test) changes
which design should be recommended compared with the standard log-rank test.

- **Full report:** [`report/report.md`](report/report.md)
- **Pre-specified design:** [`docs/design_spec.md`](docs/design_spec.md), written and approved before any scenario was run

---

## Key findings

- **Under proportional hazards (PH), log-rank is clearly the better test.** At N = 900, log-rank
  reaches 90.3% power, matching Schoenfeld's planning value. FH(0,1) reaches only 80.1%, about 10
  points lower, because it down-weights early deaths that carry real information under PH.
- **A 6-month delay costs little in this population.** Few patients die in the first six months, so
  log-rank still reaches 85.7% at N = 900 and still beats FH(0,1).
- **A 12-month delay changes the picture.** No design reaches 80% power. FH(0,1) now beats log-rank
  (by 3.0 points at N = 900 and 6.2 points at N = 1200), but the best achievable power is about 76%.
- **More patients can mean less power.** With an event-driven analysis, a larger N or faster
  enrollment brings the analysis forward. More of the 306 deaths then fall inside the delay period,
  where HR = 1, adding noise but no signal. Log-rank power falls from 73.8% at N = 800 to 63.3% at
  N = 1200.
- **Both tests control type I error** at the nominal one-sided 2.5% across all 48 null scenarios
  (log-rank 2.21–2.86%, FH(0,1) 2.02–2.77%; the two highest values are consistent with Monte Carlo
  noise).

**Recommendation.** No design meets the pre-specified target (≥ 80% power, ≤ 72 months expected
duration) in the base case. The closest design is **N = 900, analyzed with FH(0,1) at 306 deaths or
84 months**, with 74.7% power and 62.7 months expected duration. The choice of test hinges on how
confident one is in a ~1-year delay: if the effect is closer to PH, log-rank is the right choice.

<p align="center">
  <img src="figures/fig3_power_vs_n.png" alt="Power versus sample size under PH and a 12-month delay" width="85%">
</p>

---

## Background

Event-driven trials are sized so that a fixed number of deaths gives the desired power under
proportional hazards. Immunotherapies often show a *delayed* separation of survival curves, which
breaks that assumption: deaths that occur before the effect begins dilute the comparison. How badly
this hurts depends on trial conduct (how fast patients enroll, how many drop out) and on the choice
of test. This study quantifies those trade-offs with a realistic control arm.

## Data grounding

The control arm is calibrated to the observation-only arm of `survival::colon`, a 1980s adjuvant
chemotherapy trial in resected colon cancer (929 patients and 452 deaths overall, of whom 315 patients
and 168 deaths are in the observation arm; median follow-up 77 months).

| | |
|:--|:--|
| **Control hazard** | Five-interval piecewise-exponential model (breaks at 6, 12, 24, 48 months). The hazard rises to a peak in months 12–24 and then falls, a shape a Weibull cannot fit (AIC 1854.0 vs 1881.1; max distance from Kaplan–Meier 0.020 vs 0.082). |
| **Effect size** | Cox HR for Lev+5FU vs Obs is 0.689 (95% CI 0.546–0.869), with no evidence against PH (`cox.zph` p = 0.28). The study uses HR = 0.69. |
| **Not grounded by the data** | The delayed-effect patterns are hypothetical, motivated by immunotherapy trials. The data also say nothing about modern enrollment or dropout rates. |

## Study design

| Element | Specification |
|:--|:--|
| **Treatment effect patterns** | (a) null; (b) PH with HR 0.69; (c) HR 1 for 6 months, then 0.69; (d) HR 1 for 12 months, then 0.69 |
| **Sample size N** | 400, 500, 600, 700, 800, 900, 1000, 1200 (1:1 randomization) |
| **Enrollment** | Poisson arrivals at 15 or 30 patients/month |
| **Dropout** | Exponential, 0%, 5% or 15% per year, independent of arm |
| **Analysis timing** | At D = 306 deaths (Schoenfeld: HR 0.69, 90% power, one-sided α = 0.025) or 84 months after activation, whichever comes first |
| **Tests** | Log-rank and FH(0,1), applied to the same simulated datasets (paired comparison); success is one-sided p < 0.025 |
| **Replicates** | 2,000 per non-null scenario, 10,000 per null scenario (192 scenarios, 768,000 simulated trials) |
| **Precision** | Monte Carlo SE ≤ 0.011 for power, 0.0016 for type I error |
| **Base case** | 12-month delay, 30 patients/month, 5%/yr dropout |

**Decision rule (pre-specified).** A design (N, test) is eligible if its type I error is ≤ 2.81%
(0.025 + 2 MCSE) in every null scenario with that N. Among eligible designs, choose the smallest N with
≥ 80% power and ≤ 72 months expected duration in the base case. Use log-rank unless FH(0,1) is better
by more than 2 paired MCSEs. If nothing qualifies, report the closest design.

## Methods

- **Simulation.** Death times are drawn exactly by inverting the piecewise-linear cumulative hazard.
  Enrollment, dropout and death are simulated on calendar time, and each trial is cut at the analysis
  trigger.
- **Test statistics.** Log-rank and FH(0,1) (weight 1 − Ŝ(t−) from the pooled Kaplan–Meier) are coded
  directly and validated against `survival::survdiff`, `nphRCT::wlrt` and `simtrial::wlr`.
- **Estimation.** Each trial also records the Cox HR, its 95% CI, the number of deaths and the
  analysis time, so bias, coverage and expected duration are reported alongside power.
- **Reproducibility.** Every scenario has its own seed derived from a master seed, the package
  environment is pinned with `renv`, and the run log records seeds and `sessionInfo()`. A clean rerun
  reproduces `results/scenario_summary.csv` exactly.

## Limitations

- The delay patterns are hypothetical clean steps; real delayed effects are often gradual.
- Uncertainty in the fitted control hazard is not propagated into the simulation.
- Enrollment is at a constant rate with no ramp-up, and dropout is non-informative.
- The design has a single analysis with no interim looks.

Natural extensions include group-sequential designs, MaxCombo tests that hedge between log-rank and
FH weights, restricted mean survival time (RMST), and Bayesian designs that use prior information
about the delay.

---

## Reproducing the results

Requires R (developed on R 4.5.1). Package versions are pinned in `renv.lock`.

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R
```

`run_all.R` runs the five pipeline stages in order and then the test suite. A full run takes about
3.5 minutes on a 14-core machine using 12 parallel workers. Stage 4 skips scenarios whose results are
already saved in `results/reps/`; delete that folder to rerun the simulation from scratch.

| Stage | Script | What it does |
|:--|:--|:--|
| 1 | `scripts/01_explore.R` | Explores the colon OS data and calibrates the control-arm model |
| 2 | `scripts/02_design_calcs.R` | Analytic planning calculations used in the design spec |
| 3 | `scripts/03_pilot.R` | Pilot run, operating-characteristic checks and runtime projection |
| 4 | `scripts/04_run_scenarios.R` | Runs the full 192-scenario grid in parallel |
| 5 | `scripts/05_report.R` | Applies the decision rule, draws figures and knits the report |

## Design calculator (exploratory)

`app/` is a Shiny app for checking a single design beyond the pre-specified grid. You set N,
enrollment rate, dropout, HR, delay, target deaths and the calendar cap. It reports simulated
log-rank and FH(0,1) power with Monte Carlo SEs, expected analysis time and the share of analyses
triggered by the cap, next to analytic approximations from `gsDesign2` (average-HR log-rank,
weighted FH(0,1)) and Schoenfeld.

It is **not** part of the pre-specified study. The study's results and recommendation above come from
the grid and decision rule in `docs/design_spec.md`. The app calls the same tested functions
(`simulate_trial`, `cut_data`, `analyze_trial`, `run_scenario`) through `R/design_tool.R`, and
`tests/testthat/test-design-tool.R` checks that it reproduces the base-case grid scenario exactly.

```sh
Rscript -e 'shiny::runApp("app")'   # run locally from the project root
Rscript app/deploy.R                # deploy to shinyapps.io (after rsconnect::setAccountInfo)
```

## Repository layout

| Path | Contents |
|:--|:--|
| `R/` | Core functions: control model, `simulate_trial`, `cut_data`, `analyze_trial`, `run_scenario`, `summarise`, `recommend`, design calculations, plot theme |
| `scripts/` | The five pipeline stages listed above |
| `docs/design_spec.md` | Pre-specified study design |
| `results/` | Control model, planning calculations, pilot results, per-scenario replicates (`reps/`), scenario summary, recommendation, run log |
| `figures/` | `01_*` exploration figures; `fig1`–`fig6` report figures |
| `report/` | `report.Rmd` and the rendered `report.md` |
| `app/` | Exploratory design calculator (Shiny) and its deploy script |
| `tests/testthat/` | Tests for the simulator, data cut, test statistics (against `survdiff`, `nphRCT`, `simtrial`), pilot operating characteristics and the design calculator |
