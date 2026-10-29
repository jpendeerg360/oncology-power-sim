# Power of a time-to-event oncology trial under enrollment, dropout and delayed effects

**Question:** How do enrollment rate, dropout and a delayed treatment effect change the power of a
two-arm time-to-event oncology trial? Does an analysis suited to delayed effects
(Fleming–Harrington FH(0,1)) change which design to recommend?

## How to run
```
Rscript -e 'renv::restore()'
Rscript run_all.R
```

## Status
- [x] Stage 1: explore `survival::colon`, calibrate the control-arm model (`scripts/01_explore.R`)
- [ ] Stage 2: design spec (`docs/design_spec.md`)
- [ ] Stage 3: simulator and validation
- [ ] Stage 4: full simulation
- [ ] Stage 5: report

Expected runtime: TBD after Stage 3.
