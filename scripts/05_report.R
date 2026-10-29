# Stage 5: apply the recommendation rule, draw the report figures and knit the report.
# Outputs: results/recommendation.rds, figures/fig1-6_*.png, report/report.md

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

s <- read.csv("results/scenario_summary.csv")
dc <- readRDS("results/design_calcs.rds")
control <- readRDS("results/control_model.rds")

# ---- Recommendation and sensitivity (spec section 8) --------------------------
BASE <- list(pattern = "d_delay12", enrol = 30, dropout = 0.05)
rec <- apply_rule(s, BASE, dc$type1_threshold, dc$dur_max)
sens_levels <- data.frame(factor = c("enrol", "dropout", "dropout", "pattern", "pattern"))
sens_levels$value <- list(15, 0, 0.15, "b_ph", "c_delay6")
sens <- sensitivity_table(s, rec, BASE, sens_levels, dc$type1_threshold, dc$dur_max)
elig <- design_eligibility(s, dc$type1_threshold)

null_long <- to_long(s[s$pattern == "a_null", ])
type1_exceed <- null_long[null_long$power > dc$type1_threshold, ]
# Exceedances expected by chance if every true type I error were exactly 0.025.
type1_expected_exceed <- nrow(null_long) *
  (1 - pnorm((dc$type1_threshold - 0.025) / mcse_prop(0.025, 10000)))

# Exploratory (not pre-specified): deaths needed for 80% log-rank power under the base case,
# using Schoenfeld with the time-averaged Cox log HR at the recommended N.
base_rec_row <- s[s$pattern == BASE$pattern & s$enrol == BASE$enrol &
                    s$dropout == BASE$dropout & s$N == rec$N, ]
d_needed_explore <- ceiling(4 * (qnorm(0.975) + qnorm(0.8))^2 / base_rec_row$loghr_mean^2)

saveRDS(list(base = BASE, rec = rec, sens = sens, eligibility = elig,
             type1_exceed = type1_exceed, type1_expected_exceed = type1_expected_exceed,
             d_needed_explore = d_needed_explore),
        "results/recommendation.rds")

# ---- Figure helpers -----------------------------------------------------------
METHOD_COLOURS <- c("Log-rank" = SERIES_COLOURS[[1]], "FH(0,1)" = SERIES_COLOURS[[2]])
PATTERN_LABELS <- c(a_null = "(a) Null", b_ph = "(b) PH", c_delay6 = "(c) 6-month delay",
                    d_delay12 = "(d) 12-month delay")
PATTERN_COLOURS <- setNames(c("grey55", SERIES_COLOURS[1:3]), PATTERN_LABELS)
enrol_lab <- function(x) paste0(x, " patients/month")
dropout_lab <- function(x) paste0(100 * x, "%/yr dropout")
save_fig <- function(p, name, w = 7.5, h = 4.5) ggsave(file.path("figures", name), p, width = w, height = h, dpi = 150)

long <- to_long(s) |>
  mutate(method_label = factor(METHODS[method], levels = METHODS),
         pattern_label = factor(PATTERN_LABELS[pattern], levels = PATTERN_LABELS))

# ---- Fig 1: colon KM with calibrated control curve ----------------------------
os <- subset(colon, etype == 2)
os$months <- os$time / DAYS_PER_MONTH
km <- survfit(Surv(months, status) ~ rx, data = os)
km_df <- data.frame(time = c(0, 0, 0, km$time), surv = c(1, 1, 1, km$surv),
                    arm = c(levels(os$rx), rep(sub(".*=", "", names(km$strata)), km$strata))) |>
  mutate(arm = factor(arm, levels = levels(os$rx)))
t_grid <- seq(0, max(os$months), length.out = 300)
pwe_df <- data.frame(time = t_grid, surv = pwe_surv(t_grid, control$cuts, control$rate))
fig1 <- ggplot() +
  geom_step(data = km_df, aes(time, surv, colour = arm), linewidth = 0.7) +
  geom_line(data = pwe_df, aes(time, surv, linetype = "Calibrated control (piecewise exponential)"),
            colour = "black", linewidth = 0.7) +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  scale_linetype_manual(values = "22", name = NULL) +
  scale_x_continuous(breaks = seq(0, 108, 12)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Months since randomization", y = "Overall survival",
       title = "Figure 1. Overall survival in survival::colon with the calibrated control model") +
  theme_sim() + theme(legend.position = "top", legend.box = "vertical")
save_fig(fig1, "fig1_colon_km_calibration.png")

# ---- Fig 2: HR(t) patterns and implied survival -------------------------------
hr_t <- seq(0, 84, by = 0.25)
hr_df <- bind_rows(lapply(seq_len(nrow(PATTERNS)), function(i) {
  p <- PATTERNS[i, ]
  data.frame(pattern = PATTERN_LABELS[[p$pattern]], t = hr_t,
             hr = ifelse(p$hr_on & hr_t >= p$delay, dc$hr_ph, 1))
})) |> mutate(pattern = factor(pattern, levels = PATTERN_LABELS))
fig2a <- ggplot(hr_df, aes(t, hr, colour = pattern)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~pattern, nrow = 1) +
  scale_colour_manual(values = PATTERN_COLOURS, guide = "none") +
  scale_x_continuous(breaks = seq(0, 84, 24)) +
  scale_y_continuous(limits = c(0.6, 1.05)) +
  labs(x = "Months since randomization", y = "HR(t)", title = "A. Hazard ratio over time") +
  theme_sim()
surv_df <- bind_rows(
  data.frame(curve = "Control", t = hr_t, surv = pwe_surv(hr_t, control$cuts, control$rate)),
  bind_rows(lapply(2:4, function(i) {
    p <- PATTERNS[i, ]
    th <- treatment_hazard(control$cuts, control$rate, dc$hr_ph, p$delay)
    data.frame(curve = PATTERN_LABELS[[p$pattern]], t = hr_t, surv = pwe_surv(hr_t, th$cuts, th$rate))
  }))
) |> mutate(curve = factor(curve, levels = c("Control", PATTERN_LABELS[-1])))
fig2b <- ggplot(surv_df, aes(t, surv, colour = curve)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(Control = "black", PATTERN_COLOURS[-1]), name = NULL) +
  scale_x_continuous(breaks = seq(0, 84, 12)) +
  scale_y_continuous(limits = c(0.4, 1), labels = scales::percent) +
  labs(x = "Months since randomization", y = "Overall survival",
       title = "B. Implied survival, treatment arm vs control") +
  theme_sim() + theme(legend.position = "right")
fig2 <- (fig2a / fig2b) + plot_layout(heights = c(1, 1.4)) +
  plot_annotation(title = "Figure 2. Simulated treatment-effect patterns",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(fig2, "fig2_hr_patterns.png", h = 6.5)

# ---- Fig 3: power vs N by pattern and method (base enrollment and dropout) ----
base_cell <- long |> filter(enrol == BASE$enrol, dropout == BASE$dropout, pattern != "a_null") |>
  mutate(within = ifelse(duration_mean <= dc$dur_max, "Expected duration within limit",
                         "Expected duration over limit"))
fig3 <- ggplot(base_cell, aes(N, power, colour = method_label)) +
  geom_hline(yintercept = POWER_TARGET, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = power - 2 * mcse_power, ymax = power + 2 * mcse_power), width = 25) +
  geom_point(aes(shape = within), size = 2.4, fill = "white", stroke = 0.9) +
  facet_wrap(~pattern_label, nrow = 1) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_shape_manual(values = c("Expected duration within limit" = 16,
                                "Expected duration over limit" = 21), name = NULL) +
  scale_x_continuous(breaks = c(400, 600, 800, 1000, 1200)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample size N", y = "Power (one-sided α = 0.025)",
       title = "Figure 3. Power vs N, log-rank vs FH(0,1)",
       subtitle = sprintf("%s, %s; D = %d deaths, cap %d months. Bars: ±2 MCSE. Dashed: 80%%.",
                          enrol_lab(BASE$enrol), dropout_lab(BASE$dropout), dc$D, dc$cap)) +
  theme_sim() + theme(legend.position = "top", legend.box = "vertical")
save_fig(fig3, "fig3_power_vs_n.png", h = 4.8)

# ---- Fig 4: type I error with MCSE bands --------------------------------------
null_df <- long |> filter(pattern == "a_null") |> arrange(enrol, dropout) |>
  mutate(panel = paste0(enrol_lab(enrol), "\n", dropout_lab(dropout)),
         panel = factor(panel, levels = unique(panel)))
fig4 <- ggplot(null_df, aes(factor(N), power, colour = method_label)) +
  geom_hline(yintercept = 0.025, colour = "grey30") +
  geom_hline(yintercept = dc$type1_threshold, linetype = "dashed", colour = "grey30") +
  geom_pointrange(aes(ymin = power - 2 * mcse_power, ymax = power + 2 * mcse_power),
                  position = position_dodge(width = 0.6), size = 0.25) +
  facet_wrap(~panel, nrow = 2) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  labs(x = "Sample size N", y = "Type I error",
       title = "Figure 4. Type I error under the null (10,000 replicates per scenario)",
       subtitle = sprintf("Bars: ±2 MCSE. Solid: 2.5%%. Dashed: eligibility threshold %.2f%%.",
                          100 * dc$type1_threshold)) +
  theme_sim() + theme(legend.position = "top", axis.text.x = element_text(size = 7))
save_fig(fig4, "fig4_type1_error.png", h = 5.5)

# ---- Fig 5: sensitivity of power to dropout and enrollment (12-month delay) ---
sens_df <- long |> filter(pattern == BASE$pattern) |>
  mutate(enrol_f = enrol_lab(enrol), dropout_f = factor(dropout_lab(dropout),
                                                       levels = dropout_lab(c(0, 0.05, 0.15))))
fig5 <- ggplot(sens_df, aes(N, power, colour = method_label)) +
  geom_hline(yintercept = POWER_TARGET, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = power - 2 * mcse_power, ymax = power + 2 * mcse_power, fill = method_label),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  facet_grid(enrol_f ~ dropout_f) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_fill_manual(values = METHOD_COLOURS, guide = "none") +
  scale_x_continuous(breaks = c(400, 800, 1200)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample size N", y = "Power",
       title = "Figure 5. Sensitivity of power to enrollment rate and dropout, 12-month delay",
       subtitle = "Bands: ±2 MCSE. Dashed: 80%.") +
  theme_sim() + theme(legend.position = "top")
save_fig(fig5, "fig5_sensitivity_enrol_dropout.png", h = 5.2)

# ---- Fig 6: expected duration vs N --------------------------------------------
dur_df <- s |> filter(pattern == BASE$pattern) |>
  mutate(dropout_f = factor(dropout_lab(dropout), levels = dropout_lab(c(0, 0.05, 0.15))),
         enrol_f = enrol_lab(enrol))
fig6 <- ggplot(dur_df, aes(N, duration_mean, colour = dropout_f)) +
  geom_hline(yintercept = dc$dur_max, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = dc$cap, colour = "grey40") +
  geom_text(data = data.frame(y = c(dc$dur_max, dc$cap), label = c("Duration limit", "Cap"),
                              vj = c(1.6, -0.8)),
            aes(x = 400, y = y, label = label, vjust = vj), inherit.aes = FALSE,
            hjust = 0, size = 3, colour = "grey30") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~enrol_f) +
  scale_colour_manual(values = SERIES_COLOURS[1:3], name = NULL) +
  scale_x_continuous(breaks = c(400, 600, 800, 1000, 1200)) +
  scale_y_continuous(limits = c(40, 90), breaks = seq(40, 90, 12)) +
  labs(x = "Sample size N", y = "Expected study duration (months)",
       title = "Figure 6. Expected study duration vs N, 12-month delay",
       subtitle = sprintf("Analysis at %d deaths or %d months, whichever comes first.", dc$D, dc$cap)) +
  theme_sim() + theme(legend.position = "top")
save_fig(fig6, "fig6_duration_vs_n.png")

# ---- Knit the report ----------------------------------------------------------
knitr::opts_knit$set(root.dir = getwd())
knitr::knit("report/report.Rmd", output = "report/report.md", quiet = TRUE, envir = new.env())
cat(sprintf("Recommendation: %s, N = %d (%s); report written to report/report.md\n",
            METHODS[[rec$method]], rec$N, rec$status))
