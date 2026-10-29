# Stage 1: explore survival::colon (overall survival) and calibrate the control-arm model.
# Outputs: figures/01_*.png, results/stage1_summary.rds, results/control_model.rds

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
source("R/control_model.R")
source("R/plot_theme.R")

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- Data ---------------------------------------------------------------------
os <- colon |>
  filter(etype == 2) |>
  mutate(months = time / DAYS_PER_MONTH)

stopifnot(nrow(os) == 929, anyDuplicated(os$id) == 0)

data_summary <- os |>
  group_by(rx) |>
  summarise(n = n(), deaths = sum(status), .groups = "drop")

# ---- KM by arm ----------------------------------------------------------------
km_tidy <- function(fit) {
  strata <- if (is.null(fit$strata)) rep("All", length(fit$time)) else
    rep(sub(".*=", "", names(fit$strata)), fit$strata)
  first <- data.frame(group = unique(strata), time = 0, surv = 1, lower = 1, upper = 1)
  rbind(first, data.frame(group = strata, time = fit$time, surv = fit$surv,
                          lower = fit$lower, upper = fit$upper))
}

km_arm <- survfit(Surv(months, status) ~ rx, data = os)
km_df <- km_tidy(km_arm) |> mutate(group = factor(group, levels = levels(os$rx)))

p_km <- ggplot(km_df, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.7) +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  scale_x_continuous(breaks = seq(0, 108, 12)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Months since randomization", y = "Overall survival",
       title = "Overall survival by arm, survival::colon") +
  theme_sim() + theme(legend.position = "top")
ggsave("figures/01_km_by_arm.png", p_km, width = 7, height = 4.5, dpi = 150)

# ---- Follow-up distribution (reverse KM) --------------------------------------
rkm <- survfit(Surv(months, 1 - status) ~ 1, data = os)
fu_quantiles <- quantile(rkm, probs = c(0.25, 0.5, 0.75))$quantile
rkm_df <- km_tidy(rkm)

p_fu <- ggplot(rkm_df, aes(time, 1 - surv)) +
  geom_step(linewidth = 0.7, colour = ARM_COLOURS[[1]]) +
  geom_vline(xintercept = fu_quantiles[["50"]], linetype = "dashed", colour = "grey40") +
  annotate("text", x = fu_quantiles[["50"]], y = 0.05, hjust = 1.05, size = 3.2,
           label = sprintf("Median follow-up %.1f mo", fu_quantiles[["50"]])) +
  scale_x_continuous(breaks = seq(0, 108, 12)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Months since randomization", y = "Cumulative proportion censored",
       title = "Follow-up distribution (reverse Kaplan-Meier)") +
  theme_sim()
ggsave("figures/01_followup_reverse_km.png", p_fu, width = 7, height = 4, dpi = 150)

# ---- Events per year of follow-up ---------------------------------------------
yr_cuts <- seq(12, 108, 12)
os_yr <- survSplit(Surv(months, status) ~ rx, data = os, cut = yr_cuts,
                   start = "t0", episode = "year")
events_per_year <- os_yr |>
  mutate(pt = months - t0) |>
  group_by(rx, year) |>
  summarise(deaths = sum(status), person_years = sum(pt) / 12, .groups = "drop") |>
  filter(person_years > 0) |>
  mutate(rate_per_100py = 100 * deaths / person_years)

p_epy <- ggplot(filter(events_per_year, year <= 8), aes(factor(year), deaths, fill = rx)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = ARM_COLOURS, name = NULL) +
  labs(x = "Year of follow-up", y = "Deaths",
       title = "Deaths per year of follow-up, by arm") +
  theme_sim() + theme(legend.position = "top")
ggsave("figures/01_events_per_year.png", p_epy, width = 7, height = 4, dpi = 150)

# ---- Cox HR: Lev+5FU vs Obs, with PH check ------------------------------------
os2 <- os |> filter(rx %in% c("Obs", "Lev+5FU")) |> mutate(rx = droplevels(rx))
cox <- coxph(Surv(months, status) ~ rx, data = os2)
cox_hr <- c(hr = exp(coef(cox))[[1]], exp(confint(cox))[1, ])
names(cox_hr) <- c("hr", "lower", "upper")
cox_p <- summary(cox)$coefficients[1, "Pr(>|z|)"]
zph <- cox.zph(cox)
zph_p <- zph$table["rx", "p"]

# Period-specific HRs, to see whether the observed effect looks delayed or proportional.
hr_periods <- c(0, 12, 36)
os2_split <- survSplit(Surv(months, status) ~ rx, data = os2, cut = hr_periods[-1],
                       start = "t0", episode = "period")
os2_split <- os2_split |>
  mutate(trt = as.numeric(rx == "Lev+5FU"),
         trt_p1 = trt * (period == 1), trt_p2 = trt * (period == 2), trt_p3 = trt * (period == 3))
cox_period <- coxph(Surv(t0, months, status) ~ trt_p1 + trt_p2 + trt_p3, data = os2_split)
period_hr <- data.frame(
  period = c("0-12", "12-36", "36+"),
  deaths = as.vector(tapply(os2_split$status, os2_split$period, sum)),
  hr = exp(coef(cox_period)),
  lower = exp(confint(cox_period))[, 1],
  upper = exp(confint(cox_period))[, 2],
  row.names = NULL
)

zph_df <- data.frame(time = zph$time, beta = zph$y[, 1])
p_zph <- ggplot(zph_df, aes(time, beta)) +
  geom_point(alpha = 0.35, size = 1.2, colour = "grey45") +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = ARM_COLOURS[[3]], fill = ARM_COLOURS[[3]], alpha = 0.15, linewidth = 0.8) +
  geom_hline(yintercept = coef(cox), linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey60") +
  scale_x_continuous(breaks = seq(0, 108, 12)) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "Months since randomization", y = "Scaled Schoenfeld residual, beta(t)",
       title = "Proportional-hazards check, Lev+5FU vs Obs",
       subtitle = sprintf("Cox log HR (dashed) = %.2f; cox.zph p = %.2f", coef(cox), zph_p)) +
  theme_sim()
ggsave("figures/01_cox_zph.png", p_zph, width = 7, height = 4, dpi = 150)

# ---- Control-arm calibration (Obs) --------------------------------------------
obs <- os |> filter(rx == "Obs")
pwe_cuts <- c(0, 6, 12, 24, 48)
pwe <- pwe_fit(obs$months, obs$status, pwe_cuts)
wei <- weibull_fit(obs$months, obs$status)
aic <- c(pwe = -2 * pwe$loglik + 2 * pwe$n_par, weibull = -2 * wei$loglik + 2 * wei$n_par)

km_obs <- survfit(Surv(months, status) ~ 1, data = obs)
grid_t <- seq(0, max(obs$months), length.out = 400)
fit_curves <- bind_rows(
  data.frame(time = grid_t, model = "Piecewise exponential",
             surv = pwe_surv(grid_t, pwe$cuts, pwe$rate),
             hazard = pwe_hazard(grid_t, pwe$cuts, pwe$rate)),
  data.frame(time = grid_t, model = "Weibull",
             surv = weibull_surv(grid_t, wei$shape, wei$scale),
             hazard = weibull_hazard(pmax(grid_t, 0.01), wei$shape, wei$scale))
)
model_colours <- c("Piecewise exponential" = ARM_COLOURS[[1]], "Weibull" = ARM_COLOURS[[2]])

# Largest vertical gap between each fitted curve and the KM, over the KM's event times.
km_gap <- c(
  pwe = max(abs(km_obs$surv - pwe_surv(km_obs$time, pwe$cuts, pwe$rate))),
  weibull = max(abs(km_obs$surv - weibull_surv(km_obs$time, wei$shape, wei$scale)))
)

p_cal_s <- ggplot() +
  geom_ribbon(data = km_tidy(km_obs), aes(time, ymin = lower, ymax = upper),
              fill = "grey85", stat = "identity") +
  geom_step(data = km_tidy(km_obs), aes(time, surv), linewidth = 0.6) +
  geom_line(data = fit_curves, aes(time, surv, colour = model), linewidth = 0.8) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_x_continuous(breaks = seq(0, 108, 12)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Months since randomization", y = "Overall survival",
       title = "Obs arm: Kaplan-Meier (black, 95% CI) vs fitted models") +
  theme_sim() + theme(legend.position = "top")

# Empirical hazard in 6-month bins with exact Poisson 95% CI.
haz_cuts <- seq(0, 108, 6)
emp_haz <- pwe_fit(obs$months, obs$status, haz_cuts)
emp_haz_df <- data.frame(
  mid = haz_cuts + 3, d = emp_haz$events, pt = emp_haz$person_time
) |>
  filter(pt > 50) |>
  mutate(rate = d / pt,
         lower = qchisq(0.025, 2 * d) / 2 / pt,
         upper = qchisq(0.975, 2 * (d + 1)) / 2 / pt)

p_cal_h <- ggplot() +
  geom_pointrange(data = emp_haz_df, aes(mid, 12 * rate, ymin = 12 * lower, ymax = 12 * upper),
                  colour = "grey35", size = 0.25) +
  geom_line(data = fit_curves, aes(time, 12 * hazard, colour = model), linewidth = 0.8) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_x_continuous(breaks = seq(0, 108, 12)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  labs(caption = "Y-axis truncated at 0.4; the last bins have little person-time and wide CIs.",
       x = "Months since randomization", y = "Hazard (deaths per person-year)",
       title = "Obs arm hazard: 6-month empirical bins (95% CI) vs fitted models") +
  theme_sim() + theme(legend.position = "none")

ggsave("figures/01_obs_hazard.png", p_cal_h + theme(legend.position = "top"), width = 7, height = 4, dpi = 150)
ggsave("figures/01_control_calibration.png", p_cal_s / p_cal_h, width = 7, height = 8, dpi = 150)

# ---- Save ---------------------------------------------------------------------
control_model <- list(type = "pwe", cuts = pwe$cuts, rate = pwe$rate,
                      source = "survival::colon, Obs arm, OS (etype == 2), time in months")
saveRDS(control_model, "results/control_model.rds")

stage1 <- list(
  data_summary = data_summary,
  followup_quantiles = fu_quantiles,
  events_per_year = events_per_year,
  cox_hr = cox_hr, cox_p = cox_p, cox_zph_p = zph_p, period_hr = period_hr,
  pwe = pwe, weibull = wei, aic = aic, km_gap = km_gap,
  km_obs_at = summary(km_obs, times = c(12, 24, 36, 60))[c("time", "surv")],
  pwe_surv_at = pwe_surv(c(12, 24, 36, 60), pwe$cuts, pwe$rate),
  weibull_surv_at = weibull_surv(c(12, 24, 36, 60), wei$shape, wei$scale)
)
saveRDS(stage1, "results/stage1_summary.rds")

# ---- Console summary ----------------------------------------------------------
print(data_summary)
cat(sprintf("\nFollow-up (reverse KM), months: Q1 %.1f, median %.1f, Q3 %.1f\n",
            fu_quantiles[1], fu_quantiles[2], fu_quantiles[3]))
cat(sprintf("Cox HR Lev+5FU vs Obs: %.3f (95%% CI %.3f-%.3f), p = %.2g; cox.zph p = %.3f\n",
            cox_hr["hr"], cox_hr["lower"], cox_hr["upper"], cox_p, zph_p))
cat("\nPeriod-specific HRs:\n"); print(period_hr, digits = 3)
cat("\nPWE monthly hazard by interval:\n")
print(data.frame(from = pwe$cuts, to = c(pwe$cuts[-1], Inf), events = pwe$events,
                 rate_per_month = signif(pwe$rate, 3), rate_per_year = round(12 * pwe$rate, 3)))
cat(sprintf("\nWeibull: shape %.3f, scale %.1f months\n", wei$shape, wei$scale))
cat("AIC:\n"); print(round(aic, 1))
cat("Max |KM - fit|:\n"); print(round(km_gap, 3))
cat("\nS(t) at 12/24/36/60 months:\n")
print(round(rbind(KM = stage1$km_obs_at$surv, PWE = stage1$pwe_surv_at,
                  Weibull = stage1$weibull_surv_at), 3))
cat("\nEvents per year (all arms):\n")
print(events_per_year |> group_by(year) |>
        summarise(deaths = sum(deaths), person_years = round(sum(person_years))))
