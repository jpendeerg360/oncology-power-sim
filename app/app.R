# Exploratory design calculator: simulated and analytic power for one user-chosen design.
# Run locally from the project root with shiny::runApp("app"). Deploy with app/deploy.R.
#
# All simulation and analysis code comes from the project's R/ directory (tested in
# tests/testthat/), so the app adds only the interface.

library(shiny)
library(bslib)

# R/ and the control model sit next to app.R in a deployed bundle, or one level up locally.
root <- if (dir.exists("R")) "." else ".."
for (f in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
control <- readRDS(file.path(root, "results", "control_model.rds"))
future::plan(future::sequential)

# Pre-specified base case (docs/design_spec.md section 9) and the recommended design.
BASE <- list(N = 900, enrol = 30, dropout = 5, hr = 0.69, delay = 12, D = 306, cap = 84,
             reps = 1000, seed = MASTER_SEED)
MAX_REPS <- 2000

pct <- function(p, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * p)

ui <- page_sidebar(
  title = "Delayed-effect OS trial: design calculator",
  theme = bs_theme(version = 5, primary = "#2a78d6"),
  sidebar = sidebar(
    width = 300,
    numericInput("N", "Sample size N (even)", BASE$N, min = 100, max = 3000, step = 50),
    numericInput("enrol", "Enrollment (patients/month)", BASE$enrol, min = 1, max = 200, step = 5),
    numericInput("dropout", "Dropout (%/year)", BASE$dropout, min = 0, max = 50, step = 1),
    numericInput("hr", "Hazard ratio after delay", BASE$hr, min = 0.3, max = 1, step = 0.01),
    numericInput("delay", "Delay before effect (months)", BASE$delay, min = 0, max = 36, step = 1),
    numericInput("D", "Target deaths D", BASE$D, min = 20, max = 2000, step = 10),
    numericInput("cap", "Calendar cap (months)", BASE$cap, min = 12, max = 240, step = 6),
    selectInput("reps", "Simulated trials", c(500, 1000, 2000), selected = BASE$reps),
    numericInput("seed", "Seed", BASE$seed, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary"),
    actionButton("reset", "Reset to pre-specified base case", class = "btn-outline-secondary btn-sm")
  ),
  div(
    class = "alert alert-warning mb-0",
    strong("Exploratory tool. "),
    "The pre-specified study, its decision rule and recommendation (N = 900, FH(0,1), ",
    "306 deaths or 84 months) are in the ",
    tags$a("project report", href = "https://github.com/jpendeerg360/oncology-power-sim/blob/main/report/report.md",
           target = "_blank"),
    ". This calculator checks designs beyond that grid; its results are not part of the pre-specified analysis."
  ),
  layout_columns(
    col_widths = c(7, 5),
    card(
      card_header("Power: simulation vs analytic approximation"),
      tableOutput("power_tbl"),
      uiOutput("sim_note")
    ),
    card(
      card_header("Trial conduct"),
      tableOutput("conduct_tbl")
    )
  ),
  card(
    card_header("Assumed survival curves"),
    plotOutput("surv_plot", height = "280px")
  ),
  card(
    card_header("Methods"),
    p("Each simulated trial uses the study's generator (", code("simulate_trial()"), "), event-driven cut (",
      code("cut_data()"), ": D deaths or the cap, whichever comes first) and analyses (", code("analyze_trial()"),
      "): one-sided log-rank and Fleming–Harrington FH(0,1) at α = 0.025. The control arm is the ",
      "five-interval piecewise-exponential fit to the observation arm of ", code("survival::colon"), "."),
    p("Simulated power is shown with its Monte Carlo standard error, √(p(1 − p)/R). Results are reproducible ",
      "for a fixed seed. A grid scenario's own seed and replicate count reproduce it exactly: the base case is ",
      "seed 20261149 with 2,000 trials (results/scenario_grid.csv lists the rest)."),
    p("Analytic columns use ", code("gsDesign2"), " at the expected analysis time: the average-hazard-ratio ",
      "method for log-rank (", code("fixed_design_ahr()"), ") and its weighted analogue for FH(0,1) (",
      code("fixed_design_fh()"), "). Schoenfeld uses D deaths and the post-delay HR, so it ignores the delay ",
      "and is the planning value a PH design would assume. Accrual is uniform in the formulas and Poisson ",
      "in the simulation.")
  )
)

server <- function(input, output, session) {
  observeEvent(input$reset, {
    for (id in c("N", "enrol", "dropout", "hr", "delay", "D", "cap", "seed"))
      updateNumericInput(session, id, value = BASE[[id]])
    updateSelectInput(session, "reps", selected = BASE$reps)
  })

  design <- reactive({
    validate(
      need(is.numeric(input$N) && input$N >= 100 && input$N %% 2 == 0, "N must be an even number ≥ 100."),
      need(input$enrol > 0, "Enrollment rate must be positive."),
      need(input$dropout >= 0 && input$dropout < 100, "Dropout must be between 0 and 100%."),
      need(input$hr > 0 && input$hr <= 1, "HR must be in (0, 1]."),
      need(input$delay >= 0 && input$delay < input$cap, "Delay must be between 0 and the cap."),
      need(input$D >= 10 && input$D <= input$N, "Target deaths must be between 10 and N."),
      need(input$cap > 0, "Cap must be positive."),
      need(is.numeric(input$seed) && !is.na(input$seed), "Seed must be a number.")
    )
    list(N = input$N, enrol = input$enrol, dropout = input$dropout / 100, hr = input$hr,
         delay = input$delay, D = input$D, cap = input$cap,
         reps = min(as.integer(input$reps), MAX_REPS), seed = as.integer(input$seed))
  }) |> debounce(500)

  analytic <- reactive({
    d <- design()
    analytic_design(d$N, d$enrol, d$dropout, d$hr, d$delay, d$D, d$cap, control)
  })

  # Simulation runs on demand only; results are kept with the design that produced them.
  sim <- reactiveVal(NULL)
  run_sim <- function() {
    d <- design()
    withProgress(message = sprintf("Simulating %d trials", d$reps), value = 0.5, {
      s <- simulate_design(d$N, d$enrol, d$dropout, d$hr, d$delay, d$D, d$cap, control,
                           reps = d$reps, seed = d$seed)
    })
    sim(list(design = d, summary = s))
  }
  observeEvent(input$run, run_sim())
  observeEvent(design(), if (is.null(sim())) run_sim(), once = TRUE)

  stale <- reactive(!is.null(sim()) && !identical(sim()$design, design()))

  output$power_tbl <- renderTable({
    s <- sim()$summary; a <- analytic()
    req(s)
    is_null <- design()$hr == 1
    data.frame(
      Test = c("Log-rank", "FH(0,1)"),
      `Simulated` = c(sprintf("%s ± %s", pct(s$power_lr), pct(s$mcse_power_lr)),
                      sprintf("%s ± %s", pct(s$power_fh), pct(s$mcse_power_fh))),
      `AHR formula` = c(pct(a$power_lr_ahr), pct(a$power_fh_ahr)),
      `Schoenfeld (PH)` = c(pct(a$power_schoenfeld), "—"),
      check.names = FALSE
    ) |> setNames(c("Test", if (is_null) "Simulated type I error ± MCSE" else "Simulated ± MCSE",
                   "AHR formula", "Schoenfeld (PH)"))
  }, striped = TRUE, width = "100%", align = "lccc")

  output$sim_note <- renderUI({
    s <- sim()$summary
    req(s)
    diff <- sprintf("FH(0,1) − log-rank: %+.1f points (MCSE %.1f, paired on the same trials).",
                    100 * s$diff_fh_lr, 100 * s$mcse_diff)
    tagList(
      p(class = "text-muted small mb-1", sprintf("R = %d simulated trials, seed %d. %s", s$R, s$seed, diff)),
      if (stale()) p(class = "text-danger small mb-0",
                     "Inputs have changed since this simulation. Analytic columns are current; press Run simulation to update.")
    )
  })

  output$conduct_tbl <- renderTable({
    s <- sim()$summary; a <- analytic()
    req(s)
    data.frame(
      Quantity = c("Analysis time (months)", "Analyses triggered by cap", "Deaths at analysis",
                   "Patients enrolled at analysis", "Average HR over follow-up"),
      Simulated = c(sprintf("%.1f ± %.1f", s$duration_mean, s$mcse_duration),
                    sprintf("%s ± %s", pct(s$pct_cap), pct(s$mcse_pct_cap)),
                    sprintf("%.0f", s$deaths_mean), sprintf("%.0f", s$enrolled_mean),
                    sprintf("%.3f (Cox)", exp(s$loghr_mean))),
      Analytic = c(sprintf("%.1f", a$time), if (a$time >= design()$cap) "cap" else "—",
                   sprintf("%.0f", a$events), sprintf("%.0f", min(design()$N, design()$enrol * a$time)),
                   sprintf("%.3f (AHR)", a$ahr))
    )
  }, striped = TRUE, width = "100%", align = "lcc")

  output$surv_plot <- renderPlot({
    d <- design()
    t <- seq(0, d$cap, by = 0.5)
    trt <- treatment_hazard(control$cuts, control$rate, d$hr, d$delay)
    df <- rbind(data.frame(t = t, S = pwe_surv(t, control$cuts, control$rate), arm = "Control"),
                data.frame(t = t, S = pwe_surv(t, trt$cuts, trt$rate), arm = "Treatment"))
    ggplot2::ggplot(df, ggplot2::aes(t, S, colour = arm)) +
      ggplot2::geom_vline(xintercept = d$delay, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_colour_manual(values = SERIES_COLOURS[1:2], name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
      ggplot2::labs(x = "Months since randomization", y = "Overall survival",
                    subtitle = sprintf("HR = 1 for the first %g months, then %.2f (dashed line)", d$delay, d$hr)) +
      theme_sim(13) + ggplot2::theme(legend.position = "top")
  })
}

shinyApp(ui, server)
