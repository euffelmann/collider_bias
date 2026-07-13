# Collider Bias Explorer
#
# An interactive demo of collider bias (Berkson's paradox): a sample gets
# selected based on a third variable that is itself caused by / correlated
# with two other variables X and Y. Conditioning on that third variable can
# induce a spurious correlation between X and Y, even when none exists in
# the full population.
#
# Two tabs:
#   - Continuous variables: X and Y are correlated with a general selection
#     variable S (e.g. X = IQ, Y = chess skill, S = becoming a chess
#     grandmaster).
#   - Genetic (SNPs): X and Y are two LD-independent SNPs that each explain
#     some variance in the liability of a disease (e.g. Alzheimer's) under
#     the liability-threshold model. The disease (case status) is the
#     collider: ascertaining a case-enriched sample induces a spurious
#     SNP-SNP correlation even though the SNPs are independent in the
#     population.

library(shiny)
library(ggplot2)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Simple lm fit summary, guarding against degenerate subsets.
fit_line <- function(data) {
  if (nrow(data) < 3 || sd(data$x) == 0 || sd(data$y) == 0) {
    return(list(intercept = NA_real_, slope = NA_real_, r = NA_real_, n = nrow(data)))
  }
  m <- lm(y ~ x, data = data)
  list(intercept = unname(coef(m)[1]),
       slope = unname(coef(m)[2]),
       r = cor(data$x, data$y),
       n = nrow(data))
}

fmt_r <- function(fit) if (is.na(fit$r)) "NA" else sprintf("%.2f", fit$r)

COLS <- c("Not selected" = "grey70", "Selected" = "#D7191C")

# Build the scatter + two regression lines (population vs. selected).
# `jitter` is used for the SNP tab where x/y only take values 0/1/2.
# `df$group` must be a 2-level factor; the 2nd level is drawn in red.
make_plot <- function(stats, jitter = FALSE, xlab = "X", ylab = "Y",
                       cols = COLS, legend_title = "Selection status") {
  df <- stats$df
  pop_fit <- stats$pop_fit
  sel_fit <- stats$sel_fit
  sel_colour <- unname(cols[levels(df$group)[2]])

  p <- ggplot(df, aes(x = x, y = y))

  if (jitter) {
    p <- p + geom_jitter(aes(colour = group), width = 0.12, height = 0.12,
                          alpha = 0.55, size = 1.8)
  } else {
    p <- p + geom_point(aes(colour = group), alpha = 0.45, size = 1.6)
  }

  p <- p +
    (if (!is.na(pop_fit$slope))
      geom_abline(intercept = pop_fit$intercept, slope = pop_fit$slope,
                  colour = "black", linetype = "dashed", linewidth = 1)
     else NULL) +
    (if (!is.na(sel_fit$slope))
      geom_abline(intercept = sel_fit$intercept, slope = sel_fit$slope,
                  colour = sel_colour, linewidth = 1.1)
     else NULL) +
    scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")

  p
}

# HTML summary of the two correlations, colour-matched to the plot legend,
# shown below the plot (avoids clipping that a ggplot subtitle suffers at
# typical plot widths).
make_stats_html <- function(stats, labels = c("Full population", "Selected sample"),
                             cols = COLS) {
  sel_colour <- unname(cols[levels(stats$df$group)[2]])
  line <- function(label, fit, colour) {
    sprintf(
      '<span style="color:%s; font-weight:600;">%s</span>: r = %s (n = %d)',
      colour, label, fmt_r(fit), fit$n
    )
  }
  HTML(paste(
    line(labels[1], stats$pop_fit, "black"),
    line(labels[2], stats$sel_fit, sel_colour),
    sep = "&nbsp;&nbsp;|&nbsp;&nbsp;"
  ))
}

# ---------------------------------------------------------------------------
# Continuous-variables tab: simulation
# ---------------------------------------------------------------------------

# Correlated pair of standard normals. Using a fixed seed keyed only on `n`
# (not on rho) means the underlying draws stay put while rho/selection
# sliders are moved, so the scatter morphs smoothly instead of re-randomizing.
simulate_latent <- function(n, rho, seed) {
  set.seed(seed)
  z1 <- rnorm(n)
  noise <- rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * noise
  data.frame(x = z1, y = z2)
}

# Add a selection variable S = a*x + b*y + c*noise (all standardized),
# with weights (a, b) chosen so that, given the population correlation
# rho_xy between x and y, S achieves the target correlations rho_xs and
# rho_ys with x and y respectively. If the requested combination of
# rho_xy/rho_xs/rho_ys is not jointly achievable (correlation matrix would
# not be positive semi-definite), (a, b) are rescaled to the nearest
# feasible point (S becomes a deterministic function of x, y with no extra
# noise) rather than erroring.
add_selection_var <- function(df, rho_xy, rho_xs, rho_ys) {
  denom <- 1 - rho_xy^2
  if (denom < 1e-6) denom <- 1e-6
  a <- (rho_xs - rho_xy * rho_ys) / denom
  b <- (rho_ys - rho_xy * rho_xs) / denom
  var_ab <- a^2 + b^2 + 2 * a * b * rho_xy
  if (var_ab > 1) {
    scale <- sqrt(1 / var_ab)
    a <- a * scale
    b <- b * scale
    var_ab <- 1
  }
  c_coef <- sqrt(max(0, 1 - var_ab))
  df$s <- a * df$x + b * df$y + c_coef * rnorm(nrow(df))
  df
}

simulate_continuous <- function(n, rho, rho_xs, rho_ys, seed) {
  latent <- simulate_latent(n, rho, seed)
  add_selection_var(latent, rho, rho_xs, rho_ys)
}

# Classify points by selection status based on the selection variable S:
# top select_pct% of S (marginally) are "Selected" (e.g. became a chess
# grandmaster). Lower % = more stringent selection.
classify_selection <- function(df, select_pct) {
  thr_s <- quantile(df$s, probs = 1 - select_pct / 100, type = 1, names = FALSE)
  df$group <- factor(ifelse(df$s >= thr_s, "Selected", "Not selected"),
                      levels = c("Not selected", "Selected"))
  attr(df, "thr_s") <- thr_s
  df
}

# Classify + fit both regression lines once; shared by the plot and the
# text summary so the (identical) computation isn't duplicated.
compute_stats <- function(df, select_pct) {
  df <- classify_selection(df, select_pct)
  sel_data <- df[df$group == "Selected", , drop = FALSE]
  list(
    df = df,
    pop_fit = fit_line(df),
    sel_fit = fit_line(sel_data)
  )
}

explanation_cont <- paste(
  "Collider bias (Berkson's paradox): X and Y are simulated with a fixed,",
  "known correlation in the full population (black dashed line). A third",
  "variable S - the selection variable - is correlated with X and with Y",
  "by the amounts set below (e.g. X = IQ, Y = chess skill, S = becoming a",
  "chess grandmaster). If the sample is restricted to those with high S",
  "(red points/line), a spurious correlation can appear between X and Y",
  "even when little or none exists in the population - and it can even",
  "flip sign - purely because S depends on both."
)

# ---------------------------------------------------------------------------
# Genetic (SNPs) tab: simulation
#
# X (SNP1) and Y (SNP2) are independent (no LD) in the population. Each
# influences a binary disease via an underlying, standard-normal liability
# (Falconer's liability-threshold model): liability = b1*SNP1 + b2*SNP2 +
# environmental noise, with b_i chosen so SNP_i alone explains r2_i of the
# liability variance. Individuals with liability above a threshold t (set
# by the population prevalence K) are cases.
#
# A large background "pool" of genotypes + noise is drawn once per seed and
# reused as sliders move, so the scatter morphs smoothly rather than
# re-randomizing (same rationale as simulate_latent() above). Bigger is
# better for guaranteeing enough cases at low prevalence + high sample size
# + P near 1, but it's kept modest (rather than e.g. 2e7) to fit shinyapps.io
# memory limits - build_snp_stats() degrades gracefully (returns fewer rows
# than requested) rather than erroring if the pool runs short.
# ---------------------------------------------------------------------------

SNP_POOL_SIZE <- 2e6

simulate_disease_pool <- function(pool_size, maf1, maf2, seed) {
  set.seed(seed)
  data.frame(
    x = rbinom(pool_size, 2, maf1),
    y = rbinom(pool_size, 2, maf2),
    e = rnorm(pool_size)
  )
}

# Attach a liability column, given each SNP's target R^2 on the liability
# scale and effect direction (dir_i = +1 risk-increasing, -1 risk-decreasing
# per copy of the coded/"risk" allele). If r2_1 + r2_2 would leave no
# residual variance, they're rescaled down proportionally so the liability
# keeps unit variance.
add_liability <- function(pool, maf1, maf2, r2_1, r2_2, dir1 = 1, dir2 = 1) {
  total <- r2_1 + r2_2
  if (total > 0.98) {
    scale <- 0.98 / total
    r2_1 <- r2_1 * scale
    r2_2 <- r2_2 * scale
  }
  beta1 <- dir1 * sqrt(r2_1)
  beta2 <- dir2 * sqrt(r2_2)
  resid_sd <- sqrt(max(1 - r2_1 - r2_2, 0.001))
  mu1 <- 2 * maf1; sd1 <- sqrt(2 * maf1 * (1 - maf1))
  mu2 <- 2 * maf2; sd2 <- sqrt(2 * maf2 * (1 - maf2))
  g1_std <- (pool$x - mu1) / sd1
  g2_std <- (pool$y - mu2) / sd2
  pool$liability <- beta1 * g1_std + beta2 * g2_std + resid_sd * pool$e
  pool
}

# Build the two groups shown in the plot:
#   - "Population": n individuals drawn at random (disease at its natural
#     population prevalence K, no ascertainment).
#   - "Selected": a case-control sample of n individuals ascertained to a
#     target case fraction P (the "sample prevalence"). P = K reproduces
#     the population (no selection); P towards 1 approaches a cases-only
#     sample (maximal selection on case status).
build_snp_stats <- function(pool, n, K, P) {
  t <- -qnorm(K)
  pool$case <- pool$liability > t

  pop_df <- pool[seq_len(min(n, nrow(pool))), c("x", "y")]

  n_cases_target <- round(n * P)
  n_controls_target <- n - n_cases_target
  cases_pool <- pool[pool$case, c("x", "y"), drop = FALSE]
  controls_pool <- pool[!pool$case, c("x", "y"), drop = FALSE]
  n_cases <- min(n_cases_target, nrow(cases_pool))
  n_controls <- min(n_controls_target, nrow(controls_pool))
  sel_df <- rbind(
    cases_pool[seq_len(n_cases), , drop = FALSE],
    controls_pool[seq_len(n_controls), , drop = FALSE]
  )

  pop_df$group <- rep("Not selected", nrow(pop_df))
  sel_df$group <- rep("Selected", nrow(sel_df))
  df <- rbind(pop_df, sel_df)
  df$group <- factor(df$group, levels = c("Not selected", "Selected"))

  list(
    df = df,
    pop_fit = fit_line(pop_df),
    sel_fit = fit_line(sel_df),
    n_cases = n_cases,
    n_controls = n_controls
  )
}

# --- Liability-scale <-> observed-scale R^2 conversion (Lee et al. 2012,
# Genet Epidemiology), used purely for the informational display below the
# plot: it shows how case-control ascertainment inflates the *apparent*
# per-SNP effect size relative to its true, ascertainment-invariant
# liability-scale R^2.
prs_r2liab_to_r2obs <- function(K, P, prs_r2liab) {
  t <- -qnorm(K, mean = 0, sd = 1)
  z <- dnorm(t)
  i1 <- z / K
  i0 <- -z / (1 - K)

  theta <- i1 * (P - K) / (1 - K) * (i1 * (P - K) / (1 - K) - t)
  cv <- K * (1 - K) / z^2 * K * (1 - K) / (P * (1 - P))

  prs_r2liab / (cv - prs_r2liab * theta * cv)
}

make_r2_html <- function(K, P, r2_1, r2_2) {
  if (P <= 0 || P >= 1) {
    fmt <- function(label, r2) sprintf("%s: R²(liability) = %.2f", label, r2)
    return(paste0(
      fmt("SNP1", r2_1), "&nbsp;&nbsp;|&nbsp;&nbsp;", fmt("SNP2", r2_2),
      "&nbsp;&nbsp;(R²-observed is undefined for an all-case or all-control sample)"
    ))
  }
  obs1 <- prs_r2liab_to_r2obs(K, P, r2_1)
  obs2 <- prs_r2liab_to_r2obs(K, P, r2_2)
  sprintf(
    paste0(
      "SNP1: R²(liability) = %.2f → R²(observed in this sample) ≈ %.2f",
      "&nbsp;&nbsp;|&nbsp;&nbsp;",
      "SNP2: R²(liability) = %.2f → R²(observed in this sample) ≈ %.2f"
    ),
    r2_1, obs1, r2_2, obs2
  )
}

explanation_snp <- paste(
  "Here X (SNP1) and Y (SNP2) are dosages (0/1/2 copies of the coded",
  "allele) for two genetic variants that are independent in the population",
  "(no LD) but that each explain some variance in the liability of a",
  "disease (e.g. Alzheimer's), following the liability-threshold model.",
  "Each SNP's coded allele can be set to be risk-increasing or",
  "risk-decreasing. The disease is the collider: population prevalence (K)",
  "sets the true disease rate, while sample prevalence (P) sets how",
  "case-enriched your study sample is. P = K means no ascertainment",
  "(matches the population, black dashed line); P = 1 is a cases-only",
  "sample, P = 0 is a controls-only sample. As P moves away from K, a",
  "spurious correlation between SNP1 and SNP2 can appear in the selected",
  "sample (red points/line) even though the SNPs are independent in the",
  "population - and its sign depends on whether both SNPs' coded alleles",
  "move risk in the same direction or opposite directions."
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

continuous_tab <- tabPanel(
  "Continuous variables",
  sidebarLayout(
    sidebarPanel(
      helpText(explanation_cont),
      sliderInput("rho_cont", "True correlation between IQ (X) and chess skill (Y)",
                  min = -0.9, max = 0.9, value = 0, step = 0.05),
      sliderInput("rho_xs_cont", "Correlation between IQ and becoming a chess grandmaster",
                  min = -0.9, max = 0.9, value = 0.6, step = 0.05),
      sliderInput("rho_ys_cont", "Correlation between chess skill and becoming a chess grandmaster",
                  min = -0.9, max = 0.9, value = 0.6, step = 0.05),
      sliderInput("select_pct_cont",
                  "Selection strength (top % who become grandmasters)",
                  min = 1, max = 50, value = 10, step = 1),
      helpText("Lower % = more stringent (stronger) selection."),
      sliderInput("n_cont", "Sample size", min = 500, max = 5000,
                  value = 2000, step = 500),
      actionButton("resample_cont", "Draw new random sample")
    ),
    mainPanel(
      plotOutput("plot_cont", height = "550px"),
      div(style = "margin-top: 10px; font-size: 14px;", htmlOutput("stats_cont"))
    )
  )
)

snp_tab <- tabPanel(
  "Genetic (SNPs)",
  sidebarLayout(
    sidebarPanel(
      helpText(explanation_snp),
      sliderInput("maf1", "SNP1 coded allele frequency",
                  min = 0.05, max = 0.5, value = 0.3, step = 0.05),
      radioButtons("dir1", "SNP1 coded allele effect direction",
                   choices = c("Risk-increasing" = 1, "Risk-decreasing" = -1),
                   selected = 1, inline = TRUE),
      sliderInput("maf2", "SNP2 coded allele frequency",
                  min = 0.05, max = 0.5, value = 0.3, step = 0.05),
      radioButtons("dir2", "SNP2 coded allele effect direction",
                   choices = c("Risk-increasing" = 1, "Risk-decreasing" = -1),
                   selected = 1, inline = TRUE),
      sliderInput("r2_1", "SNP1: variance explained in disease liability (R²)",
                  min = 0, max = 0.5, value = 0.25, step = 0.01),
      sliderInput("r2_2", "SNP2: variance explained in disease liability (R²)",
                  min = 0, max = 0.5, value = 0.25, step = 0.01),
      sliderInput("K_snp", "Population disease prevalence (K)",
                  min = 0.01, max = 0.5, value = 0.1, step = 0.01),
      sliderInput("P_snp",
                  "Sample disease prevalence (P) — case fraction in the selected sample",
                  min = 0, max = 1, value = 1, step = 0.05),
      helpText(paste(
        "P close to K = little/no selection. P = 1 = a cases-only sample.",
        "P = 0 = a controls-only sample."
      )),
      sliderInput("n_snp", "Sample size (per group)", min = 0, max = 50000,
                  value = 5000, step = 1000),
      actionButton("resample_snp", "Draw new random sample")
    ),
    mainPanel(
      plotOutput("plot_snp", height = "550px"),
      div(style = "margin-top: 10px; font-size: 14px;", htmlOutput("stats_snp")),
      div(style = "margin-top: 6px; font-size: 13px; color: #555;", htmlOutput("r2_snp"))
    )
  )
)

ui <- navbarPage(
  title = "Collider Bias Explorer",
  continuous_tab,
  snp_tab
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  seed_cont <- reactiveVal(1)
  observeEvent(input$resample_cont, seed_cont(seed_cont() + 1))

  seed_snp <- reactiveVal(1)
  observeEvent(input$resample_snp, seed_snp(seed_snp() + 1))

  stats_cont <- reactive({
    df <- simulate_continuous(input$n_cont, input$rho_cont,
                               input$rho_xs_cont, input$rho_ys_cont, seed_cont())
    compute_stats(df, input$select_pct_cont)
  })

  # Redrawn only when genotype-generating inputs (MAF) or the seed change;
  # liability/threshold/ascertainment are derived from this pool on the fly
  # as the R^2/prevalence sliders move, so the scatter morphs smoothly.
  snp_pool <- reactive({
    simulate_disease_pool(SNP_POOL_SIZE, input$maf1, input$maf2, seed_snp())
  })

  snp_pool_liab <- reactive({
    add_liability(snp_pool(), input$maf1, input$maf2, input$r2_1, input$r2_2,
                  as.numeric(input$dir1), as.numeric(input$dir2))
  })

  stats_snp <- reactive({
    build_snp_stats(snp_pool_liab(), input$n_snp, input$K_snp, input$P_snp)
  })

  output$plot_cont <- renderPlot({
    make_plot(stats_cont(), jitter = FALSE, xlab = "IQ", ylab = "Chess skill")
  }, width = 800, height = 550)

  output$stats_cont <- renderUI({
    make_stats_html(stats_cont(), labels = c("Full population", "Chess grandmasters"))
  })

  output$plot_snp <- renderPlot({
    make_plot(stats_snp(), jitter = TRUE,
              xlab = "SNP1 dosage (0/1/2 coded alleles)",
              ylab = "SNP2 dosage (0/1/2 coded alleles)")
  }, width = 800, height = 550)

  output$stats_snp <- renderUI({
    stats <- stats_snp()
    make_stats_html(stats, labels = c(
      "Population",
      sprintf("Selected sample (n cases = %d, n controls = %d)", stats$n_cases, stats$n_controls)
    ))
  })

  output$r2_snp <- renderUI({
    HTML(make_r2_html(input$K_snp, input$P_snp, input$r2_1, input$r2_2))
  })
}

shinyApp(ui, server)
