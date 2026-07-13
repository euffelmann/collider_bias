# Collider Bias Explorer
#
# An interactive demo of collider bias (Berkson's paradox): a sample gets
# selected based on a third variable S (the "collider") that is itself
# caused by / correlated with two other variables X and Y. Conditioning on
# S can induce a spurious correlation between X and Y, even when none
# exists in the full population.
#
# Classic example: X = intelligence, Y = creativity, S = becoming a
# professional chess player. Intelligence and creativity may be uncorrelated
# in the general population, but if both raise your odds of becoming a
# professional chess player, then among professional chess players
# intelligence and creativity can appear negatively correlated.
#
# Two tabs:
#   - Continuous variables: two (possibly correlated) continuous traits X, Y
#   - Genetic (SNPs): two correlated SNP dosages (0/1/2), i.e. two variants
#     in LD, with S representing e.g. a disease liability / risk score that
#     depends on both SNPs, illustrating how case ascertainment can distort
#     apparent SNP-SNP correlation

library(shiny)
library(ggplot2)

# ---------------------------------------------------------------------------
# Simulation helpers
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

# Two SNP dosages (0/1/2) with correlation (LD) rho, via a Gaussian copula:
# correlated liabilities -> uniform ranks -> binomial genotype at each MAF.
# The selection variable S is computed on the underlying liability scale
# (e.g. a disease risk score), then attached to the dosage data.
simulate_snps <- function(n, rho, maf1, maf2, rho_xs, rho_ys, seed) {
  latent <- simulate_latent(n, rho, seed)
  latent <- add_selection_var(latent, rho, rho_xs, rho_ys)
  u1 <- pnorm(latent$x)
  u2 <- pnorm(latent$y)
  g1 <- qbinom(u1, 2, maf1)
  g2 <- qbinom(u2, 2, maf2)
  data.frame(x = g1, y = g2, s = latent$s)
}

# Classify points by selection status based on the selection variable S:
# top select_pct% of S (marginally) are "Selected" (e.g. became a
# professional chess player / a case). Lower % = more stringent selection.
classify_selection <- function(df, select_pct) {
  thr_s <- quantile(df$s, probs = 1 - select_pct / 100, type = 1, names = FALSE)
  df$group <- factor(ifelse(df$s >= thr_s, "Selected", "Not selected"),
                      levels = c("Not selected", "Selected"))
  attr(df, "thr_s") <- thr_s
  df
}

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

COLS <- c("Not selected" = "grey70",
          "Selected" = "#D7191C")

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

# Build the scatter + two regression lines. `jitter` is used for the SNP
# tab where x/y only take values 0/1/2.
make_plot <- function(stats, jitter = FALSE, xlab = "X", ylab = "Y") {
  df <- stats$df
  pop_fit <- stats$pop_fit
  sel_fit <- stats$sel_fit

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
                  colour = COLS[["Selected"]], linewidth = 1.1)
     else NULL) +
    scale_colour_manual(values = COLS, name = "Selection status", drop = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")

  p
}

# HTML summary of the two correlations, colour-matched to the plot legend,
# shown below the plot (avoids clipping that a ggplot subtitle suffers at
# typical plot widths).
make_stats_html <- function(stats) {
  line <- function(label, fit, colour) {
    sprintf(
      '<span style="color:%s; font-weight:600;">%s</span>: r = %s (n = %d)',
      colour, label, fmt_r(fit), fit$n
    )
  }
  HTML(paste(
    line("Full population", stats$pop_fit, "black"),
    line("Selected sample", stats$sel_fit, COLS[["Selected"]]),
    sep = "&nbsp;&nbsp;|&nbsp;&nbsp;"
  ))
}

explanation <- paste(
  "Collider bias (Berkson's paradox): X and Y are simulated with a fixed,",
  "known correlation in the full population (black dashed line). A third",
  "variable S - the selection variable - is correlated with X and with Y",
  "by the amounts set below (e.g. X = intelligence, Y = creativity, S =",
  "becoming a professional chess player). If the sample is restricted to",
  "those with high S (red points/line), a spurious correlation can appear",
  "between X and Y even when little or none exists in the population - and",
  "it can even flip sign - purely because S depends on both."
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

continuous_tab <- tabPanel(
  "Continuous variables",
  sidebarLayout(
    sidebarPanel(
      helpText(explanation),
      sliderInput("rho_cont", "True correlation between X and Y",
                  min = -0.9, max = 0.9, value = 0, step = 0.05),
      sliderInput("rho_xs_cont", "Correlation between X and the selection variable",
                  min = -0.9, max = 0.9, value = 0.6, step = 0.05),
      sliderInput("rho_ys_cont", "Correlation between Y and the selection variable",
                  min = -0.9, max = 0.9, value = 0.6, step = 0.05),
      sliderInput("select_pct_cont",
                  "Selection strength (top % selected on the selection variable)",
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
      helpText(explanation),
      helpText(paste(
        "Here X and Y are dosages (0/1/2 copies of the risk allele) for two",
        "SNPs, and the selection variable represents e.g. a disease",
        "liability or polygenic risk score that depends on both SNPs.",
        "Ascertaining a cohort (e.g. cases) based on that score can distort",
        "the apparent SNP-SNP correlation (LD) even if the true LD is weak."
      )),
      sliderInput("rho_snp", "LD (correlation) between SNP1 and SNP2",
                  min = -0.9, max = 0.9, value = 0, step = 0.05),
      sliderInput("maf1", "SNP1 risk allele frequency",
                  min = 0.05, max = 0.5, value = 0.3, step = 0.05),
      sliderInput("maf2", "SNP2 risk allele frequency",
                  min = 0.05, max = 0.5, value = 0.3, step = 0.05),
      sliderInput("rho_xs_snp", "Correlation between SNP1 and the selection variable",
                  min = -0.9, max = 0.9, value = 0.6, step = 0.05),
      sliderInput("rho_ys_snp", "Correlation between SNP2 and the selection variable",
                  min = -0.9, max = 0.9, value = 0.6, step = 0.05),
      sliderInput("select_pct_snp",
                  "Selection strength (top % selected on the selection variable)",
                  min = 1, max = 50, value = 10, step = 1),
      helpText("Lower % = more stringent (stronger) selection."),
      sliderInput("n_snp", "Sample size", min = 500, max = 5000,
                  value = 2000, step = 500),
      actionButton("resample_snp", "Draw new random sample")
    ),
    mainPanel(
      plotOutput("plot_snp", height = "550px"),
      div(style = "margin-top: 10px; font-size: 14px;", htmlOutput("stats_snp"))
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

  stats_snp <- reactive({
    df <- simulate_snps(input$n_snp, input$rho_snp, input$maf1, input$maf2,
                         input$rho_xs_snp, input$rho_ys_snp, seed_snp())
    compute_stats(df, input$select_pct_snp)
  })

  output$plot_cont <- renderPlot({
    make_plot(stats_cont(), jitter = FALSE, xlab = "Variable X", ylab = "Variable Y")
  }, width = 800, height = 550)

  output$stats_cont <- renderUI({
    make_stats_html(stats_cont())
  })

  output$plot_snp <- renderPlot({
    make_plot(stats_snp(), jitter = TRUE,
              xlab = "SNP1 dosage (0/1/2 risk alleles)",
              ylab = "SNP2 dosage (0/1/2 risk alleles)")
  }, width = 800, height = 550)

  output$stats_snp <- renderUI({
    make_stats_html(stats_snp())
  })
}

shinyApp(ui, server)
