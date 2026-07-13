# Collider Bias Explorer
#
# An interactive demo of collider bias (Berkson's paradox): selecting a
# sample based on two variables jointly can induce a spurious correlation
# between them, even when no such correlation exists (or a different one
# exists) in the full population. Selecting on only one of the two variables
# does not introduce this bias.
#
# Two tabs:
#   - Continuous variables: two correlated continuous traits X and Y
#   - Genetic (SNPs): two correlated SNP dosages (0/1/2), i.e. two variants
#     in LD, illustrating how case/cohort ascertainment can distort apparent
#     SNP-SNP correlation

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

simulate_continuous <- function(n, rho, seed) {
  simulate_latent(n, rho, seed)
}

# Two SNP dosages (0/1/2) with correlation (LD) rho, via a Gaussian copula:
# correlated liabilities -> uniform ranks -> binomial genotype at each MAF.
simulate_snps <- function(n, rho, maf1, maf2, seed) {
  latent <- simulate_latent(n, rho, seed)
  u1 <- pnorm(latent$x)
  u2 <- pnorm(latent$y)
  g1 <- qbinom(u1, 2, maf1)
  g2 <- qbinom(u2, 2, maf2)
  data.frame(x = g1, y = g2)
}

# Classify points by selection status:
#   - "Not selected"       : x below its threshold
#   - "Selected: X only"   : x above threshold, y below its threshold
#   - "Selected: X and Y"  : both x and y above their thresholds
# select_pct = percent of the population selected on each variable
# (marginally), so smaller values = more stringent (stronger) selection.
classify_selection <- function(df, select_pct) {
  thr_x <- quantile(df$x, probs = 1 - select_pct / 100, type = 1, names = FALSE)
  thr_y <- quantile(df$y, probs = 1 - select_pct / 100, type = 1, names = FALSE)
  df$group <- "Not selected"
  df$group[df$x >= thr_x & df$y < thr_y] <- "Selected: X only"
  df$group[df$x >= thr_x & df$y >= thr_y] <- "Selected: X and Y"
  df$group <- factor(df$group,
                      levels = c("Not selected", "Selected: X only", "Selected: X and Y"))
  attr(df, "thr_x") <- thr_x
  attr(df, "thr_y") <- thr_y
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
          "Selected: X only" = "#2C7FB8",
          "Selected: X and Y" = "#D7191C")

# Classify + fit all three regression lines once; shared by the plot and the
# text summary so the (identical) computation isn't duplicated.
compute_stats <- function(df, select_pct) {
  df <- classify_selection(df, select_pct)
  x_data    <- df[df$group %in% c("Selected: X only", "Selected: X and Y"), , drop = FALSE]
  both_data <- df[df$group == "Selected: X and Y", , drop = FALSE]
  list(
    df = df,
    pop_fit  = fit_line(df),
    x_fit    = fit_line(x_data),
    both_fit = fit_line(both_data)
  )
}

# Build the scatter + three regression lines. `jitter` is used for the SNP
# tab where x/y only take values 0/1/2.
make_plot <- function(stats, jitter = FALSE, xlab = "X", ylab = "Y") {
  df <- stats$df
  pop_fit <- stats$pop_fit
  x_fit <- stats$x_fit
  both_fit <- stats$both_fit

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
    (if (!is.na(x_fit$slope))
      geom_abline(intercept = x_fit$intercept, slope = x_fit$slope,
                  colour = COLS[["Selected: X only"]], linewidth = 1.1)
     else NULL) +
    (if (!is.na(both_fit$slope))
      geom_abline(intercept = both_fit$intercept, slope = both_fit$slope,
                  colour = COLS[["Selected: X and Y"]], linewidth = 1.1)
     else NULL) +
    scale_colour_manual(values = COLS, name = "Selection status", drop = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")

  p
}

# HTML summary of the three correlations, colour-matched to the plot legend,
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
    line("Selected on X only", stats$x_fit, COLS[["Selected: X only"]]),
    line("Selected on X and Y", stats$both_fit, COLS[["Selected: X and Y"]]),
    sep = "&nbsp;&nbsp;|&nbsp;&nbsp;"
  ))
}

explanation <- paste(
  "Collider bias (Berkson's paradox): X and Y are simulated with a fixed,",
  "known correlation in the full population (black dashed line). If a",
  "sample is then selected requiring BOTH variables to be high",
  "(red points/line), a spurious correlation can appear between X and Y",
  "even when little or none exists in the population - and it can even flip",
  "sign. Selecting on only ONE variable (blue points/line, X high",
  "regardless of Y) does not introduce this bias: its slope stays close to",
  "the population slope."
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
      sliderInput("select_pct_cont",
                  "Selection strength (% selected on each variable)",
                  min = 5, max = 50, value = 20, step = 1),
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
        "SNPs. The correlation slider mimics LD between the variants;",
        "selection mimics ascertaining a cohort (e.g. cases) based on a",
        "genetic risk score that depends on one or both SNPs."
      )),
      sliderInput("rho_snp", "LD (correlation) between SNP1 and SNP2",
                  min = -0.9, max = 0.9, value = 0, step = 0.05),
      sliderInput("maf1", "SNP1 risk allele frequency",
                  min = 0.05, max = 0.5, value = 0.3, step = 0.05),
      sliderInput("maf2", "SNP2 risk allele frequency",
                  min = 0.05, max = 0.5, value = 0.3, step = 0.05),
      sliderInput("select_pct_snp",
                  "Selection strength (% selected on each SNP)",
                  min = 5, max = 50, value = 20, step = 1),
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
    df <- simulate_continuous(input$n_cont, input$rho_cont, seed_cont())
    compute_stats(df, input$select_pct_cont)
  })

  stats_snp <- reactive({
    df <- simulate_snps(input$n_snp, input$rho_snp, input$maf1, input$maf2, seed_snp())
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
