# Collider Bias Explorer

An interactive Shiny app for building intuition about **collider bias**
(a.k.a. Berkson's paradox / selection bias): conditioning a sample on a
variable that is itself caused by (or correlated with) two other variables
can induce a spurious correlation between those two variables, even when
none exists in the full population.

**Live app:** https://89mz98-emil-uffelmann.shinyapps.io/collider_bias/

## What's in the app

The app has two tabs, each simulating data live as you move the sliders.
Both show the same core plot: a scatter of two variables X and Y, a black
dashed regression line fit to the full population, and a red regression
line fit to a selected/ascertained sample. Watch how the red line's slope
and correlation diverge from the black one as selection gets stronger.

### Continuous variables

X and Y are two continuous traits with an adjustable true correlation. A
third variable S — the *selection variable* — is generated to have chosen
correlations with X and with Y (e.g. X = intelligence, Y = creativity, S =
becoming a professional chess player). Restricting the sample to
individuals with high S can create a correlation between X and Y that
doesn't exist in the population, and can even flip its sign, purely
because S depends on both.

Adjustable:
- True correlation between X and Y
- Correlation between X and S, and between Y and S
- Selection strength (top % selected on S)
- Sample size

### Genetic (SNPs)

X (SNP1) and Y (SNP2) are two genetic variants (dosages of 0/1/2 copies of
a coded allele) that are independent in the population (no LD). Each
explains an adjustable share of the variance in the liability of a disease
(e.g. Alzheimer's), under the standard liability-threshold model. The
disease is the collider: individuals with liability above a threshold
(set by the population prevalence, K) are cases. Ascertaining a
case-control study with a case fraction (sample prevalence, P) that
differs from K induces a spurious SNP1-SNP2 correlation in the selected
sample, even though the SNPs are independent in the population — and its
sign depends on whether the two SNPs' coded alleles move disease risk in
the same or in opposite directions.

Adjustable:
- Coded allele frequency and effect direction (risk-increasing /
  risk-decreasing) for each SNP
- Variance each SNP explains in disease liability (R², liability scale)
- Population disease prevalence (K) and sample/study disease prevalence
  (P); P = K means no ascertainment, P = 1 is a cases-only sample, P = 0
  is a controls-only sample
- Sample size

The panel below the plot also shows, for each SNP, how its true
liability-scale R² translates into an inflated *observed*-scale R² under
the current ascertainment, using the liability-scale conversion of
[Lee et al. 2012, *Genetic Epidemiology*](https://doi.org/10.1002/gepi.21614).

## Running locally

```r
install.packages(c("shiny", "ggplot2"))
shiny::runApp("app.R")
```

## Deploying

The app is deployed to shinyapps.io via `rsconnect`:

```r
rsconnect::deployApp(appName = "collider_bias")
```
