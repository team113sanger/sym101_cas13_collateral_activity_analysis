library(tidyverse)
library(here)
library(broom)
library(patchwork)

# Single-cell dose analysis: within a single guide, do cells carrying more of
# that guide (or expressing more Cas13d) show more ribosomal protein depletion?
# Guide identity enters every model as a fixed effect, so the dose slope is
# estimated strictly within guide. mRuby cells give the negative control slope.

cells <- read_csv(here("output", "single_cell_dose.csv"), show_col_types = FALSE) |>
  mutate(
    log_guide_umi = log1p(assigned_guide_umi),
    log_cas13 = log1p(cas13_expr),
    log_total_counts = log(total_counts),
    condition = fct_relevel(condition, "non_targeting"),
    across(c(guide, phase, cell_line), as_factor)
  ) |>
  drop_na(assigned_guide_umi, rp_fraction)

cat(str_glue("cells: {nrow(cells)}"), "\n")

# guide dose differs ~2-fold between conditions, so slopes are only ever
# compared within condition, never by contrasting dose across conditions
cells |>
  summarise(
    n = n(),
    median_guide_umi = median(assigned_guide_umi),
    mean_rp = mean(rp_fraction),
    mean_cas13 = mean(cas13_expr),
    .by = c(cell_line, condition)
  ) |>
  print()

# Cas13_expr is a raw count, so it rises with library complexity (rho 0.34-0.42
# with n_genes_detected). Without that covariate a detection gradient reads as a
# dose effect, so it is included; the sparser model is kept for comparison.
covariates <- "log_total_counts + n_genes_detected + mito_fraction + phase"
covariates_no_complexity <- "log_total_counts + mito_fraction + phase"

fit_dose <- function(data, dose, covariate_set = covariates, include_guide = TRUE) {
  terms <- c(dose, covariate_set, if (include_guide && n_distinct(data[["guide"]]) > 1) "guide")
  formula <- reformulate(terms, response = "rp_fraction")
  lm(formula, data = data) |>
    tidy() |>
    filter(term == dose)
}

# within-condition dose slopes, per cell line
slopes <- cells |>
  nest(.by = c(cell_line, condition)) |>
  expand_grid(
    dose = c("log_guide_umi", "log_cas13"),
    model = c("with_complexity", "without_complexity")
  ) |>
  mutate(
    result = pmap(list(data, dose, model), \(d, dd, mm) {
      fit_dose(d, dd, if (mm == "with_complexity") covariates else covariates_no_complexity)
    })
  ) |>
  select(-data) |>
  unnest(result) |>
  select(cell_line, condition, dose, model, estimate, std.error, statistic, p.value) |>
  arrange(dose, model, cell_line, condition)

cat("\n=== within-guide dose slopes on RP fraction ===\n")
print(slopes, n = Inf)
write_csv(slopes, here("output", "single_cell_dose_slopes.csv"))

# formal test that the slope differs between targeting and control cells
interaction_tests <- cells |>
  nest(.by = cell_line) |>
  expand_grid(dose = c("log_guide_umi", "log_cas13")) |>
  mutate(
    result = map2(data, dose, \(d, dd) {
      formula <- reformulate(
        c(str_glue("{dd} * condition"), covariates, "guide"),
        response = "rp_fraction"
      )
      lm(formula, data = d) |>
        tidy() |>
        filter(str_detect(term, str_glue("^{dd}:condition")))
    })
  ) |>
  select(-data) |>
  unnest(result) |>
  select(cell_line, dose, term, estimate, std.error, statistic, p.value)

cat("\n=== dose x condition interaction (targeting vs mRuby slope) ===\n")
print(interaction_tests, n = Inf)
write_csv(interaction_tests, here("output", "single_cell_dose_interactions.csv"))

# does collateral need both the guide and the enzyme?
joint <- cells |>
  filter(condition == "targeting") |>
  nest(.by = cell_line) |>
  mutate(
    result = map(data, \(d) {
      formula <- reformulate(
        c("log_guide_umi * log_cas13", covariates, "guide"),
        response = "rp_fraction"
      )
      lm(formula, data = d) |>
        tidy() |>
        filter(str_detect(term, "log_guide_umi|log_cas13"))
    })
  ) |>
  select(-data) |>
  unnest(result)

cat("\n=== guide x Cas13d interaction, targeting cells only ===\n")
print(joint, n = Inf)
write_csv(joint, here("output", "single_cell_guide_by_cas13.csv"))

plot_dose_bins <- function(dose_column, dose_label, n_bins = 10) {
  cells |>
    mutate(
      bin = ntile(.data[[dose_column]], n_bins),
      .by = c(cell_line, condition)
    ) |>
    summarise(
      mean_rp = mean(rp_fraction),
      se = sd(rp_fraction) / sqrt(n()),
      dose = median(.data[[dose_column]]),
      n = n(),
      .by = c(cell_line, condition, bin)
    ) |>
    ggplot(aes(dose, mean_rp, colour = condition)) +
    geom_errorbar(aes(ymin = mean_rp - se, ymax = mean_rp + se), width = 0, linewidth = 0.3) +
    geom_point(size = 1.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    facet_wrap(vars(cell_line), scales = "free") +
    scale_colour_manual(
      values = c(non_targeting = "#1B9E77", targeting = "#B2182B"),
      labels = c(non_targeting = "mRuby control", targeting = "targeting")
    ) +
    labs(
      title = str_glue("Ribosomal protein fraction vs {dose_label}"),
      subtitle = str_c(
        "decile bins within cell line and condition; mean +/- SE. ",
        "Unadjusted - see single_cell_dose_slopes.csv for covariate-adjusted slopes"
      ),
      x = dose_label, y = "ribosomal protein fraction of counts", colour = NULL
    ) +
    theme_bw(base_size = 9)
}

ggsave(
  here("output", "single_cell_dose_response.png"),
  plot_dose_bins("log_guide_umi", "log1p(assigned guide UMIs)") /
    plot_dose_bins("log_cas13", "log1p(Cas13d expression)"),
  width = 10, height = 7.5, dpi = 300
)
