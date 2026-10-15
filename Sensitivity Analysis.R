# Prior sensitivity analysis for the hierarchical patronage model
#
# This script fits the same likelihood and hierarchical structure under three
# prior settings. Only the population-level parameters are monitored because
# the purpose is to test whether the main day-type conclusions depend on the
# priors, not to recreate every station-level diagnostic plot.

required_packages <- c("runjags", "coda", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the required package(s) before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(runjags)
  library(coda)
  library(ggplot2)
})

set.seed(2269)

# Find the project directory whether the script is run with Source or Rscript.
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) == 1) {
  script_path <- sub("^--file=", "", script_argument)
  script_path <- gsub("~\\+~", " ", script_path)
  project_dir <- dirname(normalizePath(script_path))
} else {
  project_dir <- normalizePath(getwd())
}

data_path <- file.path(
  project_dir,
  "Victorian Metro Train Patronage (and Stony Point line).csv"
)
output_dir <- file.path(project_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_path)) {
  stop("Data file not found at: ", data_path)
}

patronage_data <- read.csv(data_path, check.names = FALSE)
day_columns <- c(
  "Pax_norm_weekday",
  "Pax_sch_hol_weekday",
  "Pax_Saturday",
  "Pax_Sunday"
)

if (!all(day_columns %in% names(patronage_data))) {
  stop("The dataset does not contain all four required patronage columns.")
}

patronage_matrix <- as.matrix(patronage_data[day_columns])
if (anyNA(patronage_matrix) || any(patronage_matrix < 0)) {
  stop("Patronage values must be non-missing and non-negative.")
}

N_station <- nrow(patronage_data)
N_obs <- N_station * length(day_columns)
y <- as.integer(round(c(t(patronage_matrix))))
station <- rep(seq_len(N_station), each = length(day_columns))
school_holiday <- rep(c(0, 1, 0, 0), times = N_station)
saturday <- rep(c(0, 0, 1, 0), times = N_station)
sunday <- rep(c(0, 0, 0, 1), times = N_station)

# The hierarchy-variance priors stay fixed so this analysis isolates the
# effect of changing the priors on the population means.
A_gamma <- c(1, 2, 2, 2)
B_gamma <- c(1.0, 0.1, 0.4, 0.4)

prior_settings <- list(
  Baseline = list(
    mean = c(7.00, -0.15, -0.50, -0.80),
    precision = c(0.25, 16.00, 4.00, 4.00)
  ),
  Vague = list(
    mean = c(7.00, 0.00, 0.00, 0.00),
    precision = c(0.01, 0.01, 0.01, 0.01)
  ),
  Sceptical = list(
    mean = c(7.00, 0.00, 0.00, 0.00),
    precision = c(0.25, 16.00, 16.00, 16.00)
  )
)

prior_table <- do.call(rbind, lapply(names(prior_settings), function(scenario) {
  setting <- prior_settings[[scenario]]
  data.frame(
    Scenario = scenario,
    Parameter = c("mu[1]", "mu[2]", "mu[3]", "mu[4]"),
    PriorMean = setting$mean,
    PriorSD = 1 / sqrt(setting$precision),
    row.names = NULL
  )
}))
write.csv(
  prior_table,
  file.path(output_dir, "sensitivity_prior_settings.csv"),
  row.names = FALSE
)

model_text <- "
model {
  for (n in 1:N_obs) {
    y[n] ~ dpois(lambda[n])
    log(lambda[n]) <- beta[station[n], 1] +
      beta[station[n], 2] * school_holiday[n] +
      beta[station[n], 3] * saturday[n] +
      beta[station[n], 4] * sunday[n]
  }

  for (s in 1:N_station) {
    for (k in 1:4) {
      beta[s, k] ~ dnorm(mu[k], tau[k])
    }
  }

  for (k in 1:4) {
    mu[k] ~ dnorm(mean_prior[k], mean_precision[k])
    tau[k] ~ dgamma(A_gamma[k], B_gamma[k])
    sigma[k] <- 1 / sqrt(tau[k])
  }

  rate_normal_weekday <- exp(mu[1])
  rate_school_holiday <- exp(mu[2])
  rate_saturday <- exp(mu[3])
  rate_sunday <- exp(mu[4])
}
"

model_path <- file.path(output_dir, "sensitivity_model.txt")
writeLines(model_text, model_path)

# Use SENSITIVITY_FAST_TEST=1 for a short code check. The default is the final
# analysis run used for the report.
fast_test <- identical(Sys.getenv("SENSITIVITY_FAST_TEST"), "1")
adapt_steps <- if (fast_test) 300 else 2000
burnin_steps <- if (fast_test) 500 else 5000
sample_steps <- if (fast_test) 1000 else 10000
thin_steps <- if (fast_test) 2 else 10

monitored_parameters <- c(
  "mu",
  "sigma",
  "rate_normal_weekday",
  "rate_school_holiday",
  "rate_saturday",
  "rate_sunday"
)

rng_names <- c(
  "base::Wichmann-Hill",
  "base::Marsaglia-Multicarry",
  "base::Super-Duper",
  "base::Mersenne-Twister"
)

# Data-based initial values improve numerical stability but do not change the
# priors or posterior distribution.
initial_mu <- c(
  log(mean(patronage_matrix[, 1])),
  log(mean(patronage_matrix[, 2]) / mean(patronage_matrix[, 1])),
  log(mean(patronage_matrix[, 3]) / mean(patronage_matrix[, 1])),
  log(mean(patronage_matrix[, 4]) / mean(patronage_matrix[, 1]))
)

fit_scenario <- function(scenario, setting) {
  message("Fitting sensitivity scenario: ", scenario)

  jags_data <- list(
    N_obs = N_obs,
    N_station = N_station,
    y = y,
    station = station,
    school_holiday = school_holiday,
    saturday = saturday,
    sunday = sunday,
    mean_prior = setting$mean,
    mean_precision = setting$precision,
    A_gamma = A_gamma,
    B_gamma = B_gamma
  )

  initial_values <- lapply(seq_len(4), function(chain) {
    list(
      mu = rnorm(4, initial_mu, 0.02),
      tau = rgamma(4, shape = A_gamma, rate = B_gamma),
      .RNG.name = rng_names[chain],
      .RNG.seed = 2269 + 100 * match(scenario, names(prior_settings)) + chain
    )
  })

  fitted <- run.jags(
    method = "parallel",
    model = model_path,
    monitor = monitored_parameters,
    data = jags_data,
    inits = initial_values,
    n.chains = 4,
    adapt = adapt_steps,
    burnin = burnin_steps,
    sample = sample_steps,
    thin = thin_steps,
    plots = FALSE,
    summarise = FALSE
  )

  as.mcmc.list(fitted)
}

model_fits <- Map(fit_scenario, names(prior_settings), prior_settings)
names(model_fits) <- names(prior_settings)

rate_names <- c(
  rate_normal_weekday = "Normal weekday",
  rate_school_holiday = "School-holiday weekday",
  rate_saturday = "Saturday",
  rate_sunday = "Sunday"
)

summarise_fit <- function(model_list, scenario) {
  summary_stats <- summary(model_list)
  posterior <- data.frame(
    Parameter = rownames(summary_stats$statistics),
    Mean = summary_stats$statistics[, "Mean"],
    SD = summary_stats$statistics[, "SD"],
    Lower95 = summary_stats$quantiles[, "2.5%"],
    Median = summary_stats$quantiles[, "50%"],
    Upper95 = summary_stats$quantiles[, "97.5%"],
    row.names = NULL
  )

  rate_summary <- posterior[posterior$Parameter %in% names(rate_names), ]
  rate_summary$Scenario <- scenario
  rate_summary$DayType <- unname(rate_names[rate_summary$Parameter])
  rate_summary$PercentChange <- 100 * (rate_summary$Mean - 1)

  ess <- effectiveSize(model_list)
  rhat <- gelman.diag(
    model_list,
    multivariate = FALSE,
    autoburnin = FALSE
  )$psrf[, "Point est."]

  diagnostics <- data.frame(
    Scenario = scenario,
    Parameter = names(ess),
    ESS = as.numeric(ess),
    Rhat = as.numeric(rhat[names(ess)]),
    row.names = NULL
  )

  list(rates = rate_summary, diagnostics = diagnostics)
}

scenario_summaries <- Map(summarise_fit, model_fits, names(model_fits))
rate_comparison <- do.call(rbind, lapply(scenario_summaries, `[[`, "rates"))
diagnostic_details <- do.call(
  rbind,
  lapply(scenario_summaries, `[[`, "diagnostics")
)

baseline_means <- setNames(
  rate_comparison$Mean[rate_comparison$Scenario == "Baseline"],
  rate_comparison$Parameter[rate_comparison$Scenario == "Baseline"]
)
rate_comparison$DifferenceFromBaseline <-
  rate_comparison$Mean - baseline_means[rate_comparison$Parameter]
rate_comparison$RelativeDifferencePercent <-
  100 * rate_comparison$DifferenceFromBaseline /
  baseline_means[rate_comparison$Parameter]

diagnostic_summary <- do.call(rbind, lapply(
  split(diagnostic_details, diagnostic_details$Scenario),
  function(values) {
    data.frame(
      Scenario = values$Scenario[1],
      MaximumRhat = max(values$Rhat, na.rm = TRUE),
      MinimumESS = min(values$ESS, na.rm = TRUE),
      AllRhatBelow1.01 = all(values$Rhat < 1.01, na.rm = TRUE),
      AllESSAbove400 = all(values$ESS > 400, na.rm = TRUE),
      row.names = NULL
    )
  }
))

write.csv(
  rate_comparison,
  file.path(output_dir, "sensitivity_rate_ratio_comparison.csv"),
  row.names = FALSE
)
write.csv(
  diagnostic_details,
  file.path(output_dir, "sensitivity_diagnostic_details.csv"),
  row.names = FALSE
)
write.csv(
  diagnostic_summary,
  file.path(output_dir, "sensitivity_diagnostic_summary.csv"),
  row.names = FALSE
)

rate_comparison$Scenario <- factor(
  rate_comparison$Scenario,
  levels = c("Sceptical", "Vague", "Baseline")
)
rate_comparison$DayType <- factor(
  rate_comparison$DayType,
  levels = c("Normal weekday", "School-holiday weekday", "Saturday", "Sunday")
)

sensitivity_plot <- ggplot(
  rate_comparison,
  aes(x = Mean, y = Scenario, colour = Scenario)
) +
  geom_errorbar(
    aes(xmin = Lower95, xmax = Upper95),
    orientation = "y",
    width = 0.18,
    linewidth = 0.8
  ) +
  geom_point(size = 2.7) +
  facet_wrap(~DayType, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = c(
    Baseline = "#0072B2",
    Vague = "#009E73",
    Sceptical = "#D55E00"
  )) +
  labs(
    x = "Posterior rate ratio with 95% credible interval",
    y = NULL,
    title = "Prior sensitivity of day-type effects",
    subtitle = "Intervals overlap closely across all four prior settings"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "sensitivity_rate_ratios.png"),
  sensitivity_plot,
  width = 8,
  height = 7,
  dpi = 300
)

# Save the compact population-level chains for reproducibility.
saveRDS(
  model_fits,
  file.path(output_dir, "sensitivity_population_mcmc.rds")
)

cat("\nPrior settings:\n")
print(prior_table, row.names = FALSE)
cat("\nPosterior rate-ratio comparison:\n")
print(
  rate_comparison[
    , c(
      "Scenario", "DayType", "Mean", "Lower95", "Upper95",
      "RelativeDifferencePercent"
    )
  ],
  row.names = FALSE
)
cat("\nSensitivity-run diagnostics:\n")
print(diagnostic_summary, row.names = FALSE)
cat("\nFiles saved in: ", output_dir, "\n", sep = "")
