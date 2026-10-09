# Hierarchical Bayesian model for Victorian station patronage
# Normal weekdays are the reference day type.
install.packages("glue")
library(glue)
source("~/Downloads/DBDA2E-utilities.R")
required_packages <- c("runjags", "coda")
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
})

set.seed(2269)

data_path <- "~/desktop/Victorian Metro Train Patronage (and Stony Point line).csv"
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_path)) {
  stop("Data file not found. Open the project in RStudio and run the script again.")
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

# Convert the four day-type columns to one observation per station and day type.
patronage_matrix <- as.matrix(patronage_data[day_columns])
if (anyNA(patronage_matrix) || any(patronage_matrix < 0)) {
  stop("Patronage values must be non-missing and non-negative.")
}

N_station <- nrow(patronage_data)
N_obs <- N_station * length(day_columns)
y <- as.integer(round(c(t(patronage_matrix))))
station <- rep(seq_len(N_station), each = length(day_columns))

# Indicator variables use a normal weekday as the reference category.
# Luke - how does the data interact with this???
school_holiday <- rep(c(0, 1, 0, 0), times = N_station)
saturday <- rep(c(0, 0, 1, 0), times = N_station)
sunday <- rep(c(0, 0, 0, 1), times = N_station)

day_type_summary <- data.frame(
  DayType = c("Normal weekday", "School-holiday weekday", "Saturday", "Sunday"),
  Mean = colMeans(patronage_matrix),
  Variance = apply(patronage_matrix, 2, var)
)
day_type_summary$VarianceToMean <- day_type_summary$Variance / day_type_summary$Mean
write.csv(
  day_type_summary,
  file.path(output_dir, "day_type_descriptive_summary.csv"),
  row.names = FALSE
)

# Prior means represent the expected normal-weekday log patronage and the
# expected log rate ratios for school holidays, Saturdays and Sundays.
mean_prior <- c(7.00, -0.15, -0.50, -0.80)
mean_precision <- c(0.25, 16.00, 4.00, 4.00)
A_gamma <- c(1, 2, 2, 2)
B_gamma <- c(1.0, 0.1, 0.4, 0.4)

jags_data <- list(
  N_obs = N_obs,
  N_station = N_station,
  y = y,
  station = station,
  school_holiday = school_holiday,
  saturday = saturday,
  sunday = sunday,
  mean_prior = mean_prior,
  mean_precision = mean_precision,
  A_gamma = A_gamma,
  B_gamma = B_gamma
)

# Each station has four well-observed patronage counts, so a centred hierarchy
# is used. This avoids the strong posterior dependence seen in the earlier
# log-variance formulation while allowing station-specific day-type effects.
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

  rate_school_holiday <- exp(mu[2])
  rate_saturday <- exp(mu[3])
  rate_sunday <- exp(mu[4])
}
"

model_path <- "station_patronage_model.txt"
writeLines(model_text, model_path)

# Set BAYES_FAST_TEST=1 for a short validation run before the final analysis.
fast_test <- identical(Sys.getenv("BAYES_FAST_TEST"), "1")
adapt_steps <- if (fast_test) 500 else 2000
burnin_steps <- if (fast_test) 1000 else 5000
sample_steps <- if (fast_test) 2000 else 10000
thin_steps <- 5

parameters <- c(
  "mu",
  "sigma",
  "tau",
  "rate_school_holiday",
  "rate_saturday",
  "rate_sunday",
  "beta",
  "lambda"
)

rng_names <- c(
  "base::Wichmann-Hill",
  "base::Marsaglia-Multicarry",
  "base::Super-Duper",
  "base::Mersenne-Twister"
)
initial_values <- lapply(seq_len(4), function(chain) {
  list(
    mu = rnorm(4, mean_prior, 0.05),
    tau = rgamma(4, shape = A_gamma, rate = B_gamma),
    .RNG.name = rng_names[chain],
    .RNG.seed = 2269 + chain
  )
})

start_time <- proc.time()
station_model <- run.jags(
  method = "parallel",
  model = model_path,
  monitor = parameters,
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
elapsed_time <- proc.time() - start_time
print(elapsed_time)

model_list <- as.mcmc.list(station_model)
diagMCMC(model_list, parName = "mu[1]")
diagMCMC(model_list, parName = "mu[2]")
diagMCMC(model_list, parName = "mu[3]")
diagMCMC(model_list, parName = "mu[4]")
diagMCMC(model_list, parName = "sigma[1]")
diagMCMC(model_list, parName = "sigma[2]")
diagMCMC(model_list, parName = "sigma[3]")
diagMCMC(model_list, parName = "sigma[4]")
diagMCMC(model_list, parName = "tau[1]")
diagMCMC(model_list, parName = "tau[2]")
diagMCMC(model_list, parName = "tau[3]")
diagMCMC(model_list, parName = "tau[4]")
diagMCMC(model_list, parName = "rate_school_holiday")
diagMCMC(model_list, parName = "rate_saturday")
diagMCMC(model_list, parName = "rate_sunday")
diagMCMC(model_list, parName = paste0("beta[195,1]"))
for (n in 1:N_obs) {
  diagMCMC(model_list, parName = paste0("lambda[",n,"]")) # 888
}
for (i in 1:nrow(patronage_data)) {
  for (k in 1:4) {
    diagMCMC(model_list, parName = paste0("beta[",i,",",k,"]")) # 888
  }
}
graphics.off()
dpois(y[4], lambda = 2250)
plotPost(model_list[,"lambda[1]"])
plotPost(model_list[,"beta[1,3]"])
plotPost(model_list[,"beta[1,4]"])

saveRDS(model_list, file.path(output_dir, "station_patronage_mcmc.rds"))
summary(model_list)

summary_stats <- summary(model_list)
posterior_summary <- data.frame(
  Parameter = rownames(summary_stats$statistics),
  Mean = summary_stats$statistics[, "Mean"],
  SD = summary_stats$statistics[, "SD"],
  Q2.5 = summary_stats$quantiles[, "2.5%"],
  Median = summary_stats$quantiles[, "50%"],
  Q97.5 = summary_stats$quantiles[, "97.5%"],
  row.names = NULL
)
write.csv(
  posterior_summary,
  file.path(output_dir, "posterior_parameter_summary.csv"),
  row.names = FALSE
)

effective_size <- effectiveSize(model_list)
rhat_matrix <- gelman.diag(model_list, multivariate = FALSE)$psrf
rhat <- rhat_matrix[, "Point est."]
diagnostic_summary <- data.frame(
  Parameter = names(effective_size),
  ESS = as.numeric(effective_size),
  Rhat = as.numeric(rhat[names(effective_size)]),
  row.names = NULL
)
write.csv(
  diagnostic_summary,
  file.path(output_dir, "mcmc_diagnostic_summary.csv"),
  row.names = FALSE
)

rate_parameters <- c(
  "rate_school_holiday",
  "rate_saturday",
  "rate_sunday"
)
rate_ratio_summary <- posterior_summary[
  posterior_summary$Parameter %in% rate_parameters,
]
write.csv(
  rate_ratio_summary,
  file.path(output_dir, "day_type_rate_ratios.csv"),
  row.names = FALSE
)

hyperparameter_names <- grep(
  "^(mu|sigma)\\[",
  colnames(model_list[[1]]),
  value = TRUE
)

png(
  file.path(output_dir, "trace_hyperparameters.png"),
  width = 1800,
  height = 1400,
  res = 180
)
par(mfrow = c(4, 2), mar = c(3, 3, 2, 1), oma = c(0, 0, 1, 0))
chain_colours <- c("black", "#D55E00", "#009E73", "#0072B2")
iteration_values <- as.numeric(time(model_list[[1]]))
for (parameter in hyperparameter_names) {
  parameter_range <- range(vapply(
    model_list,
    function(chain) range(chain[, parameter]),
    numeric(2)
  ))
  plot(
    iteration_values,
    model_list[[1]][, parameter],
    type = "l",
    col = chain_colours[1],
    ylim = parameter_range,
    xlab = "Iteration",
    ylab = "Value",
    main = parameter
  )
  for (chain in 2:length(model_list)) {
    lines(
      iteration_values,
      model_list[[chain]][, parameter],
      col = chain_colours[chain]
    )
  }
}
dev.off()

png(
  file.path(output_dir, "autocorrelation_hyperparameters.png"),
  width = 1800,
  height = 1400,
  res = 180
)
par(mfrow = c(4, 2), mar = c(3, 3, 2, 1), oma = c(0, 0, 1, 0))
for (parameter in hyperparameter_names) {
  chain_acf <- lapply(model_list, function(chain) {
    acf(as.numeric(chain[, parameter]), lag.max = 50, plot = FALSE)$acf
  })
  plot(
    0:50,
    chain_acf[[1]],
    type = "h",
    col = chain_colours[1],
    ylim = c(-0.2, 1),
    xlab = "Lag",
    ylab = "Autocorrelation",
    main = parameter
  )
  abline(h = 0, col = "grey60")
  for (chain in 2:length(model_list)) {
    lines(0:50, chain_acf[[chain]], col = chain_colours[chain])
  }
}
dev.off()

print(day_type_summary)
print(posterior_summary)
print(diagnostic_summary)
print(rate_ratio_summary)
