# Hierarchical Bayesian model for Victorian station patronage
# Normal weekdays are the reference day type.
source("~/Downloads/DBDA2E-utilities.R") # From Bayesian statistics module.
"*********************************************************************
  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
*********************************************************************" 
# this is where the source is from above.
required_packages <- c("runjags", "coda", "glue", "ggpubr", "tidyverse", "GGally")
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
  library(glue)
  library(ggpubr)
  library(ggplot2)
  library(GGally)
})

set.seed(2269)

diagMCMC_but_no_autoburnin <- function( codaObject , parName=varnames(codaObject)[1] ,
                                        saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("skyblue","black","royalblue","steelblue")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors, autoburnin = FALSE )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
}

join_datasets <- function (predicted_values, observed_values, filter_name) {
  weekday_patronage = predicted_values[predicted_values$patronage_type == filter_name, c(1,4)]
  colnames(weekday_patronage) <- c(filter_name, "type")
  colnames(observed_values) <- c("Normal Weekday Patronage", "School Holiday Weekday Patronage", "Saturday Patronage", "Sunday Patronage", "type")
  if (filter_name == "Normal Weekday Patronage") {
    i = 1
  } else if (filter_name == "School Holiday Weekday Patronage") {
    i = 2
  } else if (filter_name == "Saturday Patronage") {
    i = 3
  } else if (filter_name == "Sunday Patronage") {
    i = 4
  }
  weekday_observed_values <- observed_values[,c(i,5)]
  combined_dataframe <- rbind(weekday_observed_values, weekday_patronage)
  return(combined_dataframe)
}

compare_values <- function(station_1, patronage_1, station_2, patronage_2, poisson_value_dataframe) {
  comp_1 <- poisson_value_dataframe[(poisson_value_dataframe$station == station_1) & (poisson_value_dataframe$patronage_type == patronage_1),1]
  comp_2 <- poisson_value_dataframe[(poisson_value_dataframe$station == station_2) & (poisson_value_dataframe$patronage_type == patronage_2),1]
  print(comp_1)
  difference <- NULL
  difference <- data.frame(comp_1 - comp_2)
  mode_val <- as.numeric(names(sort(-table(difference)))[1]) #mode
  HDI <- HDIofMCMC(as.matrix(difference)) # this actually is the posterior from rpois, but I think it is still acceptable
  difference$comp1 <- comp_1
  difference$comp2 <- comp_2
  mode_comp_1 <- as.numeric(names(sort(-table(difference$comp1)))[1])
  mode_comp_2 <- as.numeric(names(sort(-table(difference$comp2)))[1])
  HDI_comp1 <- HDIofMCMC(as.matrix(difference$comp1))
  HDI_comp2 <- HDIofMCMC(as.matrix(difference$comp2))
  hist_compare <- ggplot(data = difference, aes(comp_1...comp_2)) + geom_histogram(color = "steelblue") + labs(x = "Difference of Patronage", title = paste0(station_1, " ", patronage_1, " vs ", station_2, " ", patronage_2), y = "Frequency") + theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'), axis.title = element_text(size = 10)) + geom_vline(xintercept=mode_val, color="orange", linetype="dashed")+ geom_vline(xintercept = 0, color = "green", linetype = "dashed") + geom_segment(aes(x = HDI[1], xend = HDI[2], y = 0, yend = 0), color = "red", size = 2)
  max_value <- max(ggplot_build(hist_compare)$data[[1]]$ymax)
  hist_compare <- hist_compare + annotate("text", x = mode_val, y = max_value + 50, label = paste0("mode = ", mode_val)) + annotate("text", x = HDI[1], y = 50, label = HDI[1], color = "white") + annotate("text", x = HDI[2], y = 50, label = HDI[2], color = "white") + annotate("text", x = (HDI[1] + HDI[2])/2, y = 100, label = "95% HDI", color = "white")
  
  scatter_compare <- ggplot(data = difference, aes(x = comp1, y = comp2)) + geom_point(color = "steelblue") + labs(x = paste0("Difference of Patronage for ", station_1, " ", patronage_1), title = paste0(station_1, " ", patronage_1, " vs ", station_2, " ", patronage_2), y = paste0("Difference of Patronage for ", station_2, " ", patronage_2)) + theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'), axis.title = element_text(size = 10)) + geom_abline(slope = 1, intercept = 0, color = "green")
  
  comp_val_1 <- ggplot(data = difference, aes(comp1)) + geom_histogram(color = "steelblue") + labs(x = "Patronage", title = paste0(station_1, " ", patronage_1), y = "Frequency") + theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'), axis.title = element_text(size = 10)) + geom_vline(xintercept=mode_comp_1, color="orange", linetype="dashed") + geom_segment(aes(x = HDI_comp1[1], xend = HDI_comp1[2], y = 0, yend = 0), color = "red", size = 2)
  max_value <- max(ggplot_build(comp_val_1)$data[[1]]$ymax)
  comp_val_1 <- comp_val_1 + annotate("text", x = mode_comp_1, y = max_value + 50, label = paste0("mode = ", mode_comp_1)) + annotate("text", x = HDI_comp1[1], y = 50, label = HDI_comp1[1], color = "white") + annotate("text", x = HDI_comp1[2], y = 50, label = HDI_comp1[2], color = "white") + annotate("text", x = (HDI_comp1[1] + HDI_comp1[2])/2, y = 100, label = "95% HDI", color = "white")
  
  comp_val_2 <- ggplot(data = difference, aes(comp2)) + geom_histogram(color = "steelblue") + labs(x = "Patronage", title = paste0(station_2, " ", patronage_2), y = "Frequency") + theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'), axis.title = element_text(size = 10)) + geom_vline(xintercept=mode_comp_2, color="orange", linetype="dashed") + geom_segment(aes(x = HDI_comp2[1], xend = HDI_comp2[2], y = 0, yend = 0), color = "red", size = 2)
  max_value <- max(ggplot_build(comp_val_2)$data[[1]]$ymax)
  comp_val_2 <- comp_val_2 + annotate("text", x = mode_comp_2, y = max_value + 50, label = paste0("mode = ", mode_comp_2)) + annotate("text", x = HDI_comp2[1], y = 50, label = HDI_comp2[1], color = "white") + annotate("text", x = HDI_comp2[2], y = 50, label = HDI_comp2[2], color = "white") + annotate("text", x = (HDI_comp2[1] + HDI_comp2[2])/2, y = 100, label = "95% HDI", color = "white")
  
  ggarrange(comp_val_1, hist_compare, scatter_compare, comp_val_2, nrow = 2, ncol = 2)
}

data_path <- "~/desktop/Victorian Metro Train Patronage (and Stony Point line).csv"
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_path)) {
  stop("Data file not found. Open the project in RStudio and run the script again.")
}

patronage_data <- read.csv(data_path, check.names = FALSE)

normal_weekday = ggplot(data = patronage_data, aes(`Pax_norm_weekday`)) + geom_histogram(bins = 100) + labs(x = "Normal Weekday Patronage", y = "Frequency", title = "Distribution of Normal Weekday Patronage") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), axis.title = element_text(size = 17.5), axis.text = element_text(size = 15))
sch_hol_weekday = ggplot(data = patronage_data, aes(`Pax_sch_hol_weekday`)) + geom_histogram(bins = 100) + labs(x = "School Holiday Weekday Patronage", y = "Frequency", title = "Distribution of School Holiday Weekday Patronage") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), axis.title = element_text(size = 17.5), axis.text = element_text(size = 15))
saturday = ggplot(data = patronage_data, aes(`Pax_Saturday`)) + geom_histogram(bins = 100) + labs(x = "Saturday Patronage", y = "Frequency", title = "Distribution of Saturday Patronage") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), axis.title = element_text(size = 17.5), axis.text = element_text(size = 15))
sunday = ggplot(data = patronage_data, aes(`Pax_Sunday`)) + geom_histogram(bins = 100) + labs(x = "Sunday Patronage", y = "Frequency", title = "Distribution of Sunday Patronage") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), axis.title = element_text(size = 17.5), axis.text = element_text(size = 15))
figure_hist <- ggarrange(normal_weekday, sch_hol_weekday, saturday, sunday, nrow = 2, ncol = 2)
figure_hist

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
options(scipen=999)
cor_plot <- ggpairs(patronage_matrix, columnLabels = c("Normal Weekday", "School Holiday Weekday", "Saturday", "Sunday")) + labs(x = "Patronage on x", y = "Patronage on y", title = "Correlations of Patronage Day Types") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), axis.title = element_text(size = 17.5), axis.text = element_text(size = 14))

cor_plot[2,1] <- cor_plot[2,1] + coord_cartesian(ylim = c(0, 57000))
cor_plot[4,2] <- cor_plot[4,2] + coord_cartesian(xlim = c(0, 57000))
cor_plot
N_station <- nrow(patronage_data)
N_obs <- N_station * length(day_columns)
y <- as.integer(round(c(t(patronage_matrix))))
station <- rep(seq_len(N_station), each = length(day_columns))

# Indicator variables use a normal weekday as the reference category.
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
mean_prior <- c(7.00, -0.15, -0.50, -0.80) #sensitivity
mean_precision <- c(0.25, 16.00, 4.00, 4.00) #sensitivity
A_gamma <- c(1, 2, 2, 2) #sensitivity
B_gamma <- c(1.0, 0.1, 0.4, 0.4) #sensitivity

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

  rate_normal_weekday <- exp(mu[1])
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
thin_steps <- 10

parameters <- c(
  "mu",
  "sigma",
  "tau",
  "rate_normal_weekday",
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
diagMCMC_but_no_autoburnin(model_list, parName = "mu[1]")
diagMCMC_but_no_autoburnin(model_list, parName = "mu[2]")
diagMCMC_but_no_autoburnin(model_list, parName = "mu[3]")
diagMCMC_but_no_autoburnin(model_list, parName = "mu[4]")
diagMCMC_but_no_autoburnin(model_list, parName = "sigma[1]")
diagMCMC_but_no_autoburnin(model_list, parName = "sigma[2]")
diagMCMC_but_no_autoburnin(model_list, parName = "sigma[3]")
diagMCMC_but_no_autoburnin(model_list, parName = "sigma[4]")
diagMCMC_but_no_autoburnin(model_list, parName = "tau[1]")
diagMCMC_but_no_autoburnin(model_list, parName = "tau[2]")
diagMCMC_but_no_autoburnin(model_list, parName = "tau[3]")
diagMCMC_but_no_autoburnin(model_list, parName = "tau[4]")
diagMCMC_but_no_autoburnin(model_list, parName = "rate_normal_weekday")
diagMCMC_but_no_autoburnin(model_list, parName = "rate_school_holiday")
diagMCMC_but_no_autoburnin(model_list, parName = "rate_saturday")
diagMCMC_but_no_autoburnin(model_list, parName = "rate_sunday")
diagMCMC_but_no_autoburnin(model_list, parName = paste0("beta[195,1]"))
for (n in 1:N_obs) {
  diagMCMC_but_no_autoburnin(model_list, parName = paste0("lambda[",n,"]")) # 888
}
for (i in 1:nrow(patronage_data)) {
  for (k in 1:4) {
    diagMCMC_but_no_autoburnin(model_list, parName = paste0("beta[",i,",",k,"]")) # 888
  }
}

graphics.off()
plotPost(model_list[,"lambda[4]"])
exp(plotPost(model_list[,"beta[1,1]"])[4] + plotPost(model_list[,"beta[1,4]"])[4])
plotPost(model_list[,"beta[1,4]"])

saveRDS(model_list, file.path(output_dir, "station_patronage_mcmc.rds"))
summary(model_list)

summary_stats <- summary(model_list)
for (i in 1:nrow(summary_stats$statistics)) {
  print(summary_stats$statistics[i,"Time-series SE"] / summary_stats$statistics[i,"SD"])
  name <- rownames(summary_stats$statistics)[i]
  if (summary_stats$statistics[i,"Time-series SE"] / summary_stats$statistics[i,"SD"] >= 0.05) {
    diagMCMC_but_no_autoburnin(model_list, parName = rownames(summary_stats$statistics)[i])
  }
}
for (i in 1:nrow(summary_stats$statistics)) {
  print(summary_stats$statistics[i,"Time-series SE"] / summary_stats$statistics[i,"Naive SE"])
  name <- rownames(summary_stats$statistics)[i]
  if (summary_stats$statistics[i,"Time-series SE"] / summary_stats$statistics[i,"Naive SE"] >= sqrt(20)) {
    diagMCMC_but_no_autoburnin(model_list, parName = rownames(summary_stats$statistics)[i])
  }
}
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
for (i in 1:length(effective_size)) {
  print(effective_size[i])
  if (effective_size[i] <= 20000) {
    diagMCMC_but_no_autoburnin(model_list, parName = names(effective_size[i]))
  }
}
graphics.off()
rhat_matrix <- gelman.diag(model_list, multivariate = FALSE)$psrf
for (i in 1:length(rhat_matrix[, "Upper C.I."])) {
  print(rhat_matrix[i, "Upper C.I."])
  if (rhat_matrix[i, "Upper C.I."] >= 1.01) {
    diagMCMC_but_no_autoburnin(model_list, parName = rownames(rhat_matrix)[i]) # shrink factor?
  }
}
rhat <- rhat_matrix[, "Point est."]
gelman.diag(model_list, multivariate = TRUE)$mpsrf # don't use this
nChain = length(model_list)
check = "placeholder"
for (i in 1:length(rhat_matrix[, "Upper C.I."])) {
  name = rownames(rhat_matrix)[i]
  for ( cIdx in 1:nChain ) {
    print(acf(model_list[,name][[cIdx]], lag.max = 1, plot=FALSE)[10]$acf)
    if (acf(model_list[,name][[cIdx]], lag.max = 1, plot=FALSE)[10]$acf >= 0.5 && check != name) {
      diagMCMC_but_no_autoburnin(model_list, parName = name)
      check = name
    }
  }
}

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

# Generating the posterior data and labelling it
a <- NULL
for (n in 1:N_obs) {
  a <- rbind(a, summarizePost(as.matrix(model_list, chains = TRUE)[,paste0("lambda[",n,"]")])[3])
  print(n) # showing how many done out of 888. Takes about 3 minutes.
}
a = as.vector(t(as.matrix(a)))
a
poisson_values <- rpois(888 * 10000, a) # 1: Normal Weekday, 2: School Holiday Weekday, 3: Saturday, 4: Sunday
poisson_value_dataframe <- data.frame(poisson_values)
poisson_value_dataframe$patronage_type <- "placeholder"
poisson_value_dataframe$patronage_type[seq(1, 8880000, by = 4)] <- "Normal Weekday Patronage"
poisson_value_dataframe$patronage_type[seq(2, 8880000, by = 4)] <- "School Holiday Weekday Patronage"
poisson_value_dataframe$patronage_type[seq(3, 8880000, by = 4)] <- "Saturday Patronage"
poisson_value_dataframe$patronage_type[seq(4, 8880000, by = 4)] <- "Sunday Patronage"

# station labelling
station_names <- patronage_data$Stop_name
poisson_value_dataframe$station <- "placeholder" # 1-4 = Sunbury, 5-8 = Digger's Rest ... 888, then repeat to start

for (i in 1:N_station) { # 222, each station is selected here
  indexes <- seq(1 + (4*i - 4) ,8880000, by = 888) # this checks over every combo (i starts at 1)
  indexes_plus_1 <- indexes + 1
  indexes_plus_2 <- indexes + 2
  indexes_plus_3 <- indexes + 3 # this is the by 4 part, to get each station.
  poisson_value_dataframe$station[indexes] <- station_names[i]
  poisson_value_dataframe$station[indexes_plus_1] <- station_names[i]
  poisson_value_dataframe$station[indexes_plus_2] <- station_names[i]
  poisson_value_dataframe$station[indexes_plus_3] <- station_names[i]
}

# combining actual distribution to observed distribution:
predicted_values <- poisson_value_dataframe[c((888 * 4117 + 1):(888 * 4118)),]
observed_values <- NULL
observed_values <- rbind(observed_values, patronage_data$Pax_norm_weekday)
observed_values <- rbind(observed_values, patronage_data$Pax_sch_hol_weekday)
observed_values <- rbind(observed_values, patronage_data$Pax_Saturday)
observed_values <- rbind(observed_values, patronage_data$Pax_Sunday)
observed_values <- data.frame(t(observed_values))
observed_values$type <- "Observed"
predicted_values$type <- "Predicted"

combined_dataframe <- join_datasets(predicted_values, observed_values, "Normal Weekday Patronage")
Weekday_Patronage = ggplot(combined_dataframe, aes(`Normal Weekday Patronage`, fill = type)) + geom_density(alpha = 0.2) + labs(x = "Patronage", title = "Weekday Predicted vs Observed data", fill = "Data") + theme(plot.title = element_text(hjust = 0.5, size = 40, face = 'bold'), axis.title = element_text(size = 30), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key.spacing.y = unit(0.3, 'cm'), axis.text = element_text(size = 18))

combined_dataframe <- join_datasets(predicted_values, observed_values, "School Holiday Weekday Patronage")
School_Weekday_Patronage = ggplot(combined_dataframe, aes(`School Holiday Weekday Patronage`, fill = type)) + geom_density(alpha = 0.2) + labs(x = "Patronage", title = "School Weekday Predicted vs Observed data", fill = "Data") + theme(plot.title = element_text(hjust = 0.5, size = 40, face = 'bold'), axis.title = element_text(size = 30), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key.spacing.y = unit(0.3, 'cm'), axis.text = element_text(size = 18))

combined_dataframe <- join_datasets(predicted_values, observed_values, "Saturday Patronage")
Saturday_Patronage = ggplot(combined_dataframe, aes(`Saturday Patronage`, fill = type)) + geom_density(alpha = 0.2) + labs(x = "Patronage", title = "Saturday Predicted vs Observed data", fill = "Data") + theme(plot.title = element_text(hjust = 0.5, size = 40, face = 'bold'), axis.title = element_text(size = 30), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key.spacing.y = unit(0.3, 'cm'), axis.text = element_text(size = 18))

combined_dataframe <- join_datasets(predicted_values, observed_values, "Sunday Patronage")
Sunday_Patronage = ggplot(combined_dataframe, aes(`Sunday Patronage`, fill = type)) + geom_density(alpha = 0.2) + labs(x = "Patronage", title = "Sunday Predicted vs Observed data", fill = "Data") + theme(plot.title = element_text(hjust = 0.5, size = 40, face = 'bold'), axis.title = element_text(size = 30), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.key.spacing.y = unit(0.3, 'cm'), axis.text = element_text(size = 18))

Weekday_Patronage # better 1 at a time
School_Weekday_Patronage
Saturday_Patronage
Sunday_Patronage
# Now, comparing station time to a different time and different stations.

compare_values("Stony Point", "Normal Weekday Patronage", "Hastings", "Normal Weekday Patronage", poisson_value_dataframe)

print(day_type_summary)
print(posterior_summary)
print(diagnostic_summary)
print(rate_ratio_summary)
