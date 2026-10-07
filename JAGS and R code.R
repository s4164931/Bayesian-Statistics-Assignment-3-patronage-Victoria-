install.packages(c("runjags", "coda", "tidyverse"))
library(runjags)
library(coda)
library(ggplot2)
library(ggpubr)
source("plot_dist.R") # from https://github.com/rasmusab/distribution_diagrams (thankyou)
source("DBDA2E-utilities.R")
data = read.csv("~/Desktop/Victorian Metro Train Patronage (and Stony Point line).csv") #FIXME!!! GET A CORRECT DIRECTORY!!!

normal_weekday = ggplot(data = data, aes(`Pax_norm_weekday`)) + geom_density()
sch_hol_weekday = ggplot(data = data, aes(`Pax_sch_hol_weekday`)) + geom_density()
saturday = ggplot(data = data, aes(`Pax_Saturday`)) + geom_density()
sunday = ggplot(data = data, aes(`Pax_Sunday`)) + geom_density()
figure <- ggarrange(normal_weekday, sch_hol_weekday, saturday, sunday, nrow = 2, ncol = 2)
figure

normal_weekday_1 = ggplot(data = data, aes(`Pax_norm_weekday`)) + geom_histogram(bins = 100)
sch_hol_weekday_1 = ggplot(data = data, aes(`Pax_sch_hol_weekday`)) + geom_histogram(bins = 100)
saturday_1 = ggplot(data = data, aes(`Pax_Saturday`)) + geom_histogram(bins = 100)
sunday_1 = ggplot(data = data, aes(`Pax_Sunday`)) + geom_histogram(bins = 100)
figure_hist <- ggarrange(normal_weekday_1, sch_hol_weekday_1, saturday_1, sunday_1, nrow = 2, ncol = 2)
figure_hist

x = as.matrix(data[c("Pax_norm_weekday", "Pax_sch_hol_weekday", "Pax_Saturday", "Pax_Sunday")])
cor(x)

summary(x)
# s = stations - easy, also make a factor of this
# nrows
# ncols - 5?
# x - can be used as "j" for each column type for data conversion
# day types -> i -> factor this by using x col names, reference each one to the 'j' part
# data transformation?
station <- as.numeric(as.factor(data$Stop_name))
N_cols <- 4
N_rows <- dim(data)[1]
y <- data$Pax_annual
Data_List <- list(y = y, x = x, N_rows = N_rows, station = station, N_cols = N_cols)

dists$poisson$ddist_params <- 1.5
plot_dist(dists$poisson, labels = c(lambda = "λ"), scale = 5)
plot_dist(dists$normal, labels = c(mean = expression("μ"[βj]), right_sd = expression("τ"[βj])), scale = 4)
plot_dist(dists$normal, labels = c(mean = expression("M"[j]), right_sd = expression("S"[j])), scale = 4)
plot_dist(dists$gamma, labels = c(params = "A, B"), scale = 4)

model_text <- "
data {
  for (i in 1:N_rows) {
    log_y[i] <- log(y[i])
  }
  for (j in 1:N_cols) {
    for (i in 1:N_rows) {
      log_x[i, j] <- log(x[i, j])
    }
  }
  # prior settings
  mean_priors[1] <- 7
  mean_priors[2] <- -0.15 # FIXME!!! 6.85 or -0.15? (negatives first, if this fails, try the other)
  mean_priors[3] <- -0.5
  mean_priors[4] <- -0.8
  variance_priors[1] <- 0.25
  variance_priors[2] <- 16
  variance_priors[3] <- 4
  variance_priors[4] <- 4
  A_Gamma_prior[1] <- 1
  A_Gamma_prior[2] <- 2
  A_Gamma_prior[3] <- 2
  A_Gamma_prior[4] <- 2
  B_Gamma_prior[1] <- 1
  B_Gamma_prior[2] <- 0.1
  B_Gamma_prior[3] <- 0.4
  B_Gamma_prior[4] <- 0.4
}
model { #??? How do I do heirachial level???
  for (j in 1:N_rows) {
    log_y[j] ~ dpois(log_lambda[j])
    log_lambda[j] <- log_B[1, station[j]] + sum(log_B[station[j], 2:N_cols] * log_x[j, 2:N_cols]) # CAN JAGS ACCEPT THIS???
  }
  for (i in 1:N_cols) {
    log_B[i] ~ dnorm(log_mean_B[i], log_var_B[i])
    log_mean_B[i] ~ dnorm(mean_priors[i], variance_priors[i])
    log_var_B[i] ~ dnorm(A_Gamma_prior[i], B_Gamma_prior[i])
  }
  B[1:N_cols] <- exp(log_B[1:N_cols])
  mean_B[1:N_cols] <- exp(log_mean_B[1:N_cols])
  var_B[1:N_cols] <- exp(log_var_B[1:N_cols])
}"

writeLines(model_text, con = "station_patronage_model.txt")

failed.jags('model')

parameters = c("var_B", "mean_B", "B", "log_var_B", "log_mean_B", "log_B")
adaption_steps = 500
burnin_steps = 100
num_chains = 3
thinning_steps = 2
sample_steps = 1000
saved_steps = sample_steps * thinning_steps
start_time = proc.time()
Station_Model <- run.jags(method = "parallel", model = "station_patronage_model.txt", monitor = parameters, data = Data_List, n.chains = num_chains, adapt = adaption_steps, burnin = burnin_steps, thin = thinning_steps, sample = sample_steps, plots = FALSE)
end_time = proc.time()
time_taken = end_time - start_time
time_taken
model_list = as.mcmc.list(Station_Model)
model_list
# How does station patronage vary across day types 
# (normal weekday, school holiday weekday, Saturday, Sunday), and how much does 
# this pattern differ across Melbourne train stations? - across stations?