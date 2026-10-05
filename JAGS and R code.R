install.packages(c("runjags", "coda", "tidyverse"))
library(runjags)
library(coda)
library(ggplot2)
library(ggpubr)
source("plot_dist.R") # from https://github.com/rasmusab/distribution_diagrams (thankyou)
source("DBDA2E-utilities.R")
data = read.csv("Victorian Metro Train Patronage (and Stony Point line).csv")

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
dists$poisson$ddist_params <- 1.5
plot_dist(dists$poisson, labels = c(lambda = "λ"), scale = 5)
plot_dist(dists$normal, labels = c(mean = expression("μ"[βj]), right_sd = expression("τ"[βj])), scale = 4)
plot_dist(dists$normal, labels = c(mean = expression("M"[j]), right_sd = expression("S"[j])), scale = 4)
plot_dist(dists$gamma, labels = c(params = "A, B"), scale = 4)
# How does station patronage vary across day types 
# (normal weekday, school holiday weekday, Saturday, Sunday), and how much does 
# this pattern differ across Melbourne train stations? - across stations?