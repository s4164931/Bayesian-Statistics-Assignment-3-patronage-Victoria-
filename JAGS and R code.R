install.packages(c("runjags", "coda", "tidyverse"))
library(runjags)
library(coda)
library(ggplot2)
library(ggpubr)
source("~/Downloads/plot_dist.R") # from https://github.com/rasmusab/distribution_diagrams (thankyou)
source("~/Downloads/DBDA2E-utilities.R")
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
figure_hist <- ggarrange(normal_weekday, sch_hol_weekday, saturday, sunday, nrow = 2, ncol = 2)
figure_hist

# How does station patronage vary across day types 
# (normal weekday, school holiday weekday, Saturday, Sunday), and how much does 
# this pattern differ across Melbourne train stations? - across stations?