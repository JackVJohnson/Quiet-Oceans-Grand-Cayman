
#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: parrotfish brms model 
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: April 2023
# Last edited: September 2023

# clear workspace
rm(list=ls())

# libraries

library(tidyverse)
library(fishualize)
library(patchwork)
library(ggcorrplot)
library(ggrepel)
library(viridis)
library(brms)
library(tidybayes)


my_theme <- theme_classic() +
  theme(axis.title.x = element_text(size = 20, color = "black"), axis.title.y = element_text(size = 20, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=16, color="black"), axis.text.x = element_text(size=16, color="black", angle = 0, vjust = .0)) +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold")) 

#########################################################################################
################################## working directories ##################################

data_wd <- "C:/Users/jackv/OneDrive - Central Caribbean Marine Institute, Inc/Projects/Quiet_Oceans/Data files"

output_wd <-  "C:/Users/jackv/OneDrive - Central Caribbean Marine Institute, Inc/Projects/Quiet_Oceans/Output directory"



##########################################################################################
##################################### read in files ######################################

df <- read.csv(file=file.path(data_wd, "ccmi_quiet_ocean_survey.csv"))
df <- df[,1:17]
summary(df)

# clean data

df <- na.omit(df)
df$site <- gsub("fish_pot", "fishpot", df$site)
df$site <- gsub("bonnies_arch_deep", "bonnies_arch", df$site)
df$site <- gsub("bonnies_arch_shallow", "bonnies_arch", df$site)

# remove bonnies arch 
df <- subset(df, !(site == "bonnies_arch"))

df$site <- gsub(" ","", df$site)

df <- subset(df, !Island == "little_cayman")
df$Month <- trimws(df$Month)

# british dates not yank

df$Date <- gsub("6/13/22", "13/06/2022", df$Date)
df$Date <- gsub("6/14/22", "13/06/2022", df$Date)
df$Date <- gsub("12/02/21", "02/12/2021", df$Date)
df$Date <- gsub("12/03/21", "03/12/2021", df$Date)

df$Date <- as.Date(df$Date, format = "%d/%m/%Y")

# remove 2018 because of different survey method

df <- subset(df, !(Date == "2018-07-12"|Date == "2018-07-13"))

table(df$Date)
df$Year <- substr(df$Date, 1, 4)
df$date_index <- paste(df$Year, df$Month, sep = "_")
df$date_index <- gsub(" ", "", df$date_index)

df$count <- as.integer(df$count)

df$biomass_g <- as.numeric(df$biomass_g)


# filter the fish names that end in "terminal" or "initial"
subset_df <- df[grepl("terminal$|initial$", df$common_name), ]

summarised_df <-  subset_df %>%
  group_by(site, Month, Year, date_index) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


summarised_df$year_month <- ym(paste(summarised_df$Year, summarised_df$Month, sep = "-"))
str(summarised_df$year_month)

# create a vector of dates
date_vector <- summarised_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
summarised_df$rs_dates <- rescaled_dates

ggplot(summarised_df, aes(rs_dates, mean_biomass)) +
  geom_point() +
  geom_smooth(method=glm)

hist(summarised_df$mean_biomass)

m1 <- glm(mean_biomass~rs_dates, data = summarised_df, family = quasipoisson())
summary(m1)

options(mc.cores=parallel::detectCores())

fit1 <- brm(round(mean_biomass) ~ rs_dates,
            data = summarised_df,
            family = negbinomial(),
            warmup = 1500,
            iter = 3000, 
            chains = 4,
            cores = 4, 
            control = list(adapt_delta = 0.99))


summary(fit1) #mean=0.06,u95=0.09, l95=0.04
pp_check(fit1)
mcmc_plot(fit1, type = "trace")
mcmc_plot(fit1, var = "b_rs_dates")

models <- "models"
saveRDS(file = file.path(output_wd, models, 'All_parrotfish_TIME.rds'), fit1)
fit1 <- readRDS(file = file.path(output_wd, models, 'All_parrotfish_TIME.rds'))

model_fit <-summarised_df %>%
  add_predicted_draws(fit1)

p_allscar  <-  ggplot(model_fit,aes(x = rs_dates, y = mean_biomass)) +  
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                  alpha = 0.5, colour = "black") +
  geom_point(data = summarised_df, colour = "#b73779", size = 3) +   # raw data     
  scale_fill_brewer(palette = "Greys") +
  labs(x="Time period", y="Biomass (g)", title="", subtitle="All Parrotfish") +
  scale_x_continuous(breaks = c(2,7,12), labels = c("2020-09", "2021-08", "2022-04")) +
  my_theme +
  theme(legend.position = "none")

p_allscar

# initial phase only 

initial_df <- df[grepl("initial$", df$common_name), ]

summarised_initial_df <- initial_df %>%
  group_by(site, Month, Year, date_index) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


summarised_initial_df$year_month <- ym(paste(summarised_initial_df$Year, summarised_initial_df$Month, sep = "-"))
str(summarised_initial_df$year_month)

# create a vector of dates
date_vector <- summarised_initial_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
summarised_initial_df$rs_dates <- rescaled_dates

ggplot(summarised_initial_df, aes(rs_dates, mean_biomass)) +
  geom_point() +
  geom_smooth(method=glm)

hist(summarised_initial_df$mean_biomass)

m1 <- glm(mean_biomass~rs_dates, data = summarised_df, family = quasipoisson())
summary(m1)

options(mc.cores=parallel::detectCores())

fit2 <- brm(round(mean_biomass) ~ rs_dates + (1|site),
            data = summarised_initial_df,
            family = negbinomial(),
            warmup = 1500,
            iter = 3000, 
            chains = 4,
            cores = 4, 
            control = list(adapt_delta = 0.99))


summary(fit2) #mean=0.08,u95=0.11, l95=0.05
pp_check(fit2)
mcmc_plot(fit2, type = "trace")
mcmc_plot(fit2, var = "b_rs_dates")

models <- "models"
saveRDS(file = file.path(output_wd, models, 'Juvenile_parrotfish_TIME.rds'), fit2)
fit2 <- readRDS(file = file.path(output_wd, models, 'Juvenile_parrotfish_TIME.rds'))

model_fit2 <-summarised_initial_df %>%
  add_predicted_draws(fit2)

p_initialscar  <-  ggplot(model_fit2,aes(x = rs_dates, y = mean_biomass)) +  
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                  alpha = 0.5, colour = "black") +
  geom_point(data = summarised_df, colour = "#b73779", size = 3) +   # raw data     
  scale_fill_brewer(palette = "Greys") +
  labs(x="Time period", y="", title="", subtitle="Initial Phase Parrotfish") +
  scale_x_continuous(breaks = c(2,7,12), labels = c("2020-09", "2021-08", "2022-04")) +
  my_theme +
  theme(legend.position = "none")
p_initialscar


png(file=file.path(output_wd, "Parrotfish_bayes_regression_FIGURE.png"), height = 2000, width = 4000, res = 350)
p_allscar + p_initialscar + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 30, face = "bold"))
dev.off()

tapply(summarised_initial_df$mean_biomass, summarised_initial_df$year_month, mean)
tapply(summarised_initial_df$mean_biomass, summarised_initial_df$year_month, sd)

######################################################################################
############################ extract trace and ppc ###################################

m1 <- fit1

a <- pp_check(m1)
b <- mcmc_plot(m1, type = "trace", var = c("b_rs_dates"))
c <- mcmc_plot(m1, type = "dens", var = "b_rs_dates")

m2 <- fit2

x <- pp_check(m2)
y <- mcmc_plot(m2, type = "trace", var = c("b_rs_dates"))
z <- mcmc_plot(m2, type = "dens", var = "b_rs_dates")


png(file=file.path(output_wd, "brms_out_Parrotfish.png"), height = 2000, width = 3500, res = 350)
(a+b+c)/(x+y+z) + plot_annotation(
  title = 'Parrotfish model and Initial Parrotfish') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()
