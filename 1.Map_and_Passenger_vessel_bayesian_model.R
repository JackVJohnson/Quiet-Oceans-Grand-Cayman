#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: Map and Passenger vessel bayesian models 
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: April 2023
# Last edited: September 2023

# clear workspace
rm(list=ls())

# libraries

library(tidyverse)
library(fishualize)
library(patchwork)
library(sf)
library(ggmap)
library(ggsn)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggcorrplot)
library(rstan)
library(brms)
library(ggrepel)

#########################################################################################
################################## working directories ##################################

data_wd <- "C:/Users/jackv/OneDrive - Central Caribbean Marine Institute, Inc/Projects/Quiet_Oceans/Data files"

shapefile_wd <- "C:/Users/jackv/OneDrive - Central Caribbean Marine Institute, Inc/Shapefiles/Cayman Islands"

output_wd <-  "C:/Users/jackv/OneDrive - Central Caribbean Marine Institute, Inc/Projects/Quiet_Oceans/Output directory"

my_theme <- theme_classic() +
  theme(axis.title.x = element_text(size = 20, color = "black"), axis.title.y = element_text(size = 20, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=16, color="black"), axis.text.x = element_text(size=16, color="black", angle = 0, vjust = .0)) +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold")) 

##########################################################################################
##################################### read in files ######################################

cruise_data <- read.csv(file=file.path(data_wd, "cruise_ship_data.csv"))
df <- read.csv(file=file.path(data_wd, "ccmi_quiet_ocean_survey.csv"))
df <- df[,1:17]
summary(df)

# clean data

df <- na.omit(df)
df$site <- gsub("fish_pot", "fishpot", df$site)
df$site <- gsub("bonnies_arch_deep", "bonnies_arch", df$site)
df$site <- gsub("bonnies_arch_shallow", "bonnies_arch", df$site)
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
df$site <- gsub("bonnies_arch_deep", "bonnies_arch", df$site)
df$site <- gsub("bonnies_arch_shallow", "bonnies_arch", df$site)
df <- subset(df, !(site == "bonnies_arch"))

# remove bonnies arch because it's not near the cruise ships 


table(df$Date)

tapply(df$latin_name, df$trophic_guild, FUN = unique)
length(unique(df$latin_name))

tapply(df$latin_family, df$trophic_guild, FUN = unique)
length(unique(df$latin_family))

# summarise by dates

dates_df <- df %>%
  group_by(site) %>%
  summarise(n=max(transect))

sites_df <- read.csv(file=file.path(data_wd, "site_coords.csv"))
sites_df <- subset(sites_df, !Island == "Little_Cayman")

sites_df$site <- gsub("fish_pot", "fishpot", sites_df$site)

sites_df <- sites_df[!duplicated(sites_df), ]

sites_df <- left_join(dates_df, sites_df, by = "site")
sites_df$code <- c("DF", "ER", "FP", "WF")

# read in shape file 

setwd(shapefile_wd)

gc_map <- st_read("cym_admbnda_adm0_2020.shp")
gc_map
gc_map <- st_transform(gc_map, crs = 4326)

p_map1 <- ggplot() +
  geom_sf(data=gc_map, fill = "#A9DA3F") +
  coord_sf(xlim=c(-81.478113,-81.031176), ylim=c(19.20,19.438963), expand = F) +
  geom_point(data = sites_df, aes(x = Longitude, y = Latitude, color = code), size = 3) +
  geom_text_repel(aes(x=Longitude, y=Latitude, label=code, color = code), data=sites_df, hjust=1, vjust = -0, size = 5) +
  scale_color_viridis_d(option = "A") +
  annotation_scale() +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm") ,
                         style = north_arrow_fancy_orienteering) +
  theme_classic() +
  xlab("") + ylab("") +
  #theme(panel.background = element_rect(fill = "#39a786")) +
  #geom_rect(aes(xmin = -81.472893, xmax = -81.378822, ymin = 19.250968, ymax = 19.413184), color = "black", fill = NA, linewidth = 1) +
  annotate("text", x=-78, y=15, label = "Caribbean Sea", size = 8) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black", angle = 0, vjust = .0)) +
  theme(panel.background = element_rect(fill = "#5ECFFA")) +
  theme(legend.position = "none")
p_map1

ocean <- ne_download(type = "ocean", category = "physical", scale = "medium")
ocean <- st_as_sf(ocean)
p_map2 <- ggplot() +
  geom_sf(data = ocean, fill = "#5ECFFA") +
  coord_sf(xlim=c(-92.478113,-62.031176), ylim=c(8.8,27), expand = T) +
  scale_x_continuous(breaks = c(-85, -75, -65)) + 
  annotation_scale() +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm") ,
                         style = north_arrow_fancy_orienteering) +
  theme_classic() +
  xlab("") + ylab("") +
  theme(panel.background = element_rect(fill = "#A9DA3F")) +
  geom_rect(aes(xmin = -82.25, xmax = -80.2, ymin = 18.25, ymax = 20.25), color = "black", fill = NA, linewidth = 1) +
  annotate("text", x=-78, y=15, label = "Caribbean Sea", size = 8) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black", angle = 0, vjust = .0)) 
p_map2

# combine maps and export

png(file=file.path(output_wd, "Quiet_oceans_survey_sites_FIGURE.png"), height = 3000, width = 6000, res = 350)
(p_map1 + p_map2 + plot_layout(ncol=2, widths = c(2, 1))) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 30))
dev.off()


# summary plot of cruise ships and passengers 
cruise_data$date <- paste(cruise_data$Year, cruise_data$Month, 1, sep = "-")

cruise_data$date <- as.Date(cruise_data$date, format ="%Y-%B-%d")


p_cruise <- ggplot(cruise_data, aes(as.factor(date), Passenger_vessel_calls)) +
  geom_point() +
  geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  ylab("Number of cruise ships") + xlab("Time period") +
  scale_x_discrete(breaks = c("2018-02-01", "2020-04-01", "2022-06-01"), labels=c("February 2018", "April 2020", "June 2022")) +
  my_theme

p4 <- ggplot() +
  theme_void()

png(file=file.path(output_wd, "Survey_sites_and_cruise_plotspace_FIGURE.png"), height = 3500, width = 5000, res = 350)
(p_map1 + p_map2 + plot_layout(widths = c(2, 1))) / (p_cruise + p4) + plot_layout(heights = unit(c(10,1), c('cm', 'null'))) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 30, face = "bold"))
dev.off()




####################################################################################
################################## Bayesian model ##################################

df$Year <- substr(df$Date, 1, 4)
df$date_index <- paste(df$Year, df$Month, sep = "_")
df$date_index <- gsub(" ", "", df$date_index)

df$count <- as.integer(df$count)

df$biomass_g <- as.numeric(df$biomass_g)

overall_df <- df %>%
  group_by(site, Month, Year, date_index) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


overall_df$year_month <- ym(paste(overall_df$Year, overall_df$Month, sep = "-"))
str(overall_df$year_month)

# create a vector of dates
date_vector <- overall_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
overall_df$rs_dates <- rescaled_dates


# responses

hist(overall_df$mean_richness)
hist(log1p(overall_df$mean_richness))
shapiro.test(log1p(overall_df$mean_richness))

hist(overall_df$mean_biomass)
hist(log1p(overall_df$mean_biomass))
shapiro.test(log1p(overall_df$mean_biomass))

hist(overall_df$mean_abundance)
hist(log1p(overall_df$mean_abundance))
shapiro.test(log1p(overall_df$mean_abundance))

# correlation matrix

cor_matrix <- subset(overall_df, select=c(mean_richness, mean_biomass, mean_abundance, rs_dates))

cor_matrix <- cor(cor_matrix, method = "spearman")
p_cor_mat <- cor_pmat(cor_matrix, method = "spearman")

png(file=file.path(output_wd, "Correlation_matrix_SUPP.png"), height = 2000, width = 2000, res = 350)
ggcorrplot(cor_matrix, hc.order = TRUE,
           type = "lower", lab = TRUE) 
dev.off()

# model multiple cores

bform2 <- bf(mvbind(round(mean_richness), round(mean_biomass), round(mean_abundance)) ~ rs_dates + (1|site))
options(mc.cores=parallel::detectCores())

fit_multi_time <- brm(bform2, 
                      data = overall_df,
                      family = negbinomial(),
                      warmup = 1500,
                      iter = 3000, 
                      chains = 4,
                      cores = 4, 
                      control = list(adapt_delta = 0.99))

prior_summary(fit_multi_time)

summary(fit_multi_time)
pp_check(fit_multi_time, resp='roundmeanrichness')
pp_check(fit_multi_time, resp='roundmeanbiomass')
pp_check(fit_multi_time, resp='roundmeanabundance')

mcmc_plot(fit_multi_time)
mcmc_plot(fit_multi_time, type = "trace")
mcmc_plot(fit_multi_time, var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

models <- "models"
saveRDS(file = file.path(output_wd, models, 'Allfish_TIME.rds'), fit_multi_time)
#fit_multi_time <- readRDS(file=file.path(output_wd, models, "Allfish_TIME.rds"))

# coefficient plot

coeffs_timeall <- posterior_summary(fit_multi_time, digits = 0.6, probs=c(0.95, 0.05,0.8, 0.2,0.5))
coeffs_timeall <- as.data.frame(coeffs_timeall[4:6,])
coeffs_timeall <- rownames_to_column(coeffs_timeall)
coeffs_timeall$predictor <- "Time period"

p1 <- ggplot(coeffs_timeall, aes(Estimate, predictor, color = rowname, fill = rowname)) +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +         
  geom_errorbarh(aes(xmax = Q95, xmin =Q5), size = 1.3, height = 0, color = "light grey", position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmax = Q80, xmin =Q20), size = 2, height = 0, color = "Dark grey", position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width=0.5)) +
  #scale_y_discrete(labels=c("Depth", "Turbidity", "Distance", "DHW*Distance", "DHW")) +
  scale_color_viridis_d(breaks=c("b_roundmeanrichness_rs_dates","b_roundmeanbiomass_rs_dates",   "b_roundmeanabundance_rs_dates"), labels=c("Mean species richness", "Mean biomass", "Mean abundance"), option = "A", begin = 0.35, end = 0.8) +
  scale_fill_viridis_d(breaks=c("b_roundmeanrichness_rs_dates","b_roundmeanbiomass_rs_dates",   "b_roundmeanabundance_rs_dates"), labels=c("Mean species richness", "Mean biomass", "Mean abundance"), option = "A", begin = 0.35, end = 0.8) +
  labs(color = "Response", fill = "Response", x = "Model coefficient", y = "") +
  theme_classic() +
  my_theme
p1

# export figure 
png(file=file.path(output_wd, "Overall_mod_coeffs_FIGURE.png"), height = 2000, width = 3000, res = 350)
p1 
dev.off()

# compare to a gaussian model 

#bform1 <- bf(mvbind(log1p(mean_richness), log1p(mean_biomass), log1p(mean_abundance)) ~ Passenger_vessel_calls +  (1|site))

#fit_multi_gaussian <- brm(bform1, 
#                          data = overall_df,
#                         family = gaussian(),
#                        warmup = 1500,
#                       iter = 3000, 
#                      chains = 4,
#                     cores = 4, 
#                    control = list(adapt_delta = 0.99))

#summary(fit_multi_gaussian)

#loo(fit_multi, fit_multi_gaussian)


########################################################################################
####################################### Herbivores #####################################

herbivore_df <- df %>%
  filter(trophic_guild == "Herbivore") %>%
  group_by(date_index, site, Month, Year) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


herbivore_df$year_month <- ym(paste(herbivore_df$Year, herbivore_df$Month, sep = "-"))

# create a vector of dates
date_vector <- herbivore_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
herbivore_df$rs_dates <- rescaled_dates

# model

fit_herbivore_time <- brm(bform2, 
                          data = herbivore_df,
                          family = negbinomial(),
                          warmup = 1500,
                          iter = 3000, 
                          chains = 4,
                          cores = 4, 
                          control = list(adapt_delta = 0.99))


summary(fit_herbivore_time)
pp_check(fit_herbivore_time, resp='roundmeanrichness')
pp_check(fit_herbivore_time, resp='roundmeanbiomass')
pp_check(fit_herbivore_time, resp='roundmeanabundance')

mcmc_plot(fit_herbivore_time)
mcmc_plot(fit_herbivore_time, type = "trace")
mcmc_plot(fit_herbivore_time, var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))


########################################################################################
##################################### Macrocarnivore ###################################

macrocarnivore_df <- df %>%
  filter(trophic_guild == "Macrocarnivore") %>%
  group_by(date_index, site, Month, Year) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


macrocarnivore_df$year_month <- ym(paste(macrocarnivore_df$Year, macrocarnivore_df$Month, sep = "-"))

# create a vector of dates
date_vector <- macrocarnivore_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
macrocarnivore_df$rs_dates <- rescaled_dates


# model

fit_macrocarnivore_time <- brm(bform2, 
                               data = macrocarnivore_df,
                               family = negbinomial(),
                               warmup = 1500,
                               iter = 3000, 
                               chains = 4,
                               cores = 4, 
                               control = list(adapt_delta = 0.99))


summary(fit_macrocarnivore_time)
pp_check(fit_macrocarnivore_time, resp='roundmeanrichness')
pp_check(fit_macrocarnivore_time, resp='roundmeanbiomass')
pp_check(fit_macrocarnivore_time, resp='roundmeanabundance')

mcmc_plot(fit_macrocarnivore_time)
mcmc_plot(fit_macrocarnivore_time, type = "trace")
mcmc_plot(fit_macrocarnivore_time, var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))


########################################################################################
###################################### Invertivores ####################################

invertivore_df <- df %>%
  filter(trophic_guild == "Invertivore") %>%
  group_by(date_index, site, Month, Year) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


invertivore_df$year_month <- ym(paste(invertivore_df$Year, invertivore_df$Month, sep = "-"))

# create a vector of dates
date_vector <- invertivore_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
invertivore_df$rs_dates <- rescaled_dates


# model

fit_invertivore_time <- brm(bform2, 
                            data = invertivore_df,
                            family = negbinomial(),
                            warmup = 1500,
                            iter = 3000, 
                            chains = 4,
                            cores = 4, 
                            control = list(adapt_delta = 0.99))


summary(fit_invertivore_time)
pp_check(fit_invertivore_time, resp='roundmeanrichness')
pp_check(fit_invertivore_time, resp='roundmeanbiomass')
pp_check(fit_invertivore_time, resp='roundmeanabundance')

mcmc_plot(fit_invertivore_time)
mcmc_plot(fit_invertivore_time, type = "trace")
mcmc_plot(fit_invertivore_time, var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))


########################################################################################
###################################### Omnivores ####################################

omnivore_df <- df %>%
  filter(trophic_guild == "Omnivore") %>%
  group_by(date_index, site, Month, Year) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))

omnivore_df$year_month <- ym(paste(omnivore_df$Year, omnivore_df$Month, sep = "-"))

# create a vector of dates
date_vector <- omnivore_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
omnivore_df$rs_dates <- rescaled_dates


# model

fit_omnivore_time <- brm(bform2, 
                         data = omnivore_df,
                         family = negbinomial(),
                         warmup = 1500,
                         iter = 3000, 
                         chains = 4,
                         cores = 4, 
                         control = list(adapt_delta = 0.99))


summary(fit_omnivore_time)
pp_check(fit_omnivore_time, resp='roundmeanrichness')
pp_check(fit_omnivore_time, resp='roundmeanbiomass')
pp_check(fit_omnivore_time, resp='roundmeanabundance')

mcmc_plot(fit_omnivore_time)
mcmc_plot(fit_omnivore_time, type = "trace")
mcmc_plot(fit_omnivore_time, var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))


########################################################################################
###################################### Planktivores ####################################

planktivore_df <- df %>%
  filter(trophic_guild == "Planktivore") %>%
  group_by(date_index, site, Month, Year) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))


planktivore_df$year_month <- ym(paste(planktivore_df$Year, planktivore_df$Month, sep = "-"))

# create a vector of dates
date_vector <- planktivore_df$year_month
# Create a unique set of dates
unique_dates <- unique(date_vector)
# Sort the unique dates
sorted_dates <- sort(unique_dates)
# Assign sequential numbers to the sorted dates
rescaled_dates <- match(date_vector, sorted_dates)
# Check the resulting rescaled dates
print(rescaled_dates)
planktivore_df$rs_dates <- rescaled_dates


# model

fit_planktivore_time <- brm(bform2, 
                            data = planktivore_df,
                            family = negbinomial(),
                            warmup = 1500,
                            iter = 3000, 
                            chains = 4,
                            cores = 4, 
                            control = list(adapt_delta = 0.99))


summary(fit_planktivore_time)
pp_check(fit_planktivore_time, resp='roundmeanrichness')
pp_check(fit_planktivore_time, resp='roundmeanbiomass')
pp_check(fit_planktivore_time, resp='roundmeanabundance')

mcmc_plot(fit_planktivore_time)
mcmc_plot(fit_planktivore_time, type = "trace")
mcmc_plot(fit_planktivore_time, var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

# save models 

saveRDS(file = file.path(output_wd, models, 'Herbivore_TIME.rds'), fit_herbivore_time)
saveRDS(file = file.path(output_wd, models, 'Macrocarnivore_TIME.rds'), fit_macrocarnivore_time)
saveRDS(file = file.path(output_wd, models, 'Invertivore_TIME.rds'), fit_invertivore_time)
saveRDS(file = file.path(output_wd, models, 'Omnivore_TIME.rds'), fit_omnivore_time)
saveRDS(file = file.path(output_wd, models, 'Planktivore_TIME.rds'), fit_planktivore_time)

summary(fit_herbivore_time)
summary(fit_invertivore_time)
summary(fit_macrocarnivore_time)
summary(fit_omnivore_time)
summary(fit_planktivore_time)
