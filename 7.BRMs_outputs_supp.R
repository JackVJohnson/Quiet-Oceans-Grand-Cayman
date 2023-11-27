#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: traceplots and ppchecks for review 
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: November 2023
# Last edited: November 2023

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

models <- "models"
m1 <- readRDS(file=file.path(output_wd, models, "Allfish_TIME.rds"))
a <- pp_check(m1, resp='roundmeanrichness')
b <- pp_check(m1, resp='roundmeanbiomass')
c <- pp_check(m1, resp='roundmeanabundance')
d <- mcmc_plot(m1, type = "trace", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

mcmc_plot(m1)

e <- mcmc_plot(m1, type ="dens", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

png(file=file.path(output_wd, "brms_out_ALLFISH.png"), height = 2000, width = 3500, res = 350)
d/(a+b+c)/e +plot_annotation(
  title = 'All fish model') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()

# herbivore

rm(m1,a,b,c,d,e)

m1 <- readRDS(file=file.path(output_wd, models, "HERBIVORE_TIME.rds"))
a <- pp_check(m1, resp='roundmeanrichness')
b <- pp_check(m1, resp='roundmeanbiomass')
c <- pp_check(m1, resp='roundmeanabundance')
d <- mcmc_plot(m1, type = "trace", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

mcmc_plot(m1)

e <- mcmc_plot(m1, type ="dens", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

png(file=file.path(output_wd, "brms_out_HERBIVORE.png"), height = 2000, width = 3500, res = 350)
d/(a+b+c)/e +plot_annotation(
  title = 'Herbivore model') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()

# macrocarnivore

rm(m1,a,b,c,d,e)

m1 <- readRDS(file=file.path(output_wd, models, "Macrocarnivore_TIME.rds"))
a <- pp_check(m1, resp='roundmeanrichness')
b <- pp_check(m1, resp='roundmeanbiomass')
c <- pp_check(m1, resp='roundmeanabundance')
d <- mcmc_plot(m1, type = "trace", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

mcmc_plot(m1)

e <- mcmc_plot(m1, type ="dens", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

png(file=file.path(output_wd, "brms_out_Macrocarnivore.png"), height = 2000, width = 3500, res = 350)
d/(a+b+c)/e +plot_annotation(
  title = 'Macrocarnivore model') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()

# invertivore

rm(m1,a,b,c,d,e)

m1 <- readRDS(file=file.path(output_wd, models, "Invertivore_TIME.rds"))
a <- pp_check(m1, resp='roundmeanrichness')
b <- pp_check(m1, resp='roundmeanbiomass')
c <- pp_check(m1, resp='roundmeanabundance')
d <- mcmc_plot(m1, type = "trace", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

mcmc_plot(m1)

e <- mcmc_plot(m1, type ="dens", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

png(file=file.path(output_wd, "brms_out_Invertivore.png"), height = 2000, width = 3500, res = 350)
d/(a+b+c)/e +plot_annotation(
  title = 'Invertivore model') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()

# omnivore 

rm(m1,a,b,c,d,e)

m1 <- readRDS(file=file.path(output_wd, models, "Omnivore_TIME.rds"))
a <- pp_check(m1, resp='roundmeanrichness')
b <- pp_check(m1, resp='roundmeanbiomass')
c <- pp_check(m1, resp='roundmeanabundance')
d <- mcmc_plot(m1, type = "trace", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

mcmc_plot(m1)

e <- mcmc_plot(m1, type ="dens", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

png(file=file.path(output_wd, "brms_out_Omnivore.png"), height = 2000, width = 3500, res = 350)
d/(a+b+c)/e +plot_annotation(
  title = 'Omnivore model') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()


# omnivore 

rm(m1,a,b,c,d,e)

m1 <- readRDS(file=file.path(output_wd, models, "Planktivore_TIME.rds"))
a <- pp_check(m1, resp='roundmeanrichness')
b <- pp_check(m1, resp='roundmeanbiomass')
c <- pp_check(m1, resp='roundmeanabundance')
d <- mcmc_plot(m1, type = "trace", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

mcmc_plot(m1)

e <- mcmc_plot(m1, type ="dens", var = c("b_roundmeanrichness_rs_dates", "b_roundmeanbiomass_rs_dates", "b_roundmeanabundance_rs_dates"))

png(file=file.path(output_wd, "brms_out_Planktivore.png"), height = 2000, width = 3500, res = 350)
d/(a+b+c)/e +plot_annotation(
  title = 'Planktivore model') & theme(plot.title = element_text(size = 20, face = "bold"))
dev.off()
