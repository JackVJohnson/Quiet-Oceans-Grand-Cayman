#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: analysis before - during - after
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: September 2023
# Last edited: 

# clear workspace
rm(list=ls())

# libraries

library(tidyverse)
library(fishualize)
library(patchwork)
library(FSA)
library(ggsignif)


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
df$site <- gsub(" ","", df$site)

df <- subset(df, !Island == "little_cayman")
df$Month <- trimws(df$Month)

# british dates not yank

df$Date <- gsub("6/13/22", "13/06/2022", df$Date)
df$Date <- gsub("6/14/22", "13/06/2022", df$Date)
df$Date <- gsub("12/02/21", "02/12/2021", df$Date)
df$Date <- gsub("12/03/21", "03/12/2021", df$Date)

df$Date <- as.Date(df$Date, format = "%d/%m/%Y")

table(df$site)

# remove bonnies arch because it's not near the cruise ships 

df <- subset(df, !(site == "bonnies_arch"))

# comparison before, during, after

df <- df %>%
  mutate(
    period = case_when(
      year(Date) <= 2018 ~ "before",
      between(Date, as.Date("2020-07-13"), as.Date("2022-02-07")) ~ "during",
      TRUE ~ "after"
    )
  )

df$period <- factor(df$period, levels = c("before", "during", "after"))
# group data 

df$count <- as.integer(df$count)

df$biomass_g <- as.numeric(df$biomass_g)

df$Year <- substr(df$Date, 1, 4)
df$date_index <- paste(df$Year, df$Month, sep = "_")
df$date_index <- gsub(" ", "", df$date_index)

df$ym <-ym(paste(df$Year, df$Month, sep = "-"))

df$ym <- substr(df$ym, 1,7)

# for comparison AGRRA protocol does not include all fish so exclude fish not included in the 2018 survey


# for comparison AGRRA protocol does not include all fish so exclude fish not included in the 2018 survey

# vector of agrra species
agrra_sp <- c("Pomacanthus_paru", "Pomacanthus_arcuatus", "Holacanthus_ciliaris", "Holacanthus_tricolor", "Chaetodon_striatus", "Chaetodon_capistratus", "Prognathodes_aculeatus", "Chaetodon_sedentarius", "Chaetodon_ocellatus", "Anisotremus_surinamensis", "Haemulon_sciurus", "Haemulon_carbonarium", "Haemulon_melanurum", "Haemulon_flavolineatum", "Anisotremus_virginicus", "Haemulon_parra", "Haemulon_chrysargyreum", "Haemulon_macrostomum", "Haemulon_aurolineatum", "Haemulon_plumierii", "Haemulon_album", "Haemulon_/Anisotremus", "Scarus_coeruleus", "Sparisoma_atomarium", "Scarus_coelestinus", "Scarus_taeniopterus", "Scarus_vetula", "Scarus_guacamaia", "Sparisoma_aurofrenatum", "Sparisoma_chrysopterum", "Sparisoma_viride", "Scarus_iseri", "Sparisoma_rubripinne", "Scarus/Sparisoma", "Mycteroperca_bonaci", "Cephalopholis_fulva", "Cephalopholis_cruentata", "Epinephelus_striatus", "Epinephelus_guttatus", "Epinephelus_adscensionis", "Mycteroperca_tigris", "Mycteroperca_venenosa", "Mycteroperca_interstitialis", "Lutjanus_cyanopterus", "Lutjanus_jocu", "Lutjanus_griseus", "Lutjanus_synagris", "Lutjanus_mahogoni", "Lutjanus_analis", "Lutjanus_apodus", "Ocyurus_chrysurus", "Acanthurus_coeruleus", "Acanthurus_chirurgus", "Acanthurus_tractus", "Melichthys_niger", "Canthidermis_sufflamen", "Balistes_vetula", "Xanthichthys_ringens", "Lachnolaimus_maximus", "Halichoeres_radiatus", "Halichoeres_bivittatus", "Bodianus_rufus", "Halichoeres_garnoti", "Cantherhines_pullus", "Aluterus_scriptus", "Cantherhines_macrocerus", "Calamus_bajonado", "Calamus_pennatula", "Calamus_calamus", "Calamus_penna", "Diodon_holocanthus", "Diodon_hystrix", "Sphoeroides_spengleri", "Caranx_ruber", "Kyphosus_spp.", "Sphyraena_barracuda", "Pterois_volitans", "Trachinotus_falcatus", "Lactophrys_bicaudalis_", "Stegastes_planifrons", "Microspathodon_chrysurus")

unique(df$latin_name)

new_df <- subset(df, latin_name %in% agrra_sp)

df <- new_df

overall_df <- df %>%
  group_by(ym, site, Month, Year, period) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))

tapply(overall_df$mean_biomass, overall_df$period, median)

comp_all <- ggplot(overall_df, aes(period, mean_biomass/1000)) +
  geom_boxplot() +
  xlab("") + ylab("Mean biomass (kg)") + ggtitle("All fish") +
  my_theme

comp_all <- comp_all + geom_signif(comparisons = list(c("before", "during"), c("during", "after"), c("before", "after")), 
                                   test = "wilcox.test", map_signif_level = T,
                                   textsize = 5)

comp_all

# make figures even sexier 

fishualize(option="Sparisoma_viride")
fish(option="Sparisoma_viride", n=6)

trophic_df <- df %>%
  group_by(ym, site, Month, Year, trophic_guild, period) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))

tapply(trophic_df$mean_biomass, trophic_df$period, median)

comp_herbs <- ggplot(trophic_df, aes(period, mean_biomass/1000, fill = period)) +
  geom_boxplot(fill="#9BE26DFF") +
  xlab("") + ylab("") + ggtitle("Herbivores") +
  my_theme
comp_herbs


comp_herbs <- comp_herbs + geom_signif(comparisons = list(c("before", "during"), c("during", "after"), c("before", "after")), 
                                       test = "wilcox.test", map_signif_level = T,
                                       textsize = 5)

comp_herbs
# all parrotfish

subset_df <- df[grepl("terminal$|initial$", df$common_name), ]

summarised_df <-  subset_df %>%
  group_by(site, Month, Year, date_index, period) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))

tapply(summarised_df$mean_biomass, summarised_df$period, median)


all_parrot <- ggplot(summarised_df, aes(period, mean_biomass/1000, fill = period)) +
  geom_boxplot(fill = "#56D6FFFF") + 
  xlab("Covid time period") + ylab("Mean biomass (kg)") + ggtitle("All Parrotfish") +
  my_theme

all_parrot <- all_parrot + geom_signif(comparisons = list(c("before", "during"), c("during", "after"), c("before", "after")), 
                                       test = "wilcox.test", map_signif_level = T,
                                       textsize = 5)
all_parrot

# filter the fish names that end in "initial"
initial_df <- df[grepl("initial$", df$common_name), ]

summarised_initial_df <- initial_df %>%
  group_by(site, Month, Year, date_index, period) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))

tapply(summarised_initial_df$mean_biomass, summarised_initial_df$period, median)

initial_parrot <- ggplot(summarised_initial_df, aes(period, mean_biomass/1000, fill = period)) +
  geom_boxplot(fill ="#FE7562FF") + 
  xlab("Covid time period") + ylab("") + ggtitle("Initital Phase Parrotfish") +
  my_theme

initial_parrot <- initial_parrot + geom_signif(comparisons = list(c("before", "during"), c("during", "after"), c("before", "after")), 
                                               test = "wilcox.test", map_signif_level = T,
                                               textsize = 5)

initial_parrot


png(file=file.path(output_wd, "before_during_after_FIGURE.png"), height = 4000, width = 3500, res = 350)
(comp_all + comp_herbs)/(all_parrot + initial_parrot) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, face = "bold"))
dev.off()

# compare biomass for all fish 
kruskal.test(overall_df$mean_biomass~overall_df$period)
dunnTest(overall_df$mean_biomass~overall_df$period)

# herbivores
herb_df <- subset(trophic_df, trophic_guild == "Herbivore")
kruskal.test(herb_df$mean_biomass~herb_df$period)
dunnTest(herb_df$mean_biomass~herb_df$period)

# all parrotfish 
kruskal.test(summarised_df$mean_biomass~summarised_df$period)
dunnTest(summarised_df$mean_biomass~summarised_df$period)

# intial parrotfish 
kruskal.test(summarised_initial_df$mean_biomass~summarised_initial_df$period)
dunnTest(summarised_initial_df$mean_biomass~summarised_initial_df$period)



