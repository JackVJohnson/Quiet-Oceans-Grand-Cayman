#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: mvabund analysis and heatmap 
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: November 2023
# Last edited: November 2023

# clear workspace
rm(list=ls())


# libraries

library(tidyverse)
library(vegan)
library(patchwork)
library(mvabund)




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

# make 2018 vector of fish
agrra_sp <- c("Pomacanthus_paru", "Pomacanthus_arcuatus", "Holacanthus_ciliaris", "Holacanthus_tricolor", "Chaetodon_striatus", "Chaetodon_capistratus", "Prognathodes_aculeatus", "Chaetodon_sedentarius", "Chaetodon_ocellatus", "Anisotremus_surinamensis", "Haemulon_sciurus", "Haemulon_carbonarium", "Haemulon_melanurum", "Haemulon_flavolineatum", "Anisotremus_virginicus", "Haemulon_parra", "Haemulon_chrysargyreum", "Haemulon_macrostomum", "Haemulon_aurolineatum", "Haemulon_plumierii", "Haemulon_album", "Haemulon_/Anisotremus", "Scarus_coeruleus", "Sparisoma_atomarium", "Scarus_coelestinus", "Scarus_taeniopterus", "Scarus_vetula", "Scarus_guacamaia", "Sparisoma_aurofrenatum", "Sparisoma_chrysopterum", "Sparisoma_viride", "Scarus_iseri", "Sparisoma_rubripinne", "Scarus/Sparisoma", "Mycteroperca_bonaci", "Cephalopholis_fulva", "Cephalopholis_cruentata", "Epinephelus_striatus", "Epinephelus_guttatus", "Epinephelus_adscensionis", "Mycteroperca_tigris", "Mycteroperca_venenosa", "Mycteroperca_interstitialis", "Lutjanus_cyanopterus", "Lutjanus_jocu", "Lutjanus_griseus", "Lutjanus_synagris", "Lutjanus_mahogoni", "Lutjanus_analis", "Lutjanus_apodus", "Ocyurus_chrysurus", "Acanthurus_coeruleus", "Acanthurus_chirurgus", "Acanthurus_tractus", "Melichthys_niger", "Canthidermis_sufflamen", "Balistes_vetula", "Xanthichthys_ringens", "Lachnolaimus_maximus", "Halichoeres_radiatus", "Halichoeres_bivittatus", "Bodianus_rufus", "Halichoeres_garnoti", "Cantherhines_pullus", "Aluterus_scriptus", "Cantherhines_macrocerus", "Calamus_bajonado", "Calamus_pennatula", "Calamus_calamus", "Calamus_penna", "Diodon_holocanthus", "Diodon_hystrix", "Sphoeroides_spengleri", "Caranx_ruber", "Kyphosus_spp.", "Sphyraena_barracuda", "Pterois_volitans", "Trachinotus_falcatus", "Lactophrys_bicaudalis_", "Stegastes_planifrons", "Microspathodon_chrysurus")

unique(df$latin_name)

new_df <- subset(df, latin_name %in% agrra_sp)

df <- new_df

# 2. Data Transformation: Create a count matrix for fish species
df2 <- df %>%
  group_by(period, site, transect, date_index, latin_name) %>%
  summarise(total_abundance = sum(count))

df3 <- df2 %>%
  pivot_wider(names_from = latin_name, values_from = total_abundance)
df3[is.na(df3)] <- 0
df3 <- subset(df3, !period == "during")

# make it suitable for mvabund

df4 <- df3[,-c(1:4)]

abund_df <- mvabund(df4)

# for the sake of this analysis only (I dunno if I can do a pairwise on this?) just compare before with after
period <- df3$period
abund_glm <- manyglm(abund_df ~ period, family="negative.binomial")
summary(abund_glm)

a <- anova.manyglm(abund_glm, resamp = "perm.resid",  p.uni= "adjusted", nBoot = 50, test = "wald")
a_pvals <- a$uni.p
#a_pvals2 <- a_pvals[2,]
a_testvals <- a$uni.test
pj.sort = sort(a_pvals, index.return=TRUE)


png(file=file.path(output_wd, "Compare abndances before and after.png"), height = 2500, width = 4000, res = 350)
plot(abund_df ~ period, var.subset=a_pvals<0.1, cex = 1)
dev.off()

plot.manyglm(abund_glm)

a_both <- rbind(a_pvals, a_testvals)
a_both <- a_both[-c(1,3),]
newnames <- c("P value", "Test value")
rownames(a_both) <- newnames

write.csv(a_both, file=file.path(output_wd, "manyglmoutput_eachspecies.csv"))


###########################################################################################
##################################### heat map ############################################

# descending df of biomass

biomass_df <- df %>%
  group_by(latin_name, period) %>%
  summarise(mean_biomass =  mean(biomass_g)) %>%
  arrange(desc(mean_biomass)) %>%
  head(50)


# Assuming your dataframe is named 'top_species_df'
# Convert 'outplant_time' to a factor for correct ordering on the heatmap
biomass_df$period <- factor(biomass_df$period, levels = c("before", "during", "after"))

hist(log1p(biomass_df$mean_biomass))
summary(log1p(biomass_df$mean_biomass))

biomass_df$latin_name <- gsub("_", " ", biomass_df$latin_name)

# Create the heatmap
heatmap_plot <- ggplot(biomass_df, aes(x = reorder(latin_name, mean_biomass), y = period, fill = log1p(mean_biomass))) +
  geom_tile(color = "white") +
  labs(x = "Fish Species", y = "Period", fill = "", title = "Mean biomass (log)") +
  scale_fill_viridis_c(option = "A", direction = -1, breaks = c(6.8, 8.889), labels = c("Low", "High")) +
  #scale_fill_gradient(low="blue", high = "red")+
  #scale_fill_gradient(colours=viridis(4))+
  coord_flip() +
  my_theme +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5, vjust = 1.5, face = "plain", size = 14), axis.text.y = element_text(face = "italic"))
#theme(legend.title = element_text(hjust = 0.5, vjust = 1.8))
heatmap_plot

length(unique(biomass_df$latin_name))

png(file=file.path(output_wd, "species_heatmap_QO.png"), height = 4000, width = 3000, res = 350)
heatmap_plot
dev.off()
