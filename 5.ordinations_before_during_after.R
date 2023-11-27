#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: Ordinations for before, during, after 
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: April 2023
# Last edited: September 2023

# clear workspace
rm(list=ls())

# for pairwise adonis
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# libraries

library(tidyverse)
library(vegan)
library(patchwork)
library(pairwiseAdonis)




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

mat1 <- as.matrix(df3[,-c(1:4)])
mat1 <- sqrt(mat1)
dist1 <- vegdist(mat1, method = "bray")

# permnova

set.seed(36)
comp <- adonis2(dist1~as.factor(df3$period), data = df3, permutations = 9999)
comp

# pairwise permanova 
paircomp <- pairwise.adonis(dist1, factors=as.factor(df3$period))
paircomp
# nmds 

set.seed(36)
m1 <- metaMDS(mat1, distance = "bray", k=2, autotransform = T, trymax = 50)
m1

stressplot(m1)

# extract scores

scores <- data.frame(m1$points)

scores$period <- df3$period

# plot

hull_before <- scores[scores$period == "before",][chull(scores[scores$period == "before", c("MDS1", "MDS2")]),]
hull_during <- scores[scores$period == "during",][chull(scores[scores$period == "during", c("MDS1", "MDS2")]),]
hull_after <- scores[scores$period == "after",][chull(scores[scores$period == "after", c("MDS1", "MDS2")]),]

hull_df <- rbind(hull_before, hull_during, hull_after)

ord_1 <- #ggplot()
  ggplot(data=scores,aes(x=MDS1, y=MDS2, colour=period, fill = period)) +
  stat_ellipse(type='t',size =1, linetype = 1)+
  #geom_polygon(data=hull_df, aes(MDS1, MDS2, fill = period, group = period), alpha = 0.4) +
  geom_point(position=position_jitter(.1), shape=19, size = 1)+
  #geom_point(data=scores,aes(x=MDS1, y=MDS2, colour=period), size=1, alpha =.6) + 
  my_theme +
  #scale_color_fish_d(option = "Halichoeres_bivittatus") +
  #scale_fill_fish_d(option = "Halichoeres_bivittatus") +
  #scale_colour_viridis_d(option = "A", end = .8) +
  #scale_fill_viridis_d(option = "A", end = .8)+
  scale_color_manual(values=c("#8c2981","#fe9f6d", "#000004")) +
  scale_fill_manual(values=c( "#8c2981","#fe9f6d", "#000004")) +
  labs(color = "Period") +
  guides(fill="none") +
  #theme(legend.position = "none") +
  xlab("MDS1") + ylab("MDS2") +
  annotate("text", x = 0.4, y = 1, label = "Stress = 0.238", color = "black", size = 4.5, fontface = "bold")
ord_1


##########################################################################################
############################## not used in MS ordinations ################################

# compare only during lockdown and after the cruise ships returned

rm(df2, df3, m1, mat1, scores, dist1, comp)

# 2. Data Transformation: Create a count matrix for fish species
df2 <- subset(df, !(period == "before"))

df2 <- df2 %>%
  group_by(period, site, transect, date_index, latin_name) %>%
  summarise(total_abundance = sum(count))

df3 <- df2 %>%
  pivot_wider(names_from = latin_name, values_from = total_abundance)
df3[is.na(df3)] <- 0

mat1 <- as.matrix(df3[,-c(1:4)])
mat1 <- sqrt(mat1)
dist1 <- vegdist(mat1, method = "bray")

# permnova

set.seed(36)
comp <- adonis2(dist1~as.factor(df3$period), data = df3, permutations = 9999)
comp

# nmds 

set.seed(36)
m1 <- metaMDS(mat1, distance = "bray", k=2, autotransform = T, trymax = 50)
m1 # stress = 0.254

stressplot(m1)

# extract scores

scores <- data.frame(m1$points)

scores$period <- df3$period

# plot

hull_during <-scores[scores$period == "during",][chull(scores[scores$period == "during", c("MDS1", "MDS2")]),]

hull_after <-scores[scores$period == "after",][chull(scores[scores$period == "after", c("MDS1", "MDS2")]),]

hull_df <- rbind(hull_during, hull_after)

ord_2 <- #ggplot() +
  ggplot(data=scores,aes(x=MDS1, y=MDS2, colour=period, fill = period)) +
  stat_ellipse(type='t',size =1, linetype = 1)+
  #geom_polygon(data=hull_df, aes(MDS1, MDS2, fill = period, group = period), alpha = 0.4) +
  geom_point(position=position_jitter(.1), shape=19, size = 1)+
  #geom_point(data=scores,aes(x=MDS1, y=MDS2, colour=period), size=1, alpha =.6) + 
  my_theme +
  #scale_color_fish_d(option = "Halichoeres_bivittatus") +
  #scale_fill_fish_d(option = "Halichoeres_bivittatus") +
  #scale_colour_viridis_d(option = "A", end = .8) +
  #scale_fill_viridis_d(option = "A", end = .8)+
  scale_color_manual(values=c("#fe9f6d", "#000004")) +
  scale_fill_manual(values=c("#fe9f6d", "#000004")) +
  #labs(color = "Period") +
  #guides(fill="none") +
  theme(legend.position = "none") +
  xlab("MDS1") + ylab("MDS2")
ord_2

png(file=file.path(output_wd, "Ordinations_FIGURE.png"), height = 2500, width = 5000, res = 350)
ord_1 + ord_2 + plot_layout(guides="collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 30))
dev.off()

png(file=file.path(output_wd, "Ordination_FIGURE.png"), height = 2000, width = 3000, res = 350)
ord_1
dev.off()