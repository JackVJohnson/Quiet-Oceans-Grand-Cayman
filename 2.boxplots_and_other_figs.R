#########################################################################################
#################################### Quiet Oceans #######################################
#########################################################################################

# objective: boxplots and other figures 
# Author: Jack Johnson (jackvjohnson@hotmail.com)
# Date created: April 2023
# Last edited: 

# clear workspace
rm(list=ls())

# libraries

library(tidyverse)
library(fishualize)
library(patchwork)
library(ggcorrplot)
library(ggrepel)
library(viridis)


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

# group data 

df$count <- as.integer(df$count)

df$biomass_g <- as.numeric(df$biomass_g)

df$Year <- substr(df$Date, 1, 4)
df$date_index <- paste(df$Year, df$Month, sep = "_")
df$date_index <- gsub(" ", "", df$date_index)

df$ym <-ym(paste(df$Year, df$Month, sep = "-"))

df$ym <- substr(df$ym, 1,7)

overall_df <- df %>%
  group_by(ym, site, Month, Year) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))




# boxplots

p1 <- ggplot(overall_df, aes(as.factor(ym), mean_richness)) +
  geom_boxplot() +
  geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="Mean species richness")

p2 <- ggplot(overall_df, aes(as.factor(ym), log1p(mean_biomass))) +
  geom_boxplot() +
  geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="Mean biomass (log)")

p3 <- ggplot(overall_df, aes(as.factor(ym), log1p(mean_abundance))) +
  geom_boxplot() +
  geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = .5, vjust=.5)) +
  labs(x="Time period", y="Mean abundance (log)")

png(file=file.path(output_wd, "allfish_boxplots_FIGURE.png"), height = 4000, width = 4000, res = 350)
p1/p2/p3 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 30, face = "bold"))
dev.off()

df2 <- df %>%
  group_by(ym, site, Month, Year, trophic_guild) %>%
  summarise(mean_biomass = sum(biomass_g)/(length(unique(transect))),
            mean_abundance = sum(count)/length(unique(transect)),
            mean_richness = n_distinct(latin_name)/length(unique(transect)))

p4 <- ggplot(df2, aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="mean Species richness") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), linetype = "dashed", color = "black")

p5 <- ggplot(df2, aes(as.factor(ym), log1p(mean_biomass), fill = trophic_guild)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="Mean biomass (log)") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), linetype = "dashed", color = "black")

p6 <- ggplot(df2, aes(as.factor(ym), log1p(mean_abundance), fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = .5, vjust =.5)) +
  labs(x="Time period", y="Mean abundance (log)") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), linetype = "dashed", color = "black")


png(file=file.path(output_wd, "trophic_boxplots_FIGURE.png"), height = 4000, width = 4000, res = 350)
p4/p5/p6 + plot_layout(guides="collect") & theme(legend.position = 'top')
dev.off()

# facet wrap is ugly and always has been ugly 

p7 <-  ggplot(df2, aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  my_theme +
  #theme(axis.text.x = element_blank()) +
  labs(x="Date", y="Mean species richness") +
  facet_wrap(~trophic_guild)
p7

# plots for each trophic guild individually 

p8 <- ggplot(subset(df2, trophic_guild == "Herbivore"), aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="Herbivores", title = "Mean species richness") 
p8 

p9 <- ggplot(subset(df2, trophic_guild == "Herbivore"),aes(as.factor(ym), log1p(mean_biomass), fill = trophic_guild)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="", title = "Mean biomass (log)") 
p9

p10 <- ggplot(subset(df2, trophic_guild == "Herbivore"), aes(as.factor(ym), log1p(mean_abundance), fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.3) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="", title = "Mean abundance (log)") 


p8 + p9 + p10

# invertivore

p11 <- ggplot(subset(df2, trophic_guild == "Invertivore"), aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.5) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="Invertivores") 

p12 <- ggplot(subset(df2, trophic_guild == "Invertivore"),aes(as.factor(ym), log1p(mean_biomass), fill = trophic_guild)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="A", begin = 0.5) +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="") 

p13 <- ggplot(subset(df2, trophic_guild == "Invertivore"), aes(as.factor(ym), log1p(mean_abundance), fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.5) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="") 


# macrocarnivore

p14 <- ggplot(subset(df2, trophic_guild == "Macrocarnivore"), aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.65) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="Macrocarnivores") 

p15 <- ggplot(subset(df2, trophic_guild == "Macrocarnivore"),aes(as.factor(ym), log1p(mean_biomass), fill = trophic_guild)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="A", begin = 0.65) +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="") 

p16 <- ggplot(subset(df2, trophic_guild == "Macrocarnivore"), aes(as.factor(ym), log1p(mean_abundance), fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.65) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="") 

# omnivore

p17 <- ggplot(subset(df2, trophic_guild == "Omnivore"), aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.8) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="Omnivores") 

p18 <- ggplot(subset(df2, trophic_guild == "Omnivore"),aes(as.factor(ym), log1p(mean_biomass), fill = trophic_guild)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="A", begin = 0.8) +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="") 

p19 <- ggplot(subset(df2, trophic_guild == "Omnivore"), aes(as.factor(ym), log1p(mean_abundance), fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.8) +
  my_theme +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(x="", y="") 

# planktivore

# Set the desired tick positions

p20 <- ggplot(subset(df2, trophic_guild == "Planktivore"), aes(as.factor(ym), mean_richness, fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.9) +
  my_theme +
  theme(legend.position = "none") +
  labs(x="Date", y="Planktivores") +
  scale_x_discrete(breaks = c("2020-09", "2021-08", "2022-04"))
p20

p21 <- ggplot(subset(df2, trophic_guild == "Planktivore"),aes(as.factor(ym), log1p(mean_biomass), fill = trophic_guild)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="A", begin = 0.9) +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  my_theme +
  theme(legend.position = "none") +
  scale_x_discrete(breaks = c("2020-09", "2021-08", "2022-04")) +
  labs(x="Date", y="") 

p22 <- ggplot(subset(df2, trophic_guild == "Planktivore"), aes(as.factor(ym), log1p(mean_abundance), fill = trophic_guild)) +
  geom_boxplot() +
  #geom_smooth(method="loess", color = "#952ea0", size = 1.5, aes(group=1)) +
  scale_fill_viridis_d(option="A", begin = 0.9) +
  my_theme +
  theme(legend.position = "none") +
  labs(x="Date", y="") +
  scale_x_discrete(breaks = c("2020-09", "2021-08", "2022-04"))


png(file=file.path(output_wd, "trophic_boxplots_individual_FIGURE.png"), height = 5100, width = 5100, res = 350)
(p8 + p9 + p10) / (p11 + p12 + p13) / (p14 + p15 + p16) / (p17 + p18 + p19)  / (p20 + p21 +p22) 
dev.off()

# summary stats overall

# biomass
tapply(overall_df$mean_biomass, overall_df$ym, mean)
tapply(overall_df$mean_biomass, overall_df$ym, sd)
# abundance
tapply(overall_df$mean_abundance, overall_df$ym, mean)
tapply(overall_df$mean_abundance, overall_df$ym, sd)


# summary stats for trophic guild 

tapply(df2$mean_biomass, list(df2$trophic_guild, df2$ym), mean)
tapply(df2$mean_biomass, list(df2$trophic_guild, df2$ym), sd)

