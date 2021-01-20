library(tidyverse)
library(ggpubr)

## Script used to generate admixture plots
## written 21 Jan 2020 by Devon O'Rourke

## data import
mylu_k2 <- read_table(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/ADMIXTURE/LUadmixture.2.Q", 
                      col_names = c('pop1', 'pop2'))

mylu_k3 <- read_table(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/ADMIXTURE/LUadmixture.3.Q", 
                      col_names = c('pop1', 'pop2', 'pop3'))

mylu_k4 <- read_table(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/ADMIXTURE/LUadmixture.4.Q",
                      col_names = c('pop1', 'pop2', 'pop3', 'pop4'))

mylu_names <- read_delim(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PLINK/mylu_proper.fam", 
                         col_names = FALSE,
                         delim = "\t") %>% 
  select(c(`X1`, `X2`)) %>% 
  rename(Pop=`X1`, Indiv=`X2`) %>% 
  mutate(Name=paste(Pop, Indiv, sep = "-"))

## create data frames for plots
data_k2 <- cbind(mylu_k2, mylu_names) %>% 
  pivot_longer(., cols = c(pop1, pop2),
               names_to = "Ancestry")

data_k3 <- cbind(mylu_k3, mylu_names) %>% 
  pivot_longer(., cols = c(pop1, pop2, pop3),
               names_to = "Ancestry")

data_k4 <- cbind(mylu_k4, mylu_names) %>%
  pivot_longer(., cols = c(pop1, pop2, pop3, pop4),
               names_to = "Ancestry")


## import and filter metadata for MYLU only
sample_metadata <- read_csv(file="https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/metadata/resolved_sample_metadata.csv")
mylu_metadat <- sample_metadata %>% 
  filter(libPrepAlias %in% mylu_names$Indiv)
mylu_metadat$Location <- str_replace(mylu_metadat$Location, "-", "\n")

plotdat_k2 <- merge(data_k2, mylu_metadat, by.x='Name', by.y='analysisAlias')
plotdat_k2$LocationNumber <- ifelse(plotdat_k2$Ancestry=="pop1", plotdat_k2$LocationNumber, NA)

plotdat_k3 <- merge(data_k3, mylu_metadat, by.x='Name', by.y='analysisAlias')
plotdat_k3$LocationNumber <- ifelse(plotdat_k3$Ancestry=="pop1", plotdat_k3$LocationNumber, NA)

plotdat_k4 <- merge(data_k4, mylu_metadat, by.x='Name', by.y='analysisAlias')
plotdat_k4$LocationNumber <- ifelse(plotdat_k4$Ancestry=="pop1", plotdat_k4$LocationNumber, NA)


################################################################################
## making the plots; option to stitch together
################################################################################

## plot K2
## reorder x-axis to group data by Locations in similar states
plotdat_k2$Location <- factor(plotdat_k2$Location,
                              levels = c("NY\nWilliamsMine","NY\nHailesCave",
                                         "VT\nStockbridge", "VT\nAeolus","VT\nNewfane",
                                         "NH\nCharlestown", "NH\nMilford", 
                                         "MA\nChester", "MA\nPrinceton", "MA\nPepperell", "MA\nLincoln"))

p2 <- plotdat_k2 %>% 
  arrange(desc(Location)) %>%
  ggplot(aes(x=reorder(Indiv, -LocationOrder), y=value, fill=Ancestry, label=LocationNumber)) +
  geom_bar(stat="identity", 
           width = 0.92) +
  scale_fill_manual(values = c("dodgerblue1", "sienna1")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 6)) +
  labs(x="", y="Admixture", fill = "", subtitle = "K = 2") +
  facet_grid (. ~ Location, scales="free", space="free")



## plot K3
## reorder x-axis to group data by Locations in similar states
plotdat_k3$Location <- factor(plotdat_k3$Location,
                              levels = c("NY\nWilliamsMine","NY\nHailesCave",
                                         "VT\nStockbridge", "VT\nAeolus","VT\nNewfane",
                                         "NH\nCharlestown", "NH\nMilford", 
                                         "MA\nChester", "MA\nPrinceton", "MA\nPepperell", "MA\nLincoln"))

p3 <- plotdat_k3 %>% 
  arrange(desc(Location)) %>%
  ggplot(aes(x=reorder(Indiv, -LocationOrder), y=value, fill=Ancestry, label=LocationNumber)) +
  geom_bar(stat="identity", 
           width = 0.92) +
  scale_fill_manual(values = c("yellow3", "orchid3", "palegreen4")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 6)) +
  labs(x="", y="Admixture", fill = "", subtitle = "K = 3") +
  facet_grid (. ~ Location, scales="free", space="free")



## plot K4
## reorder x-axis to group data by Locations in similar states
plotdat_k4$Location <- factor(plotdat_k4$Location,
                              levels = c("NY\nWilliamsMine","NY\nHailesCave",
                                         "VT\nStockbridge", "VT\nAeolus","VT\nNewfane",
                                         "NH\nCharlestown", "NH\nMilford", 
                                         "MA\nChester", "MA\nPrinceton", "MA\nPepperell", "MA\nLincoln"))
p4 <- plotdat_k4 %>% 
  arrange(desc(Location)) %>%
  ggplot(aes(x=reorder(Indiv, -LocationOrder), y=value, fill=Ancestry, label=LocationNumber)) +
  geom_bar(stat="identity", 
           width = 0.92) +
  scale_fill_manual(values = c("dodgerblue1", "red3", "orange3", "springgreen4")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 6)) +
  labs(x="", y="Admixture", fill = "", subtitle = "K = 4") +
  facet_grid (. ~ Location, scales="free", space="free")

## stitch the admixture plots together
## save plot with 1200x600 w/h dimension; save as 'LU_admixture_allSamples'
p234 <- ggarrange(p2, p3, p4, common.legend = FALSE, nrow = 3)
annotate_figure(p234,
                top = text_grob("MYLU admixture estimates for 2-4 populations"))
ggsave("~/github/wnspopgen/figures_tables/admixture_k2-4.png", width = 24, height = 16, units="cm", dpi=150)
ggsave("~/github/wnspopgen/figures_tables/admixture_k2-4.pdf", width = 24, height = 16, units="cm", dpi=300)
  ## can manually fix the few labels that get cut off (WilliamsMine, Newfane, Milford, Lincoln)