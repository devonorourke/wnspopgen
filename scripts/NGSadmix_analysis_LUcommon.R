library(tidyverse)
library(ggpubr)

## written: 18-feb-2020 by Devon O'Rourke
## see: http://www.popgen.dk/software/index.php/NgsAdmixTutorial

# theme function for custom plot style:
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## import data from ngsAdmix
## 1. file paths for import
k2_angsd <- "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/NGSadmix/LUcommon_k2.qopt"
k3_angsd <- "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/NGSadmix/LUcommon_k3.qopt"
k4_angsd <- "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/NGSadmix/LUcommon_k4.qopt"

## 2. import sample labels to data:
meta <- read_delim(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PLINK/mylu_proper.fam",
                   delim = '\t', col_names = FALSE) %>% 
  select(`X2`) %>% rename(sample = `X2`)

## 3. import each of the three K-n sizes (K=2|3|4) from each of the SNP sets
  ## importing pairs of files because of uniqe data structure between K2-4...
  ## ...(LUcommon from all Angsd, LUregions from VCF set)

k2_importFunction <- function(filepath1) {
  angsd_data <- read_delim(file = filepath1, delim = " ", col_names = FALSE) %>% 
    select(-`X3`) %>% 
    rename(pop1 = `X1`, pop2 = `X2`) %>% 
    mutate(AncPops="2K", SNPset="angsd") %>% 
    cbind(., meta) %>% 
    pivot_longer(., cols = c('pop1', 'pop2'), names_to = "AncPop")
}

k3_importFunction <- function(filepath1) {
  angsd_data <- read_delim(file = filepath1, delim = " ", col_names = FALSE) %>% 
    select(-`X4`) %>% 
    rename(pop1 = `X1`, pop2 = `X2`, pop3 = `X3`) %>% 
    mutate(AncPops="3K", SNPset="angsd") %>% 
    cbind(., meta) %>% 
    pivot_longer(., cols = c('pop1', 'pop2', 'pop3'), names_to = "AncPop")
}

k4_importFunction <- function(filepath1) {
  angsd_data <- read_delim(file = filepath1, delim = " ", col_names = FALSE) %>% 
    select(-`X5`) %>% 
    rename(pop1 = `X1`, pop2 = `X2`, pop3 = `X3`, pop4 = `X4`) %>% 
    mutate(AncPops="4K", SNPset="angsd") %>% 
    cbind(., meta) %>% 
    pivot_longer(., cols = c('pop1', 'pop2', 'pop3', 'pop4'), names_to = "AncPop")
}

k2_data <- k2_importFunction(k2_angsd)
k3_data <- k3_importFunction(k3_angsd)
k4_data <- k4_importFunction(k4_angsd)

all_data <- rbind(k2_data, k3_data, k4_data)

## add location metadata to data:
mylu_metadata <- read_csv("https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/metadata/resolved_sample_metadata.csv") %>% 
  filter(libPrepAlias %in% unique(all_data$sample))
mylu_metadata$Location <- str_replace(mylu_metadata$Location, "-", "\n")

all_data <- merge(all_data, 
                  mylu_metadata,
                  by.x="sample", by.y="libPrepAlias")

## making the structure plots
## creating individual plots for each SNPset and K size; then stitching together
all_data$Location <- factor(all_data$Location,
                              levels = c("NY\nWilliams","NY\nHailesCave","NY\nunknown",
                                         "VT\nAeolus","VT\nNewfane","VT\nStockbridge",
                                         "NH\nCharlestown", "NH\nMilford", 
                                         "MA\nChester", "MA\nPrinceton", "MA\nPepperell", "MA\nLincoln"))

p_a2 <- all_data %>% filter(AncPops=="2K" & SNPset=="angsd") %>% 
  arrange(desc(Location)) %>%
  ggplot(aes(x=reorder(sample, -LocationOrder), y=value, fill=AncPop, label=LocationNumber)) +
  geom_bar(stat="identity", 
           width = 0.92) +
  scale_fill_manual(values = c("dodgerblue1", "red3", "orange3", "springgreen4")) +
  theme_devon() +
  theme(axis.text.x = element_blank(),
        #axis.text.x = element_text(size=6, angle=45, hjust=1),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 6)) +
  labs(x="", y="Admixture", fill = "", subtitle = "Angsd SNP set, 2 Ancestral populations") +
  facet_grid (. ~ Location, scales="free", space="free")


p_a3 <- all_data %>% filter(AncPops=="3K" & SNPset=="angsd") %>% 
  arrange(desc(Location)) %>%
  ggplot(aes(x=reorder(sample, -LocationOrder), y=value, fill=AncPop, label=LocationNumber)) +
  geom_bar(stat="identity", 
           width = 0.92) +
  scale_fill_manual(values = c("dodgerblue1", "red3", "orange3", "springgreen4")) +
  theme_devon() +
  theme(axis.text.x = element_blank(),
        #axis.text.x = element_text(size=6, angle=45, hjust=1),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 6)) +
  labs(x="", y="Admixture", fill = "", subtitle = "Angsd SNP set, 3 Ancestral populations") +
  facet_grid (. ~ Location, scales="free", space="free")


p_a4 <- all_data %>% filter(AncPops=="4K" & SNPset=="angsd") %>% 
  arrange(desc(Location)) %>%
  ggplot(aes(x=reorder(sample, -LocationOrder), y=value, fill=AncPop, label=LocationNumber)) +
  geom_bar(stat="identity", 
           width = 0.92) +
  scale_fill_manual(values = c("dodgerblue1", "red3", "orange3", "springgreen4")) +
  theme_devon() +
  theme(axis.text.x = element_blank(),
        #axis.text.x = element_text(size=6, angle=45, hjust=1),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 6)) +
  labs(x="", y="Admixture", fill = "", subtitle = "Angsd SNP set, 4 Ancestral populations") +
  facet_grid (. ~ Location, scales="free", space="free")

## stitch together
ggarrange(p_a2, p_a3, p_a4, nrow=3, ncol=1)

## are any samples dominated by K2 or K3 pops? (most of plots are K1)
data_k2 <- all_data %>% 
  filter(AncPops == "2K") %>% 
  select(AncPop, value, sample) %>% 
  pivot_wider(names_from = AncPop,
              values_from = value)

data_k3 <- all_data %>% 
  filter(AncPops == "3K") %>% 
  select(AncPop, value, sample) %>% 
  pivot_wider(names_from = AncPop,
              values_from = value)

data_k4 <- all_data %>% 
  filter(AncPops == "4K") %>% 
  select(AncPop, value, sample) %>% 
  pivot_wider(names_from = AncPop,
              values_from = value)
