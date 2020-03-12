library(tidyverse)

## written: 30-jan-2020 by Devon O'Rourke
## see: https://github.com/ANGSD/NgsRelate
## scripts: see 'angsd_allLU_ngsRelate.sh'

## import data from NGSrelate
data <- read_delim(file = "https://github.com/devonorourke/wnspopgen/raw/master/data/ngsRelate/LUregionsNGSrelate.res.gz",
                   delim = "\t", col_names = TRUE)
meta <- read_delim(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PLINK/mylu_proper.fam",
                   delim = '\t', col_names = FALSE) %>% 
  select(`X2`) %>% rename(sample = `X2`)

meta <- meta %>% 
  mutate(numID=as.numeric(row.names(meta))) %>% 
  mutate(numID=numID-1)

data <- merge(data, meta, by.x = 'a', by.y = 'numID')
data <- data %>% rename(sample_A = sample)

data <- merge(data, meta, by.x = 'b', by.y = 'numID')
data <- data %>% rename(sample_B = sample)

rm(meta)

## data summaries: any significant outliers?
data %>% 
  summarise(meanSites = mean(nSites),
            maxSites = max(nSites),
            minSites = min(nSites),
            sdSites = sd(nSites))


## visualizing various stats:
## order the x/y axis according to colnames 'a' and 'b'
## save as "ngsRelate_allMylu"; export at 500x500
ggplot(data,
       aes(x=reorder(sample_A, a), y=reorder(sample_B,b), fill=rab)) +
  geom_tile() +
  #theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text.x = element_text(size = 4, angle=90, hjust=1),
        #axis.text.y = element_text(size = 4),
        legend.position = "right") +
  scale_fill_viridis_c() +
  labs(x="", y="", fill="relatedness", subtitle = "pairwise relatedness following Hedrick et al. (2015) - all MYLU")

## appears to be a handful of pairs of samples with relatedness above 0.5
filt_data <- data %>% 
  filter(rab > 0.1) %>% 
  select(rab, sample_A, sample_B)
  ## as with other estimates of relatedness, tcap006 and tcap027 are the two major outliers among all other data

