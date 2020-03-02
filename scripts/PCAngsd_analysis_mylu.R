library(RcppCNPy)
library(tidyverse)
library(ggrepel)
#library(ggpubr)

## written: 19-feb-2020 by Devon O'Rourke
## see: https://github.com/Rosemeis/pcangsd and http://www.popgen.dk/software/index.php/PCAngsd
## see original paper too: https://www.genetics.org/content/210/2/719
## scripts: see 'pcangsd_LUcommon.sh' and 'pcangsd_LUregions.sh'
## r script to produce:
## 1. an(other) admixture plot with K=2 only estimates
## 2. inbreeding coefficients estimated using parameters in ngsF (model 'inbreed 1' in pcangsd program)
## 3. PCA plot for outlier analysis

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

## import sample labels to data:
meta <- read_delim(file = "https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PLINK/mylu_proper.fam",
                   delim = '\t', col_names = FALSE) %>% 
  select(`X2`) %>% rename(sample = `X2`)

## download .npy files from GitHub repo:
  ## get 'LUcommon.admix.Q.npy' from: https://github.com/devonorourke/wnspopgen/raw/master/data/PCAngsd/LUcommon.admix.Q.npy
  ## get 'LUcommon.inbreed.npy' from: https://github.com/devonorourke/wnspopgen/raw/master/data/PCAngsd/LUcommon.inbreed.npy

# Read in estimated admixutre matrix
  ## Note that the current working directory is set to my desktop
  ## The 'npyLoad' function didn't like the full path specified for some reason within the "filename=" argument

## change this for your own script!!
setwd("~/Documents/nau_projects/bat_genomes/popgen/pcAngsd/")

## load file
lu_common_admx <- as.data.frame(npyLoad("LUcommon.admix.Q.npy")) %>% 
  rename(K1 = `V1`, K2 = `V2`) %>% 
  cbind(., meta) %>% 
  pivot_longer(., cols = c('K1', 'K2'), names_to = "AncPop")

## read in metadata; merge with admix data
mylu_metadata <- read_csv("https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/metadata/resolved_sample_metadata.csv") %>% 
  filter(libPrepAlias %in% unique(lu_common_admx$sample))
mylu_metadata$Location <- str_replace(mylu_metadata$Location, "-", "\n")
lu_common_admx <- merge(lu_common_admx, mylu_metadata, by.x="sample", by.y="libPrepAlias")

## admixture-style plot; save as 'PCAngsd_admixture_LU_allSamps'
lu_common_admx$Location <- factor(lu_common_admx$Location,
                            levels = c("NY\nWilliams","NY\nHailesCave","NY\nunknown",
                                       "VT\nAeolus","VT\nNewfane","VT\nStockbridge",
                                       "NH\nCharlestown", "NH\nMilford", 
                                       "MA\nChester", "MA\nPrinceton", "MA\nPepperell", "MA\nLincoln"))

lu_common_admx %>% 
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
  labs(x="", y="Admixture", fill = "", subtitle = "PCAngsd admixture estimates using Angsd SNP set, 2 Ancestral populations inferred") +
  facet_grid (. ~ Location, scales="free", space="free")

# Read in inbreeding matrix
lu_common_inbr <- as.data.frame(npyLoad("LUcommon.inbreed.npy"))
colnames(lu_common_inbr) <- "inbrEst"
lu_common_inbr <- cbind(lu_common_inbr, meta)
## tcap188 has large positive estimate, indicative of exess homozygosity
## a few negative samples too - indicating exess of heterozygosity

## generate PC plots from covariance matrix:
lu_common_cov <- as.matrix(read_delim("https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PCAngsd/LUcommon.cov", delim = " ", col_names = FALSE))
lu_common_ev <- eigen(lu_common_cov)
lu_common_ev_vectors <- lu_common_ev$vectors
lu_common_ev_vectors <- cbind(lu_common_ev_vectors, meta)

## scale x/y in similar proportions:
val1 <- max(lu_common_ev_vectors$`1`) - min(lu_common_ev_vectors$`1`)
val2 <- max(lu_common_ev_vectors$`2`) - min(lu_common_ev_vectors$`2`)
#span values to 1.2 total...
## set axis2 at 0.6 to -0.6
## reset axis1 to -1.0 to 0.2

ggplot(lu_common_ev_vectors, label=sample) +
  geom_point(aes(x=`1`, y=`2`)) +
  scale_x_continuous(limits=c(-1, 0.2), breaks=seq(-1, .2, 0.4)) +
  scale_y_continuous(limits=c(-0.6, 0.6), breaks=seq(-0.6, 0.6, 0.4)) +
  theme_devon() +
  labs(x="PC1", y="PC2", subtitle = "PCAngsd PCA plot generated from estimated covariance matrix of all mylu samples") +
  geom_text_repel(data = lu_common_ev_vectors %>% filter(.$`1` < -0.2),
                  aes(x=`1`, y=`2`, label=sample), color="red")



