##notrun:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

library(pcadapt)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(qvalue)

## install theme for plot
theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


## get names of samples
pop_info <- read_delim("https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PLINK/mylu_proper.fam", 
                       col_names = FALSE, delim = "\t") %>% 
  rename(Pop=`X1`, Indiv=`X2`) %>% 
  select(Pop, Indiv) %>% 
  mutate(PopNames=paste(Pop, Indiv, sep="-"))

## download VCF file from OSF site - see: https://osf.io/97mc6/
vcfPath <- "LUcommon_maf05_noRepeats_plinkLD_filtd.vcf.gz"
fullvcf <- read.pcadapt(vcfPath, type = "vcf")

## Screeplot to check for K values to input
allsites.list.k10 <- pcadapt(input = fullvcf, K = 10, min.maf = 0.05)
plot(allsites.list.k10, option = "screeplot") + theme_devon()
## save plot as 'screeplot_k10_LU'; export at w600 x h610
  ## looks like just 2 PC's are sufficient for most variance to be captured; will gather 3 for test
  ## note the drop in proportion of varance is steep, but the absolute difference is REALLY small.

## Import data with K=3
allsites.list.k3 <- pcadapt(input = fullvcf, K = 3, min.maf = 0.05)
rm(allsites.list.k10)

## collect data from this file:
scores_dat <- data.frame(allsites.list.k3$scores) %>% 
  rename(PC1=`X1`, PC2=`X2`, PC3=`X3`)
scores_dat <- cbind(scores_dat, pop_info)

maxval_pc <- max(scores_dat$PC1, scores_dat$PC2, scores_dat$PC3) + 0.05
minval_pc <- min(scores_dat$PC1, scores_dat$PC2, scores_dat$PC3) - 0.05

## plot PC1-PC2
pc12 <- ggplot(scores_dat,
       aes(x=PC1, y=PC2, color=Pop, label=Indiv, shape=Pop)) +
  geom_point(alpha=0.5) +
  geom_text_repel(data = scores_dat %>% filter(PC2 < -0.2 | PC1 > 0.5)) +
  scale_x_continuous(limits = c(minval_pc, maxval_pc)) +
  scale_y_continuous(limits = c(minval_pc, maxval_pc)) +
  scale_color_manual(values = c("brown", "blue")) +
  theme_devon()
  
## plot PC1-PC3
pc13 <- ggplot(scores_dat,
       aes(x=PC1, y=PC3, color=Pop, label=Indiv, shape=Pop)) +
  geom_point(alpha=0.5) +
  geom_text_repel(data = scores_dat %>% filter(PC3 > 0.4)) +
  geom_text_repel(data = scores_dat %>% filter(PC1 > 0.5)) + 
  scale_x_continuous(limits = c(minval_pc, maxval_pc)) +
  scale_y_continuous(limits = c(minval_pc, maxval_pc)) +
  scale_color_manual(values = c("brown", "blue")) +
  theme_devon()

## plot PC2-PC3
pc23 <- ggplot(scores_dat,
       aes(x=PC2, y=PC3, color=Pop, label=Indiv, shape=Pop)) +
  geom_point(alpha=0.5) +
  geom_text_repel(data = scores_dat %>% filter(PC2 < -0.5)) +
  geom_text_repel(data = scores_dat %>% filter(PC3 > 0.5)) +
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  scale_y_continuous(limits = c(-0.8, 0.8)) +
  scale_color_manual(values = c("brown", "blue")) +
  theme_devon()

## stitch together:
ggarrange(pc12, NULL, pc13, pc23, 
          common.legend = TRUE, ncol = 2, nrow=2)

## save as 'pcadapt_PCA_LU_allSamps'; export at w600 x h600
