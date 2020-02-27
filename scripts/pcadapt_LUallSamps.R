##notrun:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

library(pcadapt)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(qvalue)

## get names of samples
pop_info <- read_delim("https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/PLINK/mylu_proper.fam", 
                       col_names = FALSE, delim = "\t") %>% 
  rename(Pop=`X1`, Indiv=`X2`) %>% 
  select(Pop, Indiv) %>% 
  mutate(PopNames=paste(Pop, Indiv, sep="-"))

## download VCF file from OSF site: https://osf.io/z2e74/download
vcfPath <- "LUall_pruned.vcf.gz"
fullvcf <- read.pcadapt(vcfPath, type = "vcf")

## Screeplot to check for K values to input
allsites.list.k10 <- pcadapt(input = fullvcf, K = 10, min.maf = 0.05)
plot(allsites.list.k10, option = "screeplot")
  ## looks like just 3 PC's are sufficient for most variance to be captured

## Import data with K=3
allsites.list.k3 <- pcadapt(input = fullvcf, K = 3, min.maf = 0.05)

## collect data from this file:
scores_dat <- data.frame(allsites.list.k3$scores) %>% 
  rename(PC1=`X1`, PC2=`X2`, PC3=`X3`)
scores_dat <- cbind(scores_dat, pop_info)

maxval_pc <- max(scores_dat$PC1, scores_dat$PC2, scores_dat$PC3) + 0.05
minval_pc <- min(scores_dat$PC1, scores_dat$PC2, scores_dat$PC3) - 0.05

## plot PC1-PC2
pc12 <- ggplot(scores_dat,
       aes(x=PC1, y=PC2, color=Pop, label=Indiv)) +
  geom_point() +
  geom_text_repel(data = scores_dat %>% filter(PC2 < -0.2 | PC1 > 0.5)) +
  scale_x_continuous(limits = c(minval_pc, maxval_pc)) +
  scale_y_continuous(limits = c(minval_pc, maxval_pc))

## plot PC1-PC3
pc13 <- ggplot(scores_dat,
       aes(x=PC1, y=PC3, color=Pop, label=Indiv)) +
  geom_point() +
  geom_text_repel(data = scores_dat %>% filter(PC3 > 0.5 | PC3 < -0.08)) +
  geom_text_repel(data = scores_dat %>% filter(PC1 > 0.5)) + 
  scale_x_continuous(limits = c(minval_pc, maxval_pc)) +
  scale_y_continuous(limits = c(minval_pc, maxval_pc))

## plot PC2-PC3
pc23 <- ggplot(scores_dat,
       aes(x=PC2, y=PC3, color=Pop, label=Indiv)) +
  geom_point() +
  geom_text_repel(data = scores_dat %>% filter(PC3 > 0.5 | PC3 < -0.08)) +
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  scale_y_continuous(limits = c(-0.8, 0.8))

## stitch together:
ggarrange(pc12, pc13, pc23, 
          common.legend = TRUE, ncol = 3, nrow=1)