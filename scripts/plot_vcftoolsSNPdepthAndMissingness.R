## loading vcftools outputs for per-site missinginess and per-site depth of coverage
## plots show distribution of each data type for figures used in manuscript

library(tidyverse)
library(scales)
library(ggpubr)

## load missigness:
miss_df <- read_delim(file="https://github.com/devonorourke/wnspopgen/raw/master/data/vcftools/LUcommon_maf05_noRepeats_siteMissingness.txt.lmiss.gz",
                      delim="\t",
                      col_names=TRUE)
miss_df$F_WITH <- 1-(miss_df$F_MISS)  ## plotting how many samples retained SNPs at loci (seems easier to articulate)

## bin the missingness data into 10% buckets (i.e. 0.3-0.4, 0.41-0.5, 0.51-0.6...etc.)
min(miss_df$F_WITH) ## smallest bin is 0.3 (and max is 1)
nRecordsMiss <- length(miss_df$F_MISS)
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy} ## function to bin into coverages by 10x chunks (0-10X, 11-20X, etc.)
missBins <- data.frame(table(round_any(miss_df$F_WITH, 0.1, f = ceiling)))
missBins$cumFreqFromMaxMiss <- cumsum(missBins$Freq)
missBins$VarSort <- as.character(missBins$Var1)
missBins <- missBins %>% arrange(rev(VarSort)) %>% mutate(cumFreqFromMinMiss = cumsum(Freq))
missBins <- missBins %>% mutate(fracMaxMiss = cumFreqFromMaxMiss / nRecordsMiss,
                    fracMinMiss = cumFreqFromMinMiss / nRecordsMiss)
  ## so we can see that about 67% of SNPs are present in at least 70% of samples, and ...
    ##... 90% of samples have SNPs detected in at least 50% of the samples


## plot the frequency of samples per-site containing a SNP
p_miss <- ggplot(miss_df, aes(F_WITH)) + 
  geom_histogram(fill = "gray50", colour = "black", alpha=0.5, bins = 35) +
  scale_y_continuous(labels=comma) +
  labs(x="fraction of samples with data at loci", y="number of loci") +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

## save plot
p_miss
ggsave("~/github/wnspopgen/figures_tables/persite_missingness.png", height = 5, width = 6, dpi=150)
ggsave("~/github/wnspopgen/figures_tables/persite_missingness.pdf", height = 5, width = 6, dpi=300)

rm(miss_df)

## load depth data
snpdepth_df <- read_delim(file="https://github.com/devonorourke/wnspopgen/raw/master/data/vcftools/LUcommon_maf05_noRepeats_siteMeanDepth.txt.ldepth.mean.gz",
                      delim="\t",
                      col_names=TRUE)

## bin these data into groups of 10x, 20x, 30x, etc.
max(snpdepth_df$MEAN_DEPTH) ## largest bin is 300-310
nRecords <- length(snpdepth_df$MEAN_DEPTH)
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy} ## function to bin into coverages by 10x chunks (0-10X, 11-20X, etc.)
covBins <- data.frame(table(round_any(snpdepth_df$MEAN_DEPTH, 10, f = ceiling))) 
covBins$cumFreq <- cumsum(covBins$Freq)
covBins$fracFreq = (covBins$cumFreq)/nRecords
  ## so we find that between 

## plot depth data
p_depth <- ggplot(snpdepth_df, aes(MEAN_DEPTH)) + 
  geom_histogram(fill = "gray50", colour = "black", alpha=0.5, bins = 31) +
  scale_y_continuous(labels=comma) +
  xlim(10,40) +
  labs(x="mean depth", y="number of loci") +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

## save plot
p_depth
ggsave("~/github/wnspopgen/figures_tables/persite_depth.png", height = 5, width = 6, dpi=150)
ggsave("~/github/wnspopgen/figures_tables/persite_depth.pdf", height = 5, width = 6, dpi=300)

## plot both together and save:
ggarrange(p_depth, p_miss,
          ncol = 2, nrow = 1,
          labels = c("A", "B"))

ggsave("~/github/wnspopgen/figures_tables/persite_both_missANDdepth.png", height = 6, width = 12, dpi=150)
ggsave("~/github/wnspopgen/figures_tables/persite_both_missANDdepth.pdf", height = 6, width = 12, dpi=300)


## cleanup
rm(snpdepth_df, p_depth, p_miss)
