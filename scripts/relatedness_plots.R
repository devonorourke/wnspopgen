## Script used to generate plots of kinship and relatedness using PLINK outputs
## kinship from 'plink --make-king-table'
## relatedness from 'plink --pca'
## written 12 Jan 2021 by Devon O'Rourke

library(tidyverse)
library(ggpubr)
library(scico)
library(ggrepel)
library(scales)

###############################################################################
########## 1) kinship analysis
###############################################################################

## data import
kinship_df <- read_delim(file="https://github.com/devonorourke/wnspopgen/raw/master/data/PLINK/LUking.kin0.gz",
                         delim="\t", comment = "#",
                         col_names=c("FID1", "ID1", "FID2", "ID2", "NSNP", "HETHET", "IBS0", "KINSHIP"))

kinship_names <- unique(c(unique(kinship_df$ID1), unique(kinship_df$ID2)))
sample_metadata <- read_csv(file="https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/metadata/resolved_sample_metadata.csv") %>% 
  filter(libPrepAlias %in% kinship_names)

## look at distribution of kinship values, in 0.01 (1%) increments
kinship_df %>% 
  mutate(KINSHIP=round(KINSHIP,2)) %>% 
  ggplot(aes(KINSHIP)) +
  geom_histogram(fill="gray50", alpha=0.5, color="black") +
  labs(x="kinship value", y="number of sample pairs") +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

## identify those with kinship...
  ## > 0.354 (indicates duplicate or MZ twin)
  ## between 0.177-0.354 (1st-degree)
    ## between 0.0884-0.177 (2nd degree)... these would absolutely be expected
    ## between 0.0442-0.0884 (3rd degree)... these are likely the majority of our scores

kinship_d0 <- kinship_df %>% filter(KINSHIP > 0.354)
## just one pair: tcap084 and res28... both postWNS group
## tcap084: KG17MA798, Pepperell MA, female Juvenile, MADFW16997, recaptured 8/3/17
## res28: KG17MA452, Pepperell MA, female Adult, MADFW15797, recaptured 6/13/17
  ## likely a mother/daughter (can't be same animal because both captured in same year, the earlier one being the adult, and the later in year being the juvenile!)

kinship_d1 <- kinship_df %>% filter(KINSHIP <= 0.354 ) %>% filter(KINSHIP > 0.177)
nrow(kinship_d1) / nrow(kinship_df) * 100  ## so 0.05% of all data are first degree pairs
  ## 9 pairs that are also likely siblings:
  ## tcap189 & tcap176 (HailesCave, Pepperell - original MA capture in 2016)
  ## tcap190 & tcap141 (HailesCave, Princeton MA - original MA capture in 2013)
  ## sus26 & sus23 (both Stockbridge, VT)
  ## tcap015 & tcap189 (Chester, Hailes cave) **** pre/post mix
  ## tcap072 & tcap009 (Aeolus, Chester)
  ## tcap073 & tcap159 (Aeolus, Pepperell)  *** pre/post mix
  ## tcap075 & tcap183 (Aeolus, Hailes)
  ## tcap075 & tcap056 (Aeolus, Aeolus)

kinship_d2 <- kinship_df %>% filter(KINSHIP <= 0.117 ) %>% filter(KINSHIP > 0.0884)
nrow(kinship_d2) / nrow(kinship_df) * 100  ## so ~ 1% of all pairwise comps are cousins

kinship_d3 <- kinship_df %>% filter(KINSHIP <= 0.0884 ) %>% filter(KINSHIP > 0.0442)
nrow(kinship_d3) / nrow(kinship_df) * 100  ## so ~ 9% of all data are second cousins

kinship_d4toNonNeg <- kinship_df %>% filter(KINSHIP <= 0.0442 ) %>% filter(KINSHIP >= 0)
nrow(kinship_d4toNonNeg) / nrow(kinship_df) * 100  ## so ~ 87% of all data are third or less

### what about negative values? any extreme negatives indicating significant structure?
kinship_neg <- kinship_df %>% filter(KINSHIP < 0)
kinship_neg %>% mutate(KINSHIP = round(KINSHIP, 2)) %>% pull(KINSHIP) %>% table()
  ### nearly all negative values are close to zero, and are likely driven by population structure inherent in the data

### make scatterplot of kinship and proportion of zeroIBS states
ggplot(kinship_df, 
       aes(x=IBS0, y=KINSHIP)) +
  geom_point(alpha=0.5) +
  labs(x="proportion of zero IBS", y="estimated kinship coefficient") +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

ggsave(filename = "~/github/wnspopgen/figures_tables/king_kinship.png", height=15, width=15, units="cm", dpi=150)
ggsave(filename = "~/github/wnspopgen/figures_tables/king_kinship.pdf", height=15, width=15, units="cm", dpi=300)


###############################################################################
########## 2) population stratification
###############################################################################

## visualizing loadings from 'plink --pca'

## import eigenvals per PC
eigval_df <- read_tsv(col_names = 'eigenval',
                      file="https://github.com/devonorourke/wnspopgen/raw/master/data/PLINK/LUpca.eigenval")
eigval_df <- eigval_df %>% 
  mutate(PC = as.numeric(row.names(.)))
  
## make screepot of eigenvals
p_scree <- ggplot(eigval_df, aes(x=PC, y=eigenval/100)) +
  geom_line() + 
  scale_x_continuous(breaks = seq(1,10)) +
  scale_y_continuous(breaks = c(0, 0.02, 0.04), limits = c(0, 0.04)) +
  labs(y="fraction variance explained") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
p_scree
ggsave(file="~/github/wnspopgen/figures_tables/plink_pca_screeplot.png", height=10, width=10, units="cm", dpi=150)
ggsave(file="~/github/wnspopgen/figures_tables/plink_pca_screeplot.pdf", height=10, width=10, units="cm", dpi=300)


## import loadings
loadings_df <- read_tsv(file="https://github.com/devonorourke/wnspopgen/raw/master/data/PLINK/LUpca.eigenvec")
colnames(loadings_df)[1] <- "FID"

## add metadata as needed:
sample_metadata <- read_csv(file="https://raw.githubusercontent.com/devonorourke/wnspopgen/master/data/metadata/resolved_sample_metadata.csv") %>% 
  filter(libPrepAlias %in% loadings_df$IID) %>% 
  select(libPrepAlias, WNSgroup, Location) %>% 
  mutate(WNSgroup = ifelse(WNSgroup=="pre", "SUS", "RES"))

loadings_df <- merge(loadings_df, sample_metadata, by.x='IID', by.y='libPrepAlias')

## plot three figures: PC1:2, PC1:3, PC2:3
p_pc12 <- ggplot() +
  geom_point(data=loadings_df, aes(x=PC1, y=PC2, color=Location, shape=WNSgroup)) +
  theme_bw() +
  coord_fixed() +
  geom_label_repel(data=loadings_df %>% filter(PC2 > 0.6), aes(x=PC1, y=PC2, label=IID), size=3)

p_pc13 <- ggplot() +
  geom_point(data=loadings_df, aes(x=PC1, y=PC3, color=Location, shape=WNSgroup)) +
  theme_bw() +
  coord_fixed() +
  geom_label_repel(data=loadings_df %>% filter(PC3 > 0.4), aes(x=PC1, y=PC3, label=IID), size=3) +
  geom_label_repel(data=loadings_df %>% filter(PC3 < -0.1), aes(x=PC1, y=PC3, label=IID), size=3) +
  geom_label_repel(data=loadings_df %>% filter(PC1 > 0.25), aes(x=PC1, y=PC3, label=IID), size=3)

p_pc23 <- ggplot() +
  geom_point(data=loadings_df, aes(x=PC2, y=PC3, color=Location, shape=WNSgroup)) +
  theme_bw() +
  coord_fixed() +
  geom_label_repel(data=loadings_df %>% filter(PC3 > 0.4), aes(x=PC2, y=PC3, label=IID), size=3) +
  geom_label_repel(data=loadings_df %>% filter(PC3 < -0.1), aes(x=PC2, y=PC3, label=IID), size=3) +
  geom_label_repel(data=loadings_df %>% filter(PC2 > 0.4), aes(x=PC2, y=PC3, label=IID), size=3)

## stitch together:
p_allPC <- ggarrange(p_pc12, p_pc13, p_pc23, nrow=1, ncol=3,
          common.legend = TRUE, labels = c("A", "B", "C", "D"))
p_allPC
ggsave("~/github/wnspopgen/figures_tables/plink_pca_threeOrdinations_wAliases.png",
       height=17, width=25, units="cm", dpi=150)
ggsave("~/github/wnspopgen/figures_tables/plink_pca_threeOrdinations_wAliases.pdf",
       height=17, width=25, units="cm", dpi=300)

### also replot these so that we don't have any label names:
## plot three figures: PC1:2, PC1:3, PC2:3
p_pc12_nl <- ggplot() +
  geom_point(data=loadings_df, aes(x=PC1, y=PC2, color=Location, shape=WNSgroup), size=2.5) +
  theme_bw() +
  coord_fixed()

p_pc13_nl <- ggplot() +
  geom_point(data=loadings_df, aes(x=PC1, y=PC3, color=Location, shape=WNSgroup), size=2.5) +
  theme_bw() +
  coord_fixed()

p_pc23_nl <- ggplot() +
  geom_point(data=loadings_df, aes(x=PC2, y=PC3, color=Location, shape=WNSgroup), size=2.5) +
  theme_bw() +
  coord_fixed()

ggarrange(p_pc12_nl, p_pc13_nl, p_pc23_nl, nrow=1, ncol=3,
          common.legend = TRUE, labels = c("A", "B", "C", "D"))

ggsave("~/github/wnspopgen/figures_tables/plink_pca_threeOrdinations_noLabels.png",
       height=17, width=25, units="cm", dpi=150)
ggsave("~/github/wnspopgen/figures_tables/plink_pca_threeOrdinations_noLabels.pdf",
       height=17, width=25, units="cm", dpi=300)
  ## will edit to pdf to align x axes in all 3 panels along bottom and move legend into space created in top right

###############################################################################
########## 3) Fst sliding windows
###############################################################################

## generated from 'vcftools --fst' command of RES/SUS populations
fstwindows_df <- read_delim(file="https://github.com/devonorourke/wnspopgen/raw/master/data/vcftools/RES_SUS_fst_100kbWindow_25kbStep.windowed.weir.fst.gz",
                         delim="\t", col_names=TRUE)

## calculate the median weighted FST per chromosome, ...
  ### rather than sticking with just a single global value when plotting the horizontal line
redlineMedianFst_df <- fstwindows_df %>% group_by(CHROM) %>% summarise(medianFst = median(WEIGHTED_FST))

## merge together into single dataframe for plotting
fstwindows_df <- merge(fstwindows_df, redlineMedianFst_df)

## order the plot facets properly
scaffOrder <- paste0("scaffold", (c(seq(1,3), seq(5,22))))
fstwindows_df$CHROM <- factor(fstwindows_df$CHROM,
                              levels = scaffOrder)

## plot
ggplot() +
  geom_point(data = fstwindows_df,
             aes(x=(BIN_START/1000000), y=MEAN_FST),# scale x axis in millins of bp
             alpha=0.5) +
  geom_hline(data = fstwindows_df,
             aes(yintercept = medianFst), 
             color="red", alpha=0.5) +
  labs(x="position (Mb)", y=expression(F[ST])) +
  facet_wrap(~CHROM, ncol=3, scales = "free_x") +
  theme_bw()

ggsave("~/github/wnspopgen/figures_tables/fst_100kbWindow_25kbStep.png",
       height = 28, width = 22, units = "cm", dpi=150)
ggsave("~/github/wnspopgen/figures_tables/fst_100kbWindow_25kbStep.pdf",
       height = 28, width = 22, units = "cm", dpi=300)

### unused:
# if so desired, could also grab the topN windows with outliers as follows:

# top10fst_perChrom <- fstwindows_df %>% 
#   group_by(CHROM) %>% 
#   slice_max(order_by = WEIGHTED_FST, n = 10)