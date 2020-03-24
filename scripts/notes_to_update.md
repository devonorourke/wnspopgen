# Helpful resources
See these links for some additional info that might be useful:
- Want to download a newer protein fasta for all mammals? Subset something? See this post: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc136 
- Recent example of paper for genome pub: https://academic.oup.com/gigascience/article/7/10/giy116/5104371#127032238
- D. Card's wiki about running Maker: https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
  - see these scripts when Darren is referring to merging things from RepeatMasker: https://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html


# Overview

Prior to starting I unmasked all the sites from the updated genome (mylu_hic_rails.fasta file) because they were carried over from the original Myoluc2 genome. The problem was the updated genome included masking information from Myoluc2, but parts were joined/removed in a way that made retaining that information unsuitable. Thus this became the '/scratch/dro49/myluwork/annotation/input_files/mylu_hic_rails_noMasks.fa' file which was used as input for masking. 

## section on RepeatModeler and RepeatMasker
I used that fasta file to run RepeatModeler to generate the custom bat library (fasta format).
  - note about doing this in GenSAS? Or, just rerun in RepeatMaker locally, specifying the library path!
  - summarize commands for both (scripts in /scratch/dro49/myluwork/repeatModels)

In general, Maker v-3.02.beta was used to annotate the geonome by providing:
  - the combined RepeatMasker GFF files (using DFAM mammals and custom library from RepeatModeler), 
  - the protein fasta (update this section with info on it's creation)
  - transcriptome assembly mylu wing (update section with commands on where data came from, how ORP was used)
  - transcriptome assemblies for Myotis brandtii as altest evidence (note same ORP protocol used; mention SRAs for each read set)
  
# Installation woes
The global Maker installation was failing because of Perl dependencies and RepeatMasker not functioning properly. As a result, I had to:
- make an update to the ~/.condarc file to specify to install new env/pkg data to `/scratch/dro49/conda` instead of `$HOME/.conda`. This was because our $HOME directory was filling up with too many packages. 
- Installed Maker. Go to their site (https://www.yandell-lab.org/software/maker.html), enter your info, download the binary, unpack. 
- Set up a Conda environment for Maker installation to work. Do not add in Augustus with this or Maker with this Conda environment - you just want the Perl files; otherwise you'll be knee deep in error messages...  
```
conda create -n maker3env python=3.7
conda activate /scratch/dro49/maker3env
conda install -c bioconda perl-bioperl perl-io-all perl-bit-vector perl-dbd-sqlite perl-inline-c perl-perl-unsafe-signals perl-want perl-forks perl-dbi perl-lwp-simple perl-dbd-pg trf
```
- Now we install Maker _with an active conda environment_:
```
## be sure to have invoked: `conda activate maker3env`
## go to maker's root directory, then to ./src:
./Build install
## this works, but gives this warning:
`Possible precedence issue with control flow operator at Bio/DB/IndexedBase.pm line 805.`
## Apparently the warning/issue is benign? see https://github.com/Ensembl/ensembl-vep/issues/75
```
- I tried installing Augustus locally but kept getting conflicts. What worked for me was to build a copy of Augustus within Maker using the Maker script which downloads a local copy and installs it within Maker (while the same `maker3env` conda environment is active):

```
# ensure 'conda activate maker3env' done...
# within ./src:
./Build augustus
```

Once maker was done, we followed Darren Card's workflow to capture mRNA's with decent confidence to train SNAP and Augustus.
  - see: http://darencard.net/blog/2017-05-16-maker-genome-annotation/

First, collapse all GFF info into single file, with and without sequence info:
```
## GFF with sequence info
/packages/maker/3.01.02-beta/bin/gff3_merge -s \
-d /scratch/dro49/myluwork/annotation/maker_rd1/lu.maker.output/lu_master_datastore_index.log | gzip > \
/scratch/dro49/myluwork/annotation/maker_rd1/mylu_rnd1.all.maker.gff.gz

## GFF without sequence info
/packages/maker/3.01.02-beta/bin/gff3_merge -s -n \
-d /scratch/dro49/myluwork/annotation/maker_rd1/lu.maker.output/lu_master_datastore_index.log > \
/scratch/dro49/myluwork/annotation/maker_rd1/mylu_rnd1.all.maker.noseq.gff
```

We also create a single transcript fasta file:
```
/packages/maker/3.01.02-beta/bin/fasta_merge -o mylu_rnd1.all.maker \
-d /scratch/dro49/myluwork/annotation/maker_rd1/lu.maker.output/lu_master_datastore_index.log
```

I wanted to try to determine what the distribution of AED scores were for this initial annotation round. I used this little script to determine that:
```
awk '$3=="mRNA"' mylu_rnd1.all.maker.noseq.gff | cut -f 9 | cut -d ';' -f 4 | sed 's/_AED=//' | sort | uniq -c > rd1.aed.freqs
```

There are **21,970** gene transcripts reported in this initial run:
  - **21,772** gene transcripts have an AED <= 0.5
  - **14,474** gene transcripts have an AED <= 0.25
  - **10,268** gene transcripts have an AED <= 0.15

Thus most of our mRNA sequences appear to have a confident Annotation Edit Distance (AED) score (lower is better).
  - see Maker documentation for more details: http://www.yandell-lab.org/publications/pdf/maker_current_protocols.pdf

5. Trained SNAP to generate HMM's to put back into Maker.
- first, export 'confident' gene models from Maker, then rename
```
MAKER2ZFF=/packages/maker/3.01.02-beta/bin/maker2zff
FATHOM=/projects/foster_lab/pkgs/SNAP/fathom
FORGE=/projects/foster_lab/pkgs/SNAP/forge
HMMASSEMBLER=/projects/foster_lab/pkgs/SNAP/hmm-assembler.pl
$MAKER2ZFF -x 0.25 -l 50 -d ../../maker_rd1/lu.maker.output/lu_master_datastore_index.log
rename 's/genome/lu_rnd1.zff.length50_aed0.25/g' *
```

- next, gather some stats and validate
```
$FATHOM lu_rnd1.zff.length50_aed0.25.ann lu_rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
$FATHOM lu_rnd1.zff.length50_aed0.25.ann lu_rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1
```
- collect the training sequences and annotations, plus 1000 surrounding bp for training
```
$FATHOM lu_rnd1.zff.length50_aed0.25.ann lu_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
$FATHOM uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
```
- create the training parameters and generate the HMM for Maker
```
mkdir params
cd params
$FORGE ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
$HMMASSEMBLER lu_rnd1.zff.length50_aed0.25 params > lu_rnd1.zff.length50_aed0.25.hmm
```

6. Trained with Augustus to generate HMM to put back into Maker.
- First, pull out just the regions containing mRNA annotations with 1000bp flanks:
```
MKR1GFF=/scratch/dro49/myluwork/annotation/maker_rd1/mylu_rnd1.all.maker.noseq.gff
BEDTOOLS=/home/dro49/.conda/envs/batassembly/bin/bedtools
FASTA=/scratch/dro49/myluwork/annotation/input_files/mylu_hic_rails_noMasks.fa
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $MKR1GFF | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
$BEDTOOLS getfasta -fi $FASTA -bed - -fo lu_rnd1.all.maker.transcripts1000.fasta
```
  - now there are just 598 transcripts we're training on...

- We then run busco using these transcripts, but treat them like they're the genome...

```
# conda create -n busco4env python=3.7
# conda install -c bioconda biopython perl-parallel-forkmanager sepp prodigal
BUSCO=/projects/foster_lab/pkgs/busco4/bin/busco
INFASTA=/scratch/dro49/myluwork/annotation/maker_rd2/lu_rnd1.all.maker.transcripts1000.fasta
BUSCOCONFIG=/projects/foster_lab/pkgs/busco4/config/config.ini
$BUSCO -i $INFASTA --cpu 16 -o busco_rd1 --config $BUSCOCONFIG \
-l /projects/foster_lab/busco_dbs/mammalia_odb10 \
--mode genome --long --augustus_species human --augustus_parameters='--progress=true'
```

The output of this Augustus run needed to be modified because of how our system was set up. You needed to move the data into the /config/species child directory within the Augustus software itself, so that the next round of Maker can point to those HMM parameter files.

7. Reran Maker using the HMM file from SNAP and the HMM config params from Augustus along with the Maker-derived GFFs (est, altest, and protein) as well as the RepeatMasker data. Job on a single node would have taken over 20 days; job using 72 CPUs instead took just under 9 hours. Used up about 20G of RAM.

Repeated steps 3-6, gathering GFF information from second run, pulling out mRNA transcripts, retraining SNAP and Augustus. Briefly, we gathered the collective GFF and transcriptome data first:

```
## rd2 collective GFF data (no sequences)
DIRPATH=/scratch/dro49/myluwork/annotation/maker_rd2/makerRun2
gff3_merge -s -n -d "$DIRPATH"/luM.maker.output/luM_master_datastore_index.log > "$DIRPATH"/mylu_rnd2.all.maker.noseq.gff

## rd2 collective transcripts
fasta_merge -o mylu_rnd2.all.maker -d "$DIRPATH"/lu_master_datastore_index.log

## rd2 AED summary:
awk '$3=="mRNA"' /scratch/dro49/myluwork/annotation/maker_rd2/mylu_rnd2.all.maker.noseq.gff | cut -f 9 | cut -d ';' -f 4 | sed 's/_AED=//' | sort | uniq -c > rd2.aed.freqs
```

Two things to note:
  1. There were more transcripts this time because of the additional ab initio data. The first run produced **21,970** gene transcripts using just our alignment data (est, altest, and protein inferences). This new run produced **23,507** transcripts by Maker, however the ab initio evidence considered **37,482** transcripts. The GFF itself contained **28,646** mRNA strands, so clearly there are some overlapping regions that Maker is combining.
  2. AED distributions were a bit different. The earlier runs suggested about 66% of samples would have an AED less than 0.25; the output of Run2 suggested that there were **28,646** gene transcripts reported in this initial run with AED scores as follows:
    - **22,269** gene transcripts have an AED <= 0.5
    - **15,484** gene transcripts have an AED <= 0.25
    - **7,202** gene transcripts have an AED <= 0.15
  It appears that while we have gained additional transcript evidence in the second round of Maker, the additional transcripts we're gaining aren't as confident overall. In fact, there are fewer _total_ transcripts with an AED under 0.15 that in Round 1, despite having more total transcripts.

We next trained SNAP using data trimmed to include only those transcripts longer than 50 bp and an AED equal or less than 0.25.
```
point to snap.sh script or summarize here
```

We also trained Augustus via BUSCO. First, generate the selected transcripts, filtering for AED score and transcript length:
```
GFFINPUT=/scratch/dro49/myluwork/annotation/maker_rd2/mylu_rnd2.all.maker.noseq.gff
BEDTOOLS=/packages/bedtools/2.26/bin/bedtools
INFASTA=/scratch/dro49/myluwork/annotation/input_files/mylu_hic_rails_noMasks.fa

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $GFFINPUT | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
$BEDTOOLS getfasta -fi $INFASTA -bed - -fo mylu_rnd2.all.maker.transcripts1000.fasta
```

Then use that input `mylu_rnd2.all.maker.transcripts1000.fasta` for BUSCO/Augustus.
```
INFASTA=/scratch/dro49/myluwork/annotation/maker_rd3/mylu_rnd2.all.maker.transcripts1000.fasta
BUSCOCONFIG=/scratch/dro49/conda/envs/busco4env/config/config.ini

busco --in $INFASTA --cpu 28 -o BuscoOut --config $BUSCOCONFIG \
--lineage_dataset /projects/foster_lab/busco_dbs/mammalia_odb10 \
--mode genome --long --augustus_species human --augustus_parameters='--progress=true'
```

Note that the BUSCO summary was:
```
C:53.3%[S:51.0%,D:2.3%],F:11.6%,M:35.1%,n:9226
```

This was a _decrease_ in what we observed in the first version:
```
C:66.0%[S:64.2%,D:1.8%],F:4.2%,M:29.8%,n:9226
```

Why might this be?
1. Maybe  because in the first version of BUSCO/Augustus training, we used _all_ the sequences instead of the AED+length-filtered transcripts. As a result, we would expect fewer total genes to be detected.??
2. Maybe it's because the input transcript file used in the second BUSCO test used a GFF file that included data from _all_ ab initio programs? It looks like SNAP has _double_ the number of predicted transcripts, and this might be _deflating_ the overall AED scores and BUSCO scores? I don't see why having more data would decrease the overall number of missing transcripts though...


Could try a slightly more conservative set of genes to train Augustus species model with?
I ran a small R script to filter for genes with an AED score less than or equal to 0.15 and a length between 5000 to 100000 bp in hopes of selecting suitable (and smaller set of) genes to train Augustus with.
   - specify the R script and upload to Github! (currently in .../annotation)

```
GFFINPUT=/scratch/dro49/myluwork/annotation/maker_rd2/mylu_rd2_filtd_mRNA_forAugustus.gff
BEDFASTA=/packages/bedtools/2.26/bin/fastaFromBed
INFASTA=/scratch/dro49/myluwork/annotation/input_files/mylu_hic_rails_noMasks.fa

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $GFFINPUT | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
$BEDFASTA -fi $INFASTA -bed - -fo mylu_rnd2.filtd.maker.transcripts1000.fasta
```
  - this results in a set of 2,598 sequences to train Augustus with, as opposed to using the full set of about 26k gene models



mylu_rnd2.filtd.maker.transcripts1000.fasta


Prior to running Maker in the next round, we partition the individual sections of the GFF file from the previous Maker run:
```
GFFINPUT=/scratch/dro49/myluwork/annotation/maker_rd1/mylu_rnd1.all.maker.noseq.gff
awk -v OFS="\t" '{ if ($2 == "est2genome") print $0 }' $GFFINPUT > mylu_rnd1.all.maker.est2genome.gff
awk -v OFS="\t" '{ if ($2 == "cdna2genome") print $0 }' $GFFINPUT > mylu_rnd1.all.maker.cdna2genome.gff
awk -v OFS="\t" '{ if ($2 == "protein2genome") print $0 }' $GFFINPUT > mylu_rnd1.all.maker.protein2genome.gff
awk -v OFS="\t" '{ if ($2 ~ "repeat_gff") print $0 }' $GFFINPUT > mylu_rnd1.all.maker.repeats.gff
```
These files can then be imported into the next **maker_opts.ctl** file, replacing the initial alignment and repeat fasta files and gff data. Make sure to turn the "est2genome=" and "protein2genome=" OFF.

The subsequent rounds of maker require a change to the above code because the columns are no longer strictly titled `est2genome` or `cdna2genome`. Instead, these key names are substituted as follows (this example would generate the GFF files used for the third round of Maker, from the second round of Maker's master GFF output):
```
awk -v OFS="\t" '{ if ($2 == "est_gff:est2genome") print $0 }' $GFFINPUT > mylu_rnd2.all.maker.est2genome.gff
awk -v OFS="\t" '{ if ($2 == "altest_gff:cdna2genome") print $0 }' $GFFINPUT > mylu_rnd2.all.maker.cdna2genome.gff
awk -v OFS="\t" '{ if ($2 == "protein_gff:protein2genome") print $0 }' $GFFINPUT > mylu_rnd2.all.maker.protein2genome.gff
```
FYI, the repeat GFF search remains the same:
```
awk -v OFS="\t" '{ if ($2 ~ "repeat_gff") print $0 }' $GFFINPUT > mylu_rnd2.all.maker.repeats.gff
```


/scratch/dro49/myluwork/annotation/maker_rd3/snap_rd2/lu_rnd2.zff.length50_aed0.25.hmm

## Running summary of each round of Maker
We can summarize a few things about the gene models and lengths following each round of Maker:

1. The number of gene models and gene lengths
2. The AED distribution:

I used a custom script to get this information. See `aed_legnth_compviz.R` for more details.


| VarName | maker-1 | maker-2 | maker-2alt | maker-3 |
| --- | --- | --- | --- | --- |
| genes|21970|28646|23819|23232
| meanAED |0.20 |0.34 |0.25|0.25
| medianAED|0.17|0.33|0.20|0.20
| meanExonsBymRNA|5.86|5.68|7.13|7.27
| medianExonsBymRNA|3|4|5|5
| meanmRNAlength|9016.48|16266.19|17499.49|19194.66
| medianmRNAlength|4357.5|4724.5|7910|8165
| meanProteinLength|341.95|305.50|379.42|386.61
| medianProteinLength|248|207|265|272





For transcriptome assembly of Myotis lucifugus ilium tissue
## downloaded RNAseq data from:
https://www.ncbi.nlm.nih.gov/sra/SRX3752331[accn]
## paper was:
https://www.nature.com/articles/s41598-018-33975-x#Sec8
