## example script for analyzing the M. lucifugus SNP sites common among both populations
## note that the same script was applied for K=2, K=3, and K=4 to generate the data used in the figure produced

NGSADMIX=/projects/foster_lab/pkgs/angsd/misc/NGSadmix
BEAGLEFILE=/scratch/dro49/myluwork/popgen/angsd/mylu/LUcommon.beagle.gz
OUTDIR=$(pwd)

$NGSADMIX \
-P $SLURM_CPUS_PER_TASK \
-likes $BEAGLEFILE \
-o LUcommon_k2 \
-K 2 -minMaf 0.05 -seed 1
