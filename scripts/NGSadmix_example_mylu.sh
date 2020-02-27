## example script for analyzing the M. lucifugus SNP sites common among both populations

NGSADMIX=/projects/foster_lab/pkgs/angsd/misc/NGSadmix
BEAGLEFILE=/scratch/dro49/myluwork/popgen/angsd/mylu/LUcommon.beagle.gz
OUTDIR=$(pwd)

$NGSADMIX \
-P $SLURM_CPUS_PER_TASK \   ## we use a Slurm program manager to submit jobs; change this to an integer for # thread you want to use in your own run
-likes $BEAGLEFILE \
-o LUcommon_k2 \
-K 2 -minMaf 0.05 -seed 1

$NGSADMIX \
-P $SLURM_CPUS_PER_TASK \
-likes $BEAGLEFILE \
-o LUcommon_k3 \
-K 3 -minMaf 0.05 -seed 1

$NGSADMIX \
-P $SLURM_CPUS_PER_TASK \
-likes $BEAGLEFILE \
-o LUcommon_k4 \
-K 4 -minMaf 0.05 -seed 1
