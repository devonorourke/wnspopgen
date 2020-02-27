## example script to identify relatedness among Myotis lucifugus samples
## using VCF file of all SNP sites common to both populations

NGSRELATE=/projects/foster_lab/pkgs/ngsRelate/ngsRelate
VCFFILE=/scratch/dro49/myluwork/popgen/vcfs/pops/mylu/LUcommon.vcf.gz

$NGSRELATE -h $VCFFILE -O mylu_allCommon.res -T GT -c 1 -p 22
