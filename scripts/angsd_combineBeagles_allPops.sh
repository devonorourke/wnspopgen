## Example shell script used to combine the outputs from ANGSD genotype likelihood jobs for individual populations
## This example used the two pre/post populations of Myotis lucifugs

## get the SNPsite locations from each LU post/pre populations
zcat LUpost.beagle.gz | awk 'NR>1' | cut -f 1-3 | tr '\t' ':' | sort > LUpost.beagle.SNPsites.tmp
zcat LUpre.beagle.gz | awk 'NR>1' | cut -f 1-3 | tr '\t' ':' | sort > LUpre.beagle.SNPsites.tmp

## find the common sites matching to each
## note we're keeping only those SNP sites where the major and minor allele are identical in both populations
comm -12 LUpost.beagle.SNPsites.tmp LUpre.beagle.SNPsites.tmp | cut -f 1 -d ':' > LUcommon.beagle.SNPsites.tmp

## create a new merged beagle file that includes these sites only
## merging these beagle files together:
printf 'marker allele1 allele2' | tr ' ' '\t' > tmpfile1
cut -f 7 -d '/' LUpost.list | cut -f 1 -d '_' | awk '{for(i=0;i<3;i++)print}' | tr '\n' '\t' | sed 's/[[:space:]]\+$//' > tmpfile2
cut -f 7 -d '/' LUpre.list | cut -f 1 -d '_' | awk '{for(i=0;i<3;i++)print}' | tr '\n' '\t' | sed 's/[[:space:]]\+$//' > tmpfile3
paste tmpfile{1,2,3} > LUcommon.beagle

awk 'FNR==NR {arr[$0];next} $1 in arr' LUcommon.beagle.SNPsites.tmp <(zcat LUpost.beagle.gz | awk 'NR>1' | sort) > tmpfile4

awk 'FNR==NR {arr[$0];next} $1 in arr' LUcommon.beagle.SNPsites.tmp <(zcat LUpre.beagle.gz | awk 'NR>1' | sort) | cut -f 4- > tmpfile5

paste tmpfile{4,5} >> LUcommon.beagle

gzip LUcommon.beagle

#rm tmpfile{4,5}

## create a second merged beagle file that includes just sites that were present in the earlier VCF files
## filter the common SNP sites above present in the earlier VCF list
comm -12 <(sort LUall_regions_undscor.file) LUcommon.beagle.SNPsites.tmp > LUall_regions.beagle.SNPsites.tmp
paste tmpfile{1,2,3} > LUall_regions.beagle

awk 'FNR==NR {arr[$0];next} $1 in arr' LUall_regions.beagle.SNPsites.tmp <(zcat LUpost.beagle.gz | awk 'NR>1' | sort) > tmpfile6

awk 'FNR==NR {arr[$0];next} $1 in arr' LUall_regions.beagle.SNPsites.tmp <(zcat LUpre.beagle.gz | awk 'NR>1' | sort) | cut -f 4- > tmpfile7

paste tmpfile{6,7} >> LUall_regions.beagle

gzip LUall_regions.beagle

rm tmpfile{1,2,3,6,7}
rm *SNPsites.tmp
