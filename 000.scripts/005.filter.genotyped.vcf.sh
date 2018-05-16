#!/bin/bash

#PBS -P aquilegia
#PBS -N filt.gen.vcf
#PBS -j oe
#PBS -o /lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/998.logs/colxpat.snp.filter.log
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=16:mem=160GB:local_disk=60gb

cd /lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/002.vcfs

module load Java
module load GATK/3.5-Java-1.8.0_45

# set some variables here

VCFIN=colxpat.05Sept17.genotyped.vcf
VCFOUT=colxpat.05Sept17.genotyped.vcf.FILTERED.BIALLELIC.SNP.vcf


TEMPPATH=$LOCAL_DISK
REF=/lustre/scratch/projects/field_experiments/001.TAIR10.genome/TAIR10_all.fa
GATKPATH=$EBROOTGATK/GenomeAnalysisTK.jar


java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
        -T SelectVariants \
        -R $REF \
        -V $VCFIN \
        -o $VCFOUT \
        -selectType SNP \
        -restrictAllelesTo BIALLELIC \
	 --excludeFiltered \
        --excludeNonVariants \
  	-select "AF>0.4 && AF<0.6" \
       --removeUnusedAlternates 


### output this as a genotype table

java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
        -R $REF\
        -T VariantsToTable\
        -V $VCFOUT\
        -F CHROM -F POS -F REF -F ALT -F FILTER -F NO-CALL\
        -GF GT\
        --allowMissingData\
        -o $VCFOUT.genotype.table.txt


### added 19Sept17
### there seem to be WAY too few SNPs on the top of chromosome four
### are these somehow segregating in the Pat parent?
### in that case, my allele frequency filter would remove them!
### relax this filter, then reassess SNP numbers to see if this could be the case

#java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
#        -T SelectVariants \
#        -R $REF \
#        -V $VCFIN \
#        -o relaxed.filter.$VCFOUT \
#        -selectType SNP \
#        -restrictAllelesTo BIALLELIC \
#        --excludeFiltered \
#        --excludeNonVariants \
#	-select "AF>0.1 && AF<0.9" \
#	--removeUnusedAlternates


### output this as a genotype table

#java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
#        -R $REF\
#        -T VariantsToTable\
#        -V relaxed.filter.$VCFOUT\
#        -F CHROM -F POS -F REF -F ALT -F FILTER -F NO-CALL\
#        -GF GT\
#        --allowMissingData\
#        -o relaxed.filter.$VCFOUT.genotype.table.txt

