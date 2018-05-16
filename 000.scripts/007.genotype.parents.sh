#!/bin/bash

### shell script for genotyping combined GVCFs files produced by haplotype caller

#PBS -P aquilegia
#PBS -N CnP.gen.gvcf
#PBS -j oe
#PBS -o /lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/998.logs/colxpat.genotype.g.vcfs.log
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=16:mem=60GB:local_disk=60gb


DIR=/lustre/scratch/users/daniele.filiault/001.Botto_ColxPat

cd $DIR/003.ilka.col.pat.vcf

module load GATK/3.4-0-Java-1.7.0_21

# set some variables here

PREFIX=CnP.parents
DATE=08Sept17
GATKPATH=$EBROOTGATK/GenomeAnalysisTK.jar
REF=/lustre/scratch/projects/field_experiments/001.TAIR10.genome/TAIR10_all.fa
TEMPPATH=$LOCAL_DISK

# run haplotype caller to generate common vcf with both parents
## individual g.vcfs are from Ilka
#java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH \
#   -T GenotypeGVCFs \
#   -R $REF \
#   --variant $PREFIX.g.vcf.files.list \
#   -o $PREFIX.$DATE.genotyped.vcf

# filter this for positions that are called in both parents (AN=4) and for SNP variants only

java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
        -T SelectVariants \
        -R $REF \
        -V $PREFIX.$DATE.genotyped.vcf\
        -o $PREFIX.$DATE.filtered.vcf \
        -selectType SNP \
        -restrictAllelesTo BIALLELIC \
        --excludeFiltered \
        --excludeNonVariants \
        -select "AN==4" \
        --removeUnusedAlternates 

# output as genotype table

java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
        -R $REF\
        -T VariantsToTable\
        -V $PREFIX.$DATE.filtered.vcf\
        -F CHROM -F POS -F REF -F ALT -F FILTER -F NO-CALL\
        -GF GT \
        --allowMissingData \
        -o $PREFIX.$DATE.genotype.table.txt

