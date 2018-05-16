#!/bin/bash

### shell script for genotyping combined GVCFs files produced by haplotype caller

#PBS -P aquilegia
#PBS -N col.pat.gen.gvcf
#PBS -j oe
#PBS -o /lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/998.logs/colxpat.genotype.g.vcfs.log
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=16:mem=60GB:local_disk=60gb


DIR=/lustre/scratch/users/daniele.filiault/001.Botto_ColxPat

cd $DIR/002.vcfs

module load GATK/3.4-0-Java-1.7.0_21

# set some variables here

PREFIX=colxpat
DATE=05Sept17
GATKPATH=$EBROOTGATK/GenomeAnalysisTK.jar
REF=/lustre/scratch/projects/field_experiments/001.TAIR10.genome/TAIR10_all.fa
TEMPPATH=$LOCAL_DISK

# run haplotype caller
java -Xmx60G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH \
   -T GenotypeGVCFs \
   -R $REF \
   --variant $PREFIX.g.vcf.files.list \
   -o $PREFIX.$DATE.genotyped.vcf

