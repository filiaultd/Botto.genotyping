#!/bin/bash

#PBS -S /bin/bash
#PBS -P aquilegia
#PBS -N hap.call
#PBS -J 1-192
#PBS -q workq
#PBS -l select=1:ncpus=1:mem=30gb:local_disk=30gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o /lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/998.logs

DIR=/lustre/scratch/users/daniele.filiault/001.Botto_ColxPat

cd $DIR/001.outbam

########################################
# Get the names of the input bam files #
########################################

#ls *realigned.bam > realigned.bam.list.txt
#BAMS=realigned.bam.list.txt
#array_id=$PBS_ARRAY_INDEX
#line=`sed "${array_id}q;d" $BAMS`
#set -- $line
#BAM=$1


BAM=$(ls -1 *realigned.bam | sed "${PBS_ARRAY_INDEX}q;d")



##########################
# Load important modules #
##########################

module load GATK/3.4-0-Java-1.7.0_21
module load Java

# set some variables here


DATE=05Sept17
PREFIX=${BAM:0:6}
TEMPPATH=$LOCAL_DISK
GATKPATH=$EBROOTGATK/GenomeAnalysisTK.jar
REF=/lustre/scratch/projects/field_experiments/001.TAIR10.genome/TAIR10_all.fa
OUTPATH=$DIR/002.vcfs


# run haplotype caller

java -Xmx30G -Djava.io.tmpdir=$TEMPPATH -jar $GATKPATH\
       -T HaplotypeCaller\
       -R $REF\
       -I $BAM\
       --emitRefConfidence GVCF\
       -o $OUTPATH/$PREFIX.$DATE.g.vcf\
       -l INFO
