### description of scripts in 000.scripts folder

001.  Renames realigned.bam files that are output from standard mapping protocol (bwa, GATK, etc).

002.  An array job to generate gVCF files for each bam file

003.  Makes a list of all gVCF files generated in 002

004.  Genotypes list from 003 with GATK Haplotype Caller

005.  Filters genotypes from 004 for 0.4<= allele frequency <=0.6 and outputs a genotype table
	(colxpat.05Sept17.genotyped.vcf.FILTERED.BIALLELIC.SNP.vcf.genotype.table.txt)
	Script also does another implementation with more leniant filtration for troubleshooting Chr4 problems

006.  Gets coverage per sample for individual bams from 001

007.  Genotypes SNPs in the parents alone and outputs a genotype table of SNPs
	(CnP.parents.08Sept17.genotype.table)

008.  Filters parental SNPs from 007 to those that are 0/1 or 1/0 Col/Pat.
	Filters the RIL snps to only these positions.  (filter.gt.s.Rdata)
	There are also some analyses troublshooting Chr4 problems here

009.  Genotyping of variants.  First removes potential supurious hets (het in Pat, overlap TEs)
	Then does HMM using HMM in R.
	Then does some "smoothing" in windows.
	Here 100kb windows, requiring 10snps per window, and at least 50% of snps in window called as one genotype to include window.
	(ColPat.genotypes.100kb.window.csv)

07. generates an index to link our sequencing sample names to the Botto lab sample numbers
