### parse Col and Pat parent SNP set
### filter line gts for SNPs good in parental data set
### quality control and description of SNP numbers and distribution
### DLF 08Sept17

setwd("/lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/003.ilka.col.pat.vcf/")


#### reads in SNPs called in parental lines
dat <- read.table("CnP.parents.08Sept17.genotype.table.txt",header=TRUE)

dim(dat)  #524999 variant sites total
dat <- dat[,c(1:4,7:8)]
colnames(dat)[5:6]<-c("col.gt", "pat.gt")

### change letter genotypes to binary
### 0=ref.hom,0.5=het,1=alt.hom
### ultimately need to remove het gts

bin.gt <- matrix(NA, nrow=nrow(dat),ncol=2)

### this acts on each row of the dat variable
bp.to.bin <- function(up.dat){
	ref <- up.dat[,3]
	alt <- up.dat[,4]
	hom.ref <- paste(ref,ref, sep="/")
	hom.alt <- paste(alt, alt, sep="/")
	het.a <- paste(alt, ref, sep="/")
	het.b <- paste(ref, alt, sep="/")
	ind <- cbind(c(hom.ref, hom.alt,het.a, het.b),c(0,1,0.5,0.5))
	cg <- ind[ind[,1]==up.dat$col.gt,2]
	pg <- ind[ind[,1]==up.dat$pat.gt,2]
	out.gt <- c(cg, pg)
	return(out.gt)
}

for (up in 1:nrow(dat)){
	up.dat <- dat[up,]
	up.out <- bp.to.bin(up.dat)
	bin.gt[up,] <- up.out
}

bin.gt <- cbind(dat[,1:2],bin.gt)
colnames(bin.gt)[3:4]<-c("col","pat")
bin.gt$combo <- paste(bin.gt$col, bin.gt$pat, sep="_")

###   0_0.5     0_1   0.5_0 0.5_0.5   0.5_1     1_0   1_0.5     1_1 
  73974  444943    2739    1191    1602     120      32     398 

### keep only SNPs that are 0_1 (the expected pattern here)
bin.gt.g <- bin.gt[bin.gt$combo=="0_1",]  ##444943 SNPs.  Not bad.
bin.gt.s <- split(bin.gt.g, bin.gt.g$CHROM)


######## how many SNPs per line?

# reads in line SNP data
line.dat <- read.table("/lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/002.vcfs/colxpat.05Sept17.genotyped.vcf.FILTERED.BIALLELIC.SNP.vcf.genotype.table.txt", header=TRUE, comment.char="",stringsAsFactors=FALSE)  ### 325106 SNPs, 193 lines

sample.index <- cbind(7:193,colnames(line.dat)[7:193])
colnames(line.dat)[7:193] <- paste("S",7:193,sep="_")

# turn these into 0/1 coding

sample.gt <- matrix(NA,nrow=nrow(line.dat),ncol=ncol(line.dat)-6)
bp.to.bin2 <- function(up.dat){
	ref <- up.dat[,3]
	alt <- up.dat[,4]
	hom.ref <- paste(ref,ref, sep="/")
	hom.alt <- paste(alt, alt, sep="/")
	het.a <- paste(alt, ref, sep="/")
	het.b <- paste(ref, alt, sep="/")
	ind <- cbind(c(hom.ref, hom.alt,het.a, het.b),c(0,1,0.5,0.5))
	ind <-rbind(ind,c("./.",NA))
	out.gt <- sapply(7:193, function(x){
		up.gt <- up.dat[1,x]
		xg <- ind[ind[,1]==up.gt,2]
		return(xg)
	})
	return(out.gt)
}

for(up in 1:nrow(line.dat)){
	#print(up)
	up.dat <- line.dat[up,]
	up.gt <- bp.to.bin2(up.dat)
	sample.gt[up,] <- up.gt
}

sample.gt <- cbind(line.dat[,1:2], sample.gt)
sample.gt.s <- split(sample.gt,sample.gt$CHROM)

### remove SNPs not in "good" (i.e. parental) SNPs as determined above

filter.gt.s <- as.list(1:5)
for(up in 1:5){
	up.sg <- sample.gt.s[[up]]
	up.snps <- bin.gt.s[[up]]
	out.snps <- up.sg[up.sg$POS%in%up.snps$POS,]
	filter.gt.s[[up]] <- out.snps
}

filter.gt <- do.call(rbind,filter.gt.s) ###325106 SNPs

sample.t <- apply(filter.gt[,3:ncol(filter.gt)], 2, table)
sample.t <- do.call(rbind,sample.t)

sample.total <- apply(sample.t,1,sum)
sample.p <- sample.t/sample.total

sample.total[sample.total<30000] #8 samples

table(filter.gt$CHROM)
### chromosome four has MANY fewer SNPs!  why is this?
### the problem is not in the parental genotyping (lapply(bin.gt.s,nrow))
### it is in the line.dat variable (the number of SNPs called per sample that gets entered, ie the raw data for this script)
### very strange!!!  need to explore why this is failing...

#let's see if there are obvious patterns to the data along the chromosome

pdf("missing.calls.by.pos.pdf", width=10, height=8)
par(mfrow=c(2,3))
for (up.chr in unique(filter.gt$CHROM)){
	print(up.chr)
	up.dat <- filter.gt[filter.gt$CHROM==up.chr,]
	up.na <- apply(up.dat,1,function(x){sum(is.na(x))})
	plot(up.dat$POS, (up.na)/187, main=up.chr, xlab="position",ylab="proportion of samples not called", ylim=c(0,1))
	}
dev.off()

### OK, so the parentals were called with all samples lumped together - are SNPs more het here?

lapply(bin.gt.s, function(x){table(x$pat)})

### so, are there patterns in SNP calling in individual samples on chromosome four?  i.e. more positions het,nocalls, per sample on chromosome four

by.chr.gts <- lapply(filter.gt.s, function(x){
	sample.t <- apply(x[,3:ncol(filter.gt)], 2, table)
	sample.t <- do.call(rbind,sample.t)
	sample.total <- apply(sample.t,1,sum)
	sample.p <- sample.t/sample.total
	return(sample.p)
})

pdf("line.gts.by.chr.pdf",width=8, height=4)
for(up in 1:length(by.chr.gts)){
	up.dat <- by.chr.gts[[up]]
	up.dat <- up.dat[order(up.dat[,1]),]
	barplot(t(up.dat), main=up)
	
}
dev.off()

#### yes, these distributions are clearly weird.
#### I think the thing to do is to "unfilter" the samples, do the filtering with VariantFiltration properly, and then look to see how the individual Pat samples behave.  Or if there is a pattern to the filtered samples.  If something was het in Pat to start with, then it would only have an allele frequency of 25% here and I would have removed it in my filtration.  So if there is some big SV on Chromosome four that was segregating in the original Pat parent, essentially causing a lack of recombination on one arm of Chr4, this could be the cause...

### refiltered the lines SNPs data, using cutoff proportions of 0.1 and 0.9
### repeat above to see if things are better

relax.dat <- read.table("/lustre/scratch/users/daniele.filiault/001.Botto_ColxPat/002.vcfs/relaxed.filter.colxpat.05Sept17.genotyped.vcf.FILTERED.BIALLELIC.SNP.vcf.genotype.table.txt", header=TRUE, comment.char="",stringsAsFactors=FALSE)  ### 325106 SNPs, 193 lines
colnames(relax.dat)[7:193] <- paste("S",7:193,sep="_")



rsnp.t <- table(relax.dat$CHROM)  ### yes, this is much better!!!
rsnp.call.t <- table(relax.dat$CHROM, relax.dat$NO.CALL)
### so, indeed, there is something segregating here!
### check gts of these SNPs.  More hets?
chr.cols <- c("black","blue","green","red","purple")
pdf("no.calls.by.chr.relaxed.pdf")
plot(rsnp.call.t[1,],type="l")
for(up in 2:5){
points(rsnp.call.t[up,],type="l", col=chr.cols[up])}
dev.off()
### still slightly skewed, but much better!


#sample.index <- cbind(7:193,colnames(line.dat)[7:193])
#colnames(line.dat)[7:193] <- paste("S",7:193,sep="_")
# filter for SNPs in the good SNP set

r.sample.gt <- matrix(NA,nrow=nrow(relax.dat),ncol=ncol(relax.dat)-6)

for(up in 1:nrow(relax.dat)){
        #print(up)
        up.dat <- relax.dat[up,]
        up.gt <- bp.to.bin2(up.dat)
        r.sample.gt[up,] <- up.gt
}

 

#### output "bad"ly segregating SNP positions for testing

homa <- apply(r.sample.gt,1, function(x){
	x.out <- x[x==1]
	x.out <- x.out[is.na(x.out)==FALSE]
	x.l <- length(x.out)
	return(x.l)
})

homr <- apply(r.sample.gt,1, function(x){
	x.out <- x[x==0]
	x.out <- x.out[is.na(x.out)==FALSE]
	x.l <- length(x.out)
	return(x.l)
})

het <- apply(r.sample.gt,1, function(x){
	x.out <- x[x==0.5]
	x.out <- x.out[is.na(x.out)==FALSE]
	x.l <- length(x.out)
	return(x.l)
})

total <- homr + homa + (0.5*het)
p.r <- (homr + (0.5*het))/total
p.a <- (homa + (0.5*het))/total
p.ar <- p.a/p.r:w

par(mfcol=1,3)
hist(p.r)
hist(p.a)
hist(p.ar)
dev.off()

p.dat <- cbind(relax.dat[,1:2], p.a)
p.dat.s <- split(p.dat, p.dat$CHROM)

pdf("alt.allele.proportion.by.chromsome.pdf")
for(up in 1:5){
	up.dat <- p.dat.s[[up]]
	up.chr <- up.dat[1,1]
	hist(up.dat[,3], xlab="proportion of alternate alleles", main=up.chr)
	}
dev.off()

pdf("alt.alleles.along.chromosomes.pdf", width=8, height=4)
for(up in 1:5){
	up.dat <- p.dat.s[[up]]
	up.chr <- up.dat[1,1]

	colramp = colorRampPalette(c('white', 'blue', 'green', 'yellow', 'red'))
	smoothScatter(up.dat[,2], up.dat[,3], colramp=colramp, main=up.chr, xlab="position",ylab="proportion of alternate allele")
	}
dev.off()

write.csv(p.dat.s[[4]],"alt.allele.prop.chromsome.four.csv", quote=FALSE, row.names=FALSE)

################################################
### output gt variable
### for running HMM
################################################

save(filter.gt.s, file="filter.gt.s.Rdata")


