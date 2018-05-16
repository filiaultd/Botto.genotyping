#### HMM for genotyping of ColxPat F6s
#### DLF 19 Sept 17


library(HMM)
setwd("/Volumes/field_experiments/001.Botto_ColxPat/006.hmm.genotyping")


### load SNP data
load("filter.gt.s.Rdata")  ### filter.gt.s is filtered genotypes for all samples, split by chromosome

### pull out one sample, one chromosome to work with

par(mfcol=c(1,5))

for (up in 1:5){
up.chr <- up
#up.dat <- filter.gt.s[[up.chr]]
#up.dat <- ugh.clean.s[[up.chr]]
up.dat <- snps.no.te.s[[up.chr]]
up.line <- 3
up.snps <- up.dat[,c(1:2,up.line+2)]
up.snps <- up.snps[is.na(up.snps[,3])==FALSE,]

## initialize HMM
gt.hmm <- initHMM(c("ref", "het", "alt"), c("0","0.5","1"), transProbs=matrix(c(.9989,.001,.0001,.0001,.9998,.0001,.0001,.001,.9989),3), emissionProbs=matrix(c(.599,.4,.001,.1,.8,.1,.001,.4,.599),3))
#obs <- c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,0.5,0.5,0.5,1,0.5)
#obs <- rep("1",100)
obs <- as.character(up.snps[,3])

gt.vit <- viterbi(gt.hmm, obs)


out.dat <- cbind(up.snps, gt.vit)

gt.vit <- as.factor(gt.vit)
gt.vit <- as.numeric(gt.vit)
plot(up.snps$POS, gt.vit, col="blue", ylim=c(0.5,3))
points(up.snps$POS, abs(4-as.numeric(up.snps[,3]))-0.1, col="red")
legend(x="bottom", col=c("blue", "red"), legend=c("genotype","snp"), pch=19,)

#gt.test <- (gt.vit-1)/2
##gt.t.v <- viterbi(gt.hmm, as.character(gt.test))
#gt.t.v <- as.numeric(as.factor(gt.t.v))
#quartz()
#plot(up.snps$POS, gt.t.v, col="blue", ylim=c(0.5,3))
#points(up.snps$POS, abs(4-as.numeric(up.snps[,3]))-0.1, col="red")
#legend(x="bottom", col=c("blue", "red"), legend=c("genotype","snp"), pch=19,)
}

####################################################################
#### do some additional filtering to remove SNPs with "too much" heterozygosity
#### in all samples, plus in 184,186,187 (which are pure Pat samples)

### filter.gt.s[[1]] -> ugh

ugh <- do.call(rbind, filter.gt.s)

an <- ugh==1
an <- apply(an, 1, function(x){sum(x, na.rm=TRUE)})
hn <- ugh==0.5
hn <- apply(hn, 1, function(x){sum(x, na.rm=TRUE)})
rn <- ugh==0
rn <- apply(rn, 1, function(x){sum(x, na.rm=TRUE)})
gt.n <- cbind(rn, hn, an)
sn <- apply(gt.n,1, sum)
gt.n <- cbind(gt.n, sn)
gt.p <-gt.n[,1:3]/sn
par(mfcol=c(1,3))
for(up in 1:3){hist(gt.p[,up])} ## so the distributions here look OK, not a huge skew in the number of het gts per SNP.  Nonetheless, one could filter on this?

### perhaps an easier/more straightforward task is to remove any SNP that is het in Pat

pat.s <- c(184,186,187)
pat.dat <- ugh[,colnames(ugh)%in%pat.s]
pat.dat <- cbind(ugh[,1:2], pat.dat)
pat.dat[,3]<-as.numeric(as.character(pat.dat[,3]))
pat.dat[,4]<-as.numeric(as.character(pat.dat[,4]))
pat.dat[,5]<-as.numeric(as.character(pat.dat[,5]))

### get any snp that is het here
pat.h <- apply(pat.dat, 1,function(x){0.5%in%x})
pat.hs <- pat.dat[pat.h==TRUE,]

### remove from SNP set
ugh.clean <- ugh[pat.h==FALSE,]
ugh.clean.s <- split(ugh.clean, ugh.clean$CHROM)
### try the HMM for this set
### it actually helps!


### one more thing to try - remove any SNP that maps to a TE position
### got filtered TAIR10 from Robin
te.pos <- read.table("../Tair10.triplemasked.bed")
te.pos.s <- split(te.pos, te.pos[,1])


snps.no.te.s <- as.list(1:5)
for(up.chr in 1:5){
	te.up <- te.pos.s[[up.chr]]
	snps.up <- ugh.clean.s[[up.chr]]
	filt.int <- as.list(1:nrow(te.up))
	for(up in 1:nrow(te.up)){
		print(up)
		out.pos <- (te.up[up,2]+1):te.up[up,3]
		filt.int[[up]] <- out.pos
		}
	filt.int <- unlist(filt.int)
	snps.filt <- snps.up[snps.up$POS%in%filt.int==FALSE,]
	snps.no.te.s[[up.chr]]<-snps.filt
}

### now try this with the HMM...

#### write an HMM function to do all lines
#snps = snps.no.te.s
#line.number = 3

hmm.gt <- function(line.no, snps){
	out.gts.s <- as.list(1:5)
	for (up in 1:5){
		#print(up)
		
			return(out.gts.s)
	}


#ugh <- hmm.gt(line.no=3,snps=snps.no.te.s)

#ugh2 <- hmm.gt(line.no=67,snps=snps.no.te.s)

#### This works.  Run for all lines
#### doesn't work for line 61!!!  not enough SNPs.


all.gts <- as.list(:(ncol1(snps.no.te.s[[1]])-2))
for (up.l in c(1:60,62:(ncol(snps.no.te.s[[1]])-2))){
	print(up.l)
	l.gts <- hmm.gt(line.no=up.l, snps=snps.no.te.s)
	all.gts[[up.l]] <- l.gts
}


########### these need to be genotyped in windows
########### easy window-based filtering
########### variables : window size, prop of SNPs that need to be one gt to call window, min #SNPs in window

window.size <- 100000
min.snp.no <- 10
prop.gt <- .50

#### step one: make windows
index <- read.table("/Volumes/field_experiments/001.common.reference.files/001.TAIR10.genome/TAIR10_all.fa.fai")
wins <- as.list(1:5)
for(up in 1:5){
	up.m <- index[up,2]
	up.start <- seq(1,up.m,window.size)
	up.stop <- up.start+(window.size-1)
	up.stop[length(up.stop)]<-up.m
	up.win <- cbind(up.start,up.stop)
	wins[[up]] <- up.win
}

#### step two: go through each sample, each chromosome
#### each window, call gt of window

### don't use sample 61!  Not enough SNPs
line.gts <- c(1:60,62:length(all.gts))

line.gts.out <- lapply(wins, function(x){
	dat.m <- matrix(data=NA,nrow=nrow(x), ncol=187)
	return(dat.m)
})



for(up.s in line.gts){
	print(up.s)
	up.gts <- all.gts[[up.s]]
	for(up.c in 1:5){
		print(up.c)
		uu.g <- up.gts[[up.c]]
		uu.w <- wins[[up.c]]
		all.w.gts <- rep(NA,nrow(uu.w))
		for(up.w in 1:nrow(uu.w)){
			up.pos <- seq(uu.w[up.w,1], uu.w[up.w,2])
			uuu.g <- uu.g[uu.g$POS%in%up.pos,]
			if(nrow(uuu.g)<=min.snp.no){w.gt <- NA}
			else{
				uuu.table <- table(uuu.g$gt.vit)
				uuu.table <- uuu.table[order(uuu.table,decreasing=TRUE)]
				win.total <- sum(uuu.table)
				win.prop <- uuu.table[1]/win.total
				if(win.prop>=prop.gt){w.gt <- names(win.prop)} else {w.gt <- NA}
				}
			all.w.gts[up.w]<-w.gt
			}
	line.gts.out[[up.c]][,up.s]<-all.w.gts	
	}
}

##### now make this into one file to send to Javier and Daniel

### add sample number to each column
### add chromosome, window position to each
### concatenate chromosomes

### index of which numbers are which samples
### from 07.generate.index.R script

load("/Volumes/field_experiments/001.Botto_ColxPat/007.sample.index/s.short.Rdata")

test <- lapply(line.gts.out, function(x){
	colnames(x)<-s.short
	return(x)
	})

test <- do.call(rbind, test)

### get coordinates of these windows
win.co <- wins
win.pos <- lapply(as.list(1:5),function(x){
	up.dat <- win.co[[x]]
	up.chr <- rep(x,nrow(up.dat))
	up.dat <- cbind(up.chr, up.dat)
	return(up.dat)
})
win.pos <- do.call(rbind,win.pos)
colnames(win.pos)<-c("chr", "window.start","window.stop")

window.data <- cbind(win.pos, test)
write.csv(window.data, file="ColPat.genotypes.100kb.window.csv",quote=FALSE,row.names=FALSE)





###################################################################################
##### testing here
######################


hmm = initHMM(c("X","Y"),c("a","b","c"))
simHMM(hmm, 100) -> test.dat

viterbi(hmm,test.dat$observation )

forward(hmm, test.dat$observation)


a = sample(c(rep("a",100),rep("c",300)))
b = sample(c(rep("a",300),rep("c",100)))
observation = c(a,b)

viterbi(hmm, observation )


# Initialise HMM
hmm = initHMM(c("A","B"), c("L","R"), transProbs=matrix(c(.6,.4,.4,.6),2), emissionProbs=matrix(c(.6,.4,.4,.6),2)) 

print(hmm)

# Sequence of observations
observations = c("L","L","R","R")

# Calculate Viterbi path
viterbi = viterbi(hmm,observations)
print(viterbi)