#### need to send an index of which samples are which to Javier
#### 06Nov17

setwd("/Volumes/field_experiments/001.Botto_ColxPat/007.sample.index")


ss4 <- read.csv("Botte_F6ColxPat1_INDEX_S4 - Sheet1.csv")
ss3 <- read.csv("Botte_F6ColxPat2_INDEX_S3 - Sheet1.csv")


s.names <- read.table("/Volumes/field_experiments/001.Botto_ColxPat/002.vcfs/colxpat.05Sept17.genotyped.vcf.FILTERED.BIALLELIC.SNP.vcf.genotype.table.txt", nrows=1,comment.char="",stringsAsFactors=FALSE)

s.names <- s.names[,-c(1:6)]
s.names <- gsub("CAYP9ANXX_3#","",s.names)
s.list <- sapply(s.names, function(x){strsplit(x,"_")})
s.short <- lapply(s.list, function(x){paste(x[2],x[1],sep="#")})
s.short <- unlist(s.short)

### these are the names, in order, of the genotypes that I have generated
### save to use to label the data frame
save(s.short, file="s.short.Rdata")


#### generate a file that links these to the F6 sample names

sample.index <- rbind(ss3, ss4)
sample.index$sample.description <- paste("ColPat",sample.index$sample.description,sep="_")
sample.index <- sample.index[,c(5,4)]
names(sample.index) <- c("gt.name","F6.name")
write.csv(sample.index, file="ColPat.sample.index.csv",quote=FALSE,row.names=FALSE)


