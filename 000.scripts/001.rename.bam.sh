### take long realigned bam names and shorten to 6 character informative names to make everything less crazy!


UPLIST=realigned.bams.txt

while read f; do
	echo $f
	MIN_F=${f:12:6}
	echo $MIN_F
	#cp $f $MIN_F.realigned.bam
	#cp $f.bai $MIN_F.realigned.bam.bai
	rm $f
	rm $f.bai

done < $UPLIST

ls *realigned.bam > short.realigned.bams.txt

### now ready to go into generating gvcfs.  Use Envel's scripts for this
