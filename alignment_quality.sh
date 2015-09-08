#SMF 7 September 2015
#Script does initial QC on bams, identifying which bams have <80% perfectly aligned reads. 

#rm -f bad_bam.txt good_bam.txt
while read sample_type_bam; do
    sample_name=$(echo $sample_type_bam | cut -f1 -d' ')
    type=$(echo $sample_type_bam | cut -f2 -d' ')
    bam=$(echo $sample_type_bam | cut -f3 -d' ')
    if [ -f $bam ]; then 
	samtools view $bam 20 | cut -f6 > $sample_name.temp #view bam file for all of chromosome 20
	total=$(wc -l $sample_name.temp | cut -f1 -d' ') #total number of reads
	good=$(grep -c 100M $sample_name.temp) #number of perfectly matching 100M reads
	prop=$(echo 100*$good/$total | bc -l) #proportion of perfectly matching 100M reads
	if [ $(echo ${prop%%.*}) -lt 80 ]; then #grabs the digits to left of decimal point
	    echo $sample_name $type $bam $prop | tr ' ' '\t' >> bad_bam.txt #if prop < 80, bad
	else
	    echo $sample_name $type $bam $prop | tr ' ' '\t' >> good_bam.txt #if prop >= 80, good
	fi
	rm -f $sample_name.temp
    else
	echo "$sample_name $type $bam (no bam)" >> bad_bam.txt
    fi
done < bams
#input looks like
#SAMPLE_NAME CANCER_TYPE /path/to/bamfile.bam