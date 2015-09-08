#!/bin/bash

#Reyka Jayasinghe (rjayasin@genome.wustl.edu), Steven Foltz (sfoltz@genome.wustl.edu)
#4 September 2015

#QC Pipeline is a series of commands calling several scripts to extract, classify, and annotate complex events reported by Pindel-C.

#Annotated QC output is stored here: origdata/variants/ID.type.*.anno.vcf where *={germline, somatic, loh}

#The input to the pipeline consists of:
#1. sample ID
#2. cancer type
#3. absolute path to the folder of Pindel-C outputs for that sample (to the folder containing sampleID_D file(s), not the actual sampleID_D file(s))
#4. absolute path to the sample's tumor bam
#5. absolute path to the sample's normal bam
#6. project name (if needed for your application to parse sample ID, etc.)
#7. coverage minimum (default: 20, can be changed in bsub_qc.sh)
#8. steps to be completed (default: all steps, can be changed in bsub_qc.sh)

#How the script is called:
#bash /path/to/qc_pipeline.sh SAMPLE_ID CANCER_TYPE /PATH/TO/PINDEL/OUTPUT/ /PATH/TO/TUMOR/BAM.bam /PATH/TO/NORMAL/BAM.bam PROJECT_NAME

#In actually running, each individual will be sent to bsub as their own job (see bsub_qc.sh) 

###THESE VARIABLES MUST BE SET BY USER BEFORE RUNNING
gsl3="perl germline_somatic_loh_3.pl" #path to the perl script germline_somatic_loh_3.pl
p2v="pindel2vcf" #path to the executable pindel2vcf
ref="reference_sequence_used_by_pindel.fa" #path to the reference fasta used by Pindel-C
rc="read_count2.pl" #path to the executable perl script read_count2.pl
vep_dir="/your/directory/vep80" #path to your VEP directory (ex: /your/directory/vep80)
vep="$vep_dir/vep/variant_effect_predictor.pl" #path to your version of variant_effect_predictor.pl
data_dir="$vep_dir/.vep" #path to your hidden VEP directory
fasta="$data_dir/homo_sapiens/80_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa" #path to your VEP fasta
#Also, user must edit germline_somatic_loh_3.pl to set rules for identifying tumor and normal
###

if [ $# -ne 8 ]; then
        echo "Wrong number of arguments supplied. Usage: bash /path/to/pipeline.sh sample_id cancer_type /path/to/pindel/output /path/to/tumor/bam.bam path/to/normal/bam.bam project_name coverage_minimum steps"
        exit 1
fi

#Inputs
id=$1 #the sample id
type=$2 #their cancer type
echo "###" >> log/$id.$type.log
echo "New command at "$(date)": $1 $2 $3 $4 $5 $6 $7 $8" >> log/$id.$type.log
echo "###" >> log/$id.$type.log
echo "Timestamp 0: "$(date)" Sorting inputs..." >> log/$id.$type.log
pindel_path=$3 #absolute path to their pindel outputs
if [ -d $pindel_path ]; then
    echo "#Pindel output directory exists: "$pindel_path >> log/$id.$type.log
else
    echo "#Fail: Pindel output directory does not exist: "$pindel_path >> log/$id.$type.log
    exit 1
fi
tumor_bam=$4 #absolute path to their tumor bam
if [ -f $tumor_bam ]; then
    if [ -s $tumor_bam ]; then
	echo "#Tumor bam file exists and is not empty: "$tumor_bam >> log/$id.$type.log
    else
	echo "#Tumor bam file exists but is EMPTY: "$tumor_bam >> log/$id.$type.log
    fi
else
    echo "#Fail: Tumor bam file does not exist: "$tumor_bam >> log/$id.$type.log
    exit 1
fi
normal_bam=$5 #absolute path to their normal bam
if [ -f $normal_bam ]; then
    if [ -s $normal_bam ]; then
	echo "#Normal bam file exists and is not empty: "$normal_bam >> log/$id.$type.log
    else
	echo "#Normal bam file exists but is EMPTY: "$normal_bam >> log/$id.$type.log
    fi
else
    echo "#Fail: Normal bam file does not exist: "$normal_bam >> log/$id.$type.log
    exit 1
fi
project=$(echo $6 | tr '[a-z]' '[A-Z]') #project name converted to upppercase
#if your pipeline requires a specific project name:
#if [ $project == "X" ] || [ $project == "Y" ]; then
#    echo "#Project is "$project >> log/$id.$type.log
#else
#    echo "#Fail: Project must be X or Y: "$project >> log/$id.$type.log
#    exit 1
#fi
coverage_min=$7
steps=$8

#Extract complex indels
if [[ $steps == *"1"* ]]; then
    echo "Timestamp 1: "$(date)" Extracting complex indels..." >> log/$id.$type.log
    mkdir -p origdata
    #picks up complex events (insertion length is not zero, deletion and insertion lengths are not equal)
    grep ChrID $pindel_path/*_D | awk '{if($5) print}' | awk '{if ($3!=$5) print}' > origdata/$id.$type.complex
    if [ -f origdata/$id.$type.complex ]; then
	if [ -s origdata/$id.$type.complex ]; then
	    echo "#Complex indel file exists and is not empty: "origdata/$id.$type.complex >> log/$id.$type.log
	else
	    echo "#Complex indel file exists but is EMPTY: "origdata/$id.$type.complex >> log/$id.$type.log
	fi
    else
	echo "#Fail: Complex indel file does not exist: "origdata/$id.$type.complex >> log/$id.$type.log
	exit 1
    fi
fi

#Separates germline, somatic, and LOH events
#Germline is in both tumor and normal; Somatic is only in tumor; LOH is only in normal
WD="origdata/$id.$type.complex"
GERMLINE="origdata/$id.$type.germline"
SOMATIC="origdata/$id.$type.somatic"
LOH="origdata/$id.$type.loh"
if [[ $steps == *"2"* ]]; then
    echo "Timestamp 2: "$(date)" Separating germline, somatic, and LOH..." >> log/$id.$type.log
    $gsl3 $WD $project $GERMLINE $SOMATIC $LOH
    test1=$(grep "Successfully completed germline." $GERMLINE)
    test2=$(grep "Successfully completed somatic." $SOMATIC)
    test3=$(grep "Successfully completed loh." $LOH)
    if [ "$test1" == "Successfully completed germline." ]; then
	echo "#Germline successfully completed: "$GERMLINE >> log/$id.$type.log
    else
	echo "#Fail: Germline failed: "$GERMLINE >> log/$id.$type.log
	exit 1
    fi
    if [ "$test2" == "Successfully completed somatic." ]; then
	echo "#Somatic successfully completed: "$SOMATIC >> log/$id.$type.log
    else
	echo "#Fail: Somatic failed: "$SOMATIC >> log/$id.$type.log
	exit 1
    fi
    if [ "$test3" == "Successfully completed loh." ]; then
	echo "#LOH successfully completed: "$LOH >> log/$id.$type.log
    else
	echo "#Fail: LOH failed: "$LOH >> log/$id.$type.log
	exit 1
    fi
fi

#Make sure there is sufficient coverage, filter out events with low coverage
GERMLINEc="origdata/$id.$type.germline_coverage"
SOMATICc="origdata/$id.$type.somatic_coverage"
LOHc="origdata/$id.$type.loh_coverage"
if [[ $steps == *"3"* ]]; then
    echo "Timestamp 3: "$(date)" Filtering out variants with coverage less than $coverage_min reads..." >> log/$id.$type.log
    grep -v "Successfully completed germline." $GERMLINE | awk -F' ' '{a=$33+$35+$37;b=$40+$42+$44;if(a>='$coverage_min'&&b>='$coverage_min'){print}}' > $GERMLINEc
    grep -v "Successfully completed somatic." $SOMATIC | awk -F' ' '{a=$33+$35+$37;b=$40+$42+$44;if(a>='$coverage_min'&&b>='$coverage_min'){print}}' > $SOMATICc
    grep -v "Successfully completed loh." $LOH | awk -F' ' '{a=$33+$35+$37;b=$40+$42+$44;if(a>='$coverage_min'&&b>='$coverage_min'){print}}' > $LOHc
    if [ -f $GERMLINEc ]; then
	if [ -s $GERMLINEc ]; then
	    echo "#Germline file exists and is not empty: "$GERMLINEc >> log/$id.$type.log
	else
	    echo "#Germline file exists but is EMPTY: "$GERMLINEc >> log/$id.$type.log
	fi
    else
	echo "#Fail: Germline file does not exist: "$GERMLINEc >> log/$id.$type.log
	exit 1
    fi
    if [ -f $SOMATICc ]; then
	if [ -s $SOMATICc ]; then
	    echo "#Somatic file exists and is not empty: "$SOMATICc >> log/$id.$type.log
	else
	    echo "#Somatic file exists but is EMPTY: "$SOMATICc >> log/$id.$type.log
	fi
    else
	echo "#Fail: Somatic file does not exist: "$SOMATICc >> log/$id.$type.log
	exit 1
    fi
    if [ -f $LOHc ]; then
	if [ -s $LOHc ]; then
	    echo "#LOH file exists and is not empty: "$LOHc >> log/$id.$type.log
	else
	    echo "#LOH file exists but is EMPTY: "$LOHc >> log/$id.$type.log
	fi
    else
	echo "#Fail: LOH file does not exist: "$LOHc >> log/$id.$type.log
	exit 1
    fi
fi

#Makes VCF files for germline, somatic, and LOH pindel complex indel outputs
mkdir -p VCFs
if [[ $steps == *"4"* ]]; then
    for gsl in germline somatic loh; do
	echo "Timestamp 4: "$(date)" Making VCF file for $gsl..." >> log/$id.$type.log
	$p2v -p origdata/$id.$type.${gsl}_coverage -r $ref -R $ref -d 2015 -v VCFs/$id.$type.$gsl.vcf
	if [ -f VCFs/$id.$type.$gsl.vcf ]; then
	    if [ -s VCFs/$id.$type.$gsl.vcf ]; then
		echo "#VCF "$gsl" file exists and is not empty: "VCFs/$id.$type.$gsl.vcf >> log/$id.$type.log
	    else
		echo "#VCF "$gsl" file exists but is EMPTY: "VCFs/$id.$type.$gsl.vcf >> log/$id.$type.log
	    fi
	else
	    echo "#Fail: VCF "$gsl" does not exist: "VCFs/$id.$type.$gsl.vcf >> log/$id.$type.log
	    exit 1
	fi
    done
fi

#Runs read counts perl script for each germline, somatic, and LOH 
mkdir -p read_counts
#for gsl in germline somatic loh; do
for gsl in somatic loh; do
    vcf="VCFs/$id.$type.$gsl.vcf"
    outtumor="read_counts/$id.$type.$gsl.tumor.rc"
    outnormal="read_counts/$id.$type.$gsl.normal.rc"
    if [[ $steps == *"5"* ]]; then
	echo "Timestamp 5: "$(date)" Running read counts for tumor $gsl..." >> log/$id.$type.log
	perl $rc -id $id $ref $vcf $tumor_bam 20 $outtumor
	if [ -f $outtumor ]; then
	    if [ -s $outtumor ]; then
		echo "#Read counts "$gsl" tumor file exists and is not empty: "$outtumor >> log/$id.$type.log
	    else
		echo "#Read counts "$gsl" tumor file exists but is EMPTY: "$outtumor >> log/$id.$type.log
	    fi
	else
	    echo "#Fail: Read counts "$gsl" tumor does not exist: "$outtumor >> log/$id.$type.log
	    exit 1
	fi
    fi
    if [[ $steps == *"6"* ]]; then
	echo "Timestamp 6: "$(date)" Running read counts for normal $gsl..." >> log/$id.$type.log
	perl $rc -id $id $ref $vcf $normal_bam 20 $outnormal #20 refers to mapping quality minimum
	if [ -f $outtumor ]; then
	    if [ -s $outtumor ]; then
		echo "#Read counts "$gsl" normal file exists and is not empty: "$outnormal >> log/$id.$type.log
	    else
		echo "#Read counts "$gsl" normal file exists but is EMPTY: "$outnormal >> log/$id.$type.log
	    fi
	else
            echo "#Fail: Read counts "$gsl" normal does not exist: "$outnormal >> log/$id.$type.log
            exit 1
	fi
    fi
done

#Reclassify germline, somatic, and LOH based on read counts data
if [[ $steps == *"7"* ]]; then
    echo "Timestamp 7: "$(date)" Reclassifying complex indels based on read counts..." >> log/$id.$type.log
    #Support for "somatic" event in normal bam
    awk '{if ($16>0) print "\""substr($6,2)"\"\tChrID "$2"\tBP "$3}' read_counts/$id.$type.somatic.normal.rc | sort | uniq > $id.TEMP
    #No support for "somatic" event in tumor bam
    awk '{if ($16==0) print "\""substr($6,2)"\"\tChrID "$2"\tBP "$3}' read_counts/$id.$type.somatic.tumor.rc | sort | uniq > $id.TEMP_NO_SOMATIC
    #Support for "LOH" event in tumor bam 
    awk '{if ($16>0) print "\""substr($6,2)"\"\tChrID "$2"\tBP "$3}' read_counts/$id.$type.loh.tumor.rc | sort | uniq > $id.TEMP_2
    #No support for "LOH" event in normal bam
    awk '{if ($16==0) print "\""substr($6,2)"\"\tChrID "$2"\tBP "$3}' read_counts/$id.$type.loh.normal.rc | sort | uniq > $id.TEMP_2_NO_LOH
    #Grep misclassified somatic events from original "somatic" file
    grep -wf $id.TEMP $SOMATICc | sort > $id.TEMP_somatic
    #Grep misclassified LOH events from original "LOH" file
    grep -wf $id.TEMP_2 $LOHc | sort > $id.TEMP_loh
    #Sort "somatic" and "LOH" events for use in comm
    sort $SOMATICc > $id.sorted_og_somatic
    sort $LOHc > $id.sorted_og_loh
    #Save reclassified and filtered germline, somatic, and loh files here
    GERMLINEf="origdata/$id.$type.germline_filtered"
    SOMATICf="origdata/$id.$type.somatic_filtered"
    LOHf="origdata/$id.$type.loh_filtered"
    #Cat events originally classified or reclassified as germline
    cat $GERMLINEc $id.TEMP_somatic $id.TEMP_loh | sort | uniq > $GERMLINEf
    #Print lines from "somatic" file not reclassified as germline (comm -13 prints lines unique to second file)
    #Then remove "somatic" events without support in tumor bam
    comm -13 $id.TEMP_somatic $id.sorted_og_somatic | sort | uniq | grep -vf $id.TEMP_NO_SOMATIC > $SOMATICf
    #Print lines from "LOH" file not reclassified as germline (comm -13 prints lines unique to second file)
    #Then remove "LOH" events without support in normal bam
    comm -13 $id.TEMP_loh $id.sorted_og_loh | sort | uniq | grep -vf $id.TEMP_2_NO_LOH > $LOHf
    #Remove temporary files
    rm -f $id.TEMP $id.TEMP_NO_SOMATIC $id.TEMP_2 $id.TEMP_2_NO_LOH $id.TEMP_somatic $id.TEMP_loh $id.sorted_og_somatic $id.sorted_og_loh
    if [ -f $GERMLINEf ]; then
	if [ -s $GERMLINEf ]; then
	    echo "#Germline file exists and is not empty: "$GERMLINEf >> log/$id.$type.log
	else
	    echo "#Germline file exists but is EMPTY: "$GERMLINEf >> log/$id.$type.log
	fi
    else
	echo "#Fail: Germline file does not exist: "$GERMLINEf >> log/$id.$type.log
	exit 1
    fi
    if [ -f $SOMATICf ]; then
	if [ -s $SOMATICf ]; then
	    echo "#Somatic file exists and is not empty: "$SOMATICf >> log/$id.$type.log
	else
	    echo "#Somatic file exists but is EMPTY: "$SOMATICf >> log/$id.$type.log
	fi
    else
	echo "#Fail: Somatic file does not exist: "$SOMATICf >> log/$id.$type.log
	exit 1
    fi
    if [ -f $LOHf ]; then
	if [ -s $LOHf ]; then
	    echo "#LOH file exists and is not empty: "$LOHf >> log/$id.$type.log
	else
            echo "#LOH file exists but is EMPTY: "$LOHf >> log/$id.$type.log
	fi
    else
	echo "#Fail: LOH file does not exist: "$LOHf >> log/$id.$type.log
	exit 1
    fi
fi

#Final sections: Make filtered VCFs, run VEP to annotate variants
anno_input="origdata/annotate"
anno_output="origdata/variants"
mkdir -p $anno_input $anno_output 
for gsl in germline somatic loh; do
    #Makes VCF files for germline, somatic, and LOH pindel complex indel outputs
    if [[ $steps == *"8"* ]]; then
	echo "Timestamp 8: "$(date)" Making VCF file for VEP input for $gsl..." >> log/$id.$type.log
	$p2v -p origdata/$id.$type.${gsl}_filtered -r $ref -R $ref -d 2015 -v VCFs/$id.$type.$gsl.filtered.vcf
	if [ -f VCFs/$id.$type.$gsl.filtered.vcf ]; then
	    if [ -s VCFs/$id.$type.$gsl.filtered.vcf ]; then
		echo "#VCF filtered "$gsl" file exists and is not empty: "VCFs/$id.$type.$gsl.filtered.vcf >> log/$id.$type.log
	    else
		echo "#VCF filtered "$gsl" file exists but is EMPTY: "VCFs/$id.$type.$gsl.filtered.vcf >> log/$id.$type.log
	    fi
	else
            echo "#Fail: VCF filtered "$gsl" does not exist: "VCFs/$id.$type.$gsl.filtered.vcf >> log/$id.$type.log
            exit 1
	fi
    fi
    #Runs VEP to annotate the variants
    if [[ $steps == *"9"* ]]; then
	echo "Timestamp 9: "$(date)" Annotating for $gsl..."
	perl $vep --everything -i VCFs/$id.$type.$gsl.filtered.vcf --format vcf --vcf -out $anno_output/$id.$type.$gsl.anno.vcf --dir $data_dir --assembly GRCh37 --cache --offline --fork 4
	if [ -f $anno_output/$id.$type.$gsl.anno.vcf ]; then
	    if [ -s $anno_output/$id.$type.$gsl.anno.vcf ]; then
		echo "#Annotation output file exists and is not empty: "$anno_output/$id.$type.$gsl.anno.vcf >> log/$id.$type.log
	    else
		echo "#Annotation output file exists but is EMPTY: "$anno_output/$id.$type.$gsl.anno.vcf >> log/$id.$type.log
	    fi
	else
            echo "#Fail: Annotation output file does not exist: "$anno_output/$id.$type.$gsl.anno.vcf >> log/$id.$type.log
            exit 1
	fi
    fi
done

echo "End of pipeline! Ran steps: "$steps". Finished "$(date) >> log/$id.$type.log
