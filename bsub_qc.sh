#!/gsc/bin/bash

#Steven Foltz (sfoltz@genome.wustl.edu), Reyka Jayasinghe (rjayasin@genome.wustl.edu)
#4 September 2015

#Script runs complex indel QC pipeline for individuals specfied in the input file. Each individual is run in parallel using bsub. 

#Input file has six columns: SAMPLE_ID CANCER_TYPE /PATH/TO/PINDEL/OUTPUT/ /PATH/TO/TUMOR/BAM.bam /PATH/TO/NORMAL/BAM.bam PROJECT_NAME
#Each row corresponds to a unique individual. 

#How the pipeline script is called:
#bash /path/to/qc_pipeline.sh SAMPLE_ID CANCER_TYPE /PATH/TO/PINDEL/OUTPUT/ /PATH/TO/TUMOR/BAM.bam /PATH/TO/NORMAL/BAM.bam PROJECT_NAME 

if [ $# -ne 1 ]; then
        echo "Wrong number of arguments supplied. Usage: bash pipeline_bsub.sh six_column_input_file"
        exit 1
fi

coverage_min=20 #complex events called with coverage below this threhold are removed
steps="1,2,3,4,5,6,7,8,9" #user can select which parts of the QC pipeline to run
#1 Extracting complex indels
#2 Separating germline, somatic, and LOH
#3 Filtering out low coverage
#4 Making VCF file for germline/somatic/loh
#5 Running read counts for tumor germline/somatic/loh
#6 Running read counts for normal germline/somatic/loh
#7 Reclassifying complex indels based on read counts
#8 Making VCF files for VEP input for germline/somatic/loh
#9 Annotating for germline/somatic/loh

file_containing_each_sample_info=$1
#absolute path to the actual QC pipeline
pipeline="bash qc_pipeline.sh"

#screen error and output goes in these folders
mkdir -p log

while read line; do
    id=$(echo $line | cut -f1 -d' ') #Sample ID
    type=$(echo $line | cut -f2 -d' ') #Cancer type
    pindel_output=$(echo $line | cut -f3 -d' ') #Path to pindel output folder
    tumor_bam=$(echo $line | cut -f4 -d' ') #Path to tumor bam
    normal_bam=$(echo $line | cut -f5 -d' ') #Path to normal bam
    project=$(echo $line | cut -f6 -d' ') #Project name
    #Using bsub:
    bsub -e log/$id.$type.err -o log/$id.$type.out $pipeline $id $type $pindel_output $tumor_bam $normal_bam $project $coverage_min $steps
    echo $id >> Samples_$type #Adds sample ID to list of others with same cancer type
done < $file_containing_each_sample_info
