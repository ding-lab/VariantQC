# VariantQC
Variant quality checking scripts.
# How to run QC
Main QC script is run using bsub_qc.sh, which initiates the main qc_pipeline.sh. The input to bsub_qc.sh is described in the file.

#Reyka Jayasinghe (rjayasin@genome.wustl.edu) and Steven Foltz (sfoltz@genome.wustl.edu)

#The Variant QC script will go through the following steps:
1. Extract complex insertions and deletions from pindel output. 	
2. Identify somatic, germline, and loss of heterozygosity(loh) events.
3. Filter out low coverage sites (20 read min).
4. Make unfiltered VCF for germline, somatic and loh events.
5. Run readcount tool on tumor sample. Performing readcount analysis will determine if somatic and loh events are appropriately classified (Note: Not run for germline). 
6. Run readcount tool on normal sample. Performing readcount analysis will determine if somatic and loh events are appropriately classified (Note: Not run for germline).
7. Reclassify germline, somatic, and loh based on read count data of somatic events.  
8. Making VCFs for filtered pindel output for VEP input & annotate final filtered VCF using VEP.
