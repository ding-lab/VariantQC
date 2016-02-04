#Reyka Jayasinghe, Steven Foltz
#August 2015

#Initial script classifies complex events into germline, somatic, or loh based on having read support in the tumor and/or normal bam.
#Germline has support in both tumor and normal.
#Somatic has support in tumor only.
#LOH has support in normal only. 

use strict;
my $usage =<<USAGE;
 Usage: $0 <complex> <project> <germline> <somatic> <loh>
    Where <complex> is the complex indel file for this sample
    Where <project> is a project identifier, useful for parsing sample IDs into tumor or normal
    Where <germline> <somatic> and <loh> are the output files for germline, somatic, and LOH mutations.
    USER MUST DEFINE RULES FOR PARSING SAMPLE IDS INTO TUMOR OR NORMAL (see line 37-38)
USAGE

my $file=$ARGV[0];
my $project=$ARGV[1];
my $GERMLINE=$ARGV[2];
my $SOMATIC=$ARGV[3];
my $LOH=$ARGV[4];

open(COMPLEX,'<',$file) or die "Couldn't open complex $file.";
open(SOMATIC,'>',$SOMATIC) or die "Couldn't open file for writing $SOMATIC.";
open(LOH,'>',$LOH) or die "Couldn't open file for writing $LOH.";
open(GERMLINE,'>',$GERMLINE) or die "Couldn't open file for writing $GERMLINE.";
while(my $line=<COMPLEX>){
	chomp ($line);
	my @pindel=split(/\s/,$line);
	my $size=@pindel; 
	my $supsamples=$pindel[29];
	#Identify samples and store their sample type (tumor/primary/normal/germline)
	my @sample1=$pindel[31];
	my @sample2=$pindel[38];
	#USER MUST WRITE OWN RULES FOR DETERMINING THE TYPE OF EACH BAM (tumor or normal)
	my $type1=""; #must be 'tumor' or 'normal', opposite of type2
	my $type2=""; #must be 'tumor' or 'normal', opposite of type1
    #CHECK TO SEE IF NO ENTERIES WERE PROVIDED FOR BAM INFORMATION - REFER TO ABOVE TWO LINES 
    if ($type1=~/^$/){
        print STDERR "ERROR: User needs to define rules for parsing sample IDS into tumor and normal!! See germline_somatic_loh_3.pl Lines 37 and 38\n";
        die;
    }
    if ($type2=~/^$/){
        print STDERR "ERROR: User needs to define rules for parsing sample IDS into tumor and normal!! See germline_somatic_loh_3.pl Lines 37 and 38\n";
        die;
    }	
	#Set column 34 as the anchor column for later comparisons
	#Script works only if two samples are provided in pindel output file (size=45)
	my $column=34 if $size==45;

	#If the number of supporting samples is 1, run check_sup1()
	if (($supsamples==1)){
		check_sup1($column,$type1,$type2,\@pindel,$line);
	}
	#If the number of supporting samples i 2, run check_sup2()
	if (($supsamples==2)){
		check_sup2($column,$type1,$type2,\@pindel,$line);
	}
}

sub check_sup1{
    #Get passed arguments
    my ($column,$type1,$type2,$pindel1,$line)=@_;
    my @pindel=@{$pindel1};
    #Check first sample
    if ($pindel[$column]>0||$pindel[$column+1]>0||$pindel[$column+2]>0||$pindel[$column+3]>0){	
	if ($type1 eq "tumor" && $type2 eq "normal"){
	    print SOMATIC "$line\n";
	}
	if ($type1 eq "normal" && $type2 eq "tumor"){
	    print LOH "$line\n";
	}
    }
    #Check second sample
    elsif ($pindel[$column+7]>0||$pindel[$column+8]>0||$pindel[$column+9]>0||$pindel[$column+10]>0){
	if ($type1 eq "tumor" && $type2 eq "normal"){
	    print LOH "$line\n";
	}
	if ($type1 eq "normal" && $type2 eq "tumor"){
		print SOMATIC "$line\n";
	}
    }
}

sub check_sup2{
    #Get passed arguments
    my ($column,$type1,$type2,$pindel1,$line)=@_;
    my @pindel=@{$pindel1};
    #Check both samples
    if (($pindel[$column]>0||$pindel[$column+1]>0||$pindel[$column+2]>0||$pindel[$column+3]>0)&&($pindel[$column+7]>0||$pindel[$column+8]>0||$pindel[$column+9]>0||$pindel[$column+10]>0)){
	if ($type1 eq "tumor" && $type2 eq "normal"){
	    print GERMLINE "$line\n";
	}
	if ($type1 eq "normal" && $type2 eq "tumor"){
	    print GERMLINE "$line\n";
	}
    }
}

#QC steps
print GERMLINE "Successfully completed germline.\n";
print SOMATIC "Successfully completed somatic.\n";
print LOH "Successfully completed loh.\n";

close GERMLINE;
close SOMATIC;
close LOH;
