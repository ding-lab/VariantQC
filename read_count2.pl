#! /usr/bin/perl -w
use strict;
use Getopt::Long;



my($help,$bam_files,$bam_batch,$id);
my %hash;
&GetOptions('bam_file=s'=>\$bam_files,'bam_bat=s'=>\$bam_batch,'help|h'=>\$help,'id=s'=>\$id);


if($#ARGV<4 or defined $help)
{
  print "USAGE: $0 [option] fasta_file vcf_file bam_file mapping_quality_cutoff output_file.\n";
  print "       -id for multiple process.\n";
  exit;
}

my $ref_file=shift;
my $vcf_file=shift;
my $bam_file=shift;
my $MQ_cutoff=shift;
my $out_file=shift;
if( -e $out_file)
{
  unlink  $out_file;
}

my $key;
my %str;
open SEQ,"<$ref_file" or die "can not open ref sequence file:$!";
while(<SEQ>)
{
  chomp;
  if(/>/)
  {
    s/>//;
    s/\s+.*//;
    $key=$_;
  }
  else
  {
    $str{$key}.=uc $_;
  }
}
close SEQ;

my $ref;
my $mod_ref;
my $num=0;
open VCF,"<$vcf_file" or die "can not open vcf file:$!";
while(<VCF>)
{
  next if /^#/;
  if(/^\s*$/){next;}
  chomp;
  my @arr=split(/\s+/);
  $num++;
  my $tem1;
  if($arr[1]>=10001){$tem1=substr($str{$arr[0]},$arr[1]-10001,10000);}
  else{$tem1=substr($str{$arr[0]},0,$arr[1]);}
  my $tem2=substr($str{$arr[0]},$arr[1]+length($arr[3])-1,10000);
#  print substr($str{$arr[0]},$arr[1]-2001,4000+length($arr[3])),"\n";
  $ref=$tem1.$arr[3].$tem2;
  my $key;
  if(defined $id)
  {
    $key=$id.'_'.$arr[0].'_'.$num;
  }
  else{$key=$arr[0].'_'.$num;}
  my $fafile=$key.'.fa';
  if( -e $fafile)
  {
    unlink  $fafile;
  }
  print_out($key,'ref',$ref);
  $mod_ref=$tem1.$arr[4].$tem2;
  print_out($key,'alt',$mod_ref);
  my $start=$arr[1]-2000;
  my $end=$arr[1]+2000+length($arr[3]);
  print $start, "\t", $end, "\n";
  
#  system(" samtools view $bam_file $arr[0]:$start-$end| awk '{print \"@\"\$1;print \$10;print \"+\";print \$11;}' > ${key}_reads.fastq ");
  system(" samtools view $bam_file $arr[0]:$start-$end -b > ${key}_reads.bam");
  `java -jar /gsc/scripts/pkg/bio/picard/picard-tools-1.92/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT I=${key}_reads.bam F=${key}_reads_1.fastq F2=${key}_reads_2.fastq`;
  

print "***start $key bwa **********\n";
  my $bwa_index="bwa index $key.fa";
  system ("$bwa_index");
  my $samtools_faidx="samtools faidx $key.fa";
  system ("$samtools_faidx");
  my $bwa_cmd1 = "bwa aln -t4 -q 5 $key.fa ${key}_reads_1.fastq > ${key}_1.sai";
  system ("$bwa_cmd1"); 
  my $bwa_cmd2 = "bwa aln -t4 -q 5 $key.fa ${key}_reads_2.fastq > ${key}_2.sai";
  system ("$bwa_cmd2"); 
  my $bwa_cmd = "bwa sampe -a 600 $key.fa ${key}_1.sai ${key}_2.sai ${key}_reads_1.fastq ${key}_reads_2.fastq > $key.sam"; 
  system ("$bwa_cmd");

  #my $samtools_view = "samtools view -bS $key.sam > $key.bam"; 
  #system ("$samtools_view");
  #my $samtools_sort = "samtools sort $key.bam $key.sort"; 
  #system ("$samtools_sort");
  #my $samtools_index = "samtools index $key.sort.bam"; 
  #system ("$samtools_index");

  my $samfile=$key.'.sam';
  my $head_len=length($tem1);
  &analysis($samfile,$head_len,\@arr, $MQ_cutoff, $out_file, $vcf_file);
  
  my $rm_bam="rm -rf ${key}_reads.bam";
  system ("$rm_bam");
  my $rm_fa="rm -rf $key.fa*";
  system ("$rm_fa");
  my $rm_sai="rm -rf ${key}_*.sai";
  system ("$rm_sai");
  my $rm_sam="rm -rf $key.sam";
  system ("$rm_sam");
  my $rm_fastq="rm -rf ${key}_reads_*.fastq";
  system ("$rm_fastq");
}

sub analysis
{
  my $file=shift;
  my $head_len=shift;
  my $arrref=shift;
  my $MQ=shift;
  my $out_file=shift;
  my $vcf_name=shift;
  open SAM,"<$file" or die "can not open sam file:$!";
  open STAT,">>$out_file" or die "can not open stats file:$!";
  my $ref_num=0;
  my $alt_num=0;
  while(<SAM>)
  {
    my @arr=split;
    if($arr[2] eq 'alt')
    {
      if($arr[3] + length($arr[9]) <= $head_len + length(${$arrref}[4]) or $arr[3]>$head_len+1){next;}
      if(substr($arr[9],$head_len-$arr[3]+1,length(${$arrref}[4])) eq ${$arrref}[4] and $arr[4] >= $MQ) {$alt_num++;}
    }
    if($arr[2] eq 'ref')
    {
      if($arr[3]+length($arr[9]) <=$head_len+length(${$arrref}[3]) or $arr[3]>$head_len+1){next;}
      if(substr($arr[9],$head_len-$arr[3]+1,length(${$arrref}[3])) eq ${$arrref}[3] and $arr[4] >= $MQ) {$ref_num++;}   
    }
  } 

  close SAM;
  $"="\t";
  print STAT $vcf_name, "\t@{$arrref}\t";
  print STAT "\tref\t",$ref_num,"\t";
  print STAT "\talt\t",$alt_num,"\n";
  close STAT;
}

sub print_out
{
  my $k=shift;
  my $head=shift;
  my $seq=shift;
  open REFOUT,">>$k.fa" or die "can not open $k.fa file:$!"; 
  print REFOUT ">",$head,"\n";
  for(my $n=0;$n<=length($seq);$n+=60)
  {
    print REFOUT substr($seq,$n,60),"\n";
  }
  close REFOUT;
}      

exit;
