#!/usr/bin/perl
#use Data::Dumper;
print "Usage: epiread.paired.to.window.012.automatic.pl YOUR-FAVOURITE-BAM-FILE YOUR-FAVOURITE-METHYLATION-WINDOW OUTPUT-name window-size\n";
print "parameter: 1: bamfile; 2. meth bed; 3. output name; 4.window-size \n";
print "OUTPUT: window_name, a fixed-length read-type with 0 as unknown, 1 as methylated, and 2 as unmethylated\n";
print "NOTE: the default thread is 4!\n";

if(@ARGV<=3) { die "not sufficient parameter"; }

system "samtools view -t 4 -L $ARGV[1] -b -h $ARGV[0] > $ARGV[0].target.bam";
system "samtools index $ARGV[0].target.bam";
system "/gpfs/bin/biscuit-2/biscuit epiread -r /gpfs/genomedb/b37/v37_decoy_plus_phage.fasta -i $ARGV[0].target.bam > $ARGV[0].paired.epiread";

## 3. From the CG-map and epiread, output fixed length read for downstream processing


open (EPIREAD, "<$ARGV[0].paired.epiread");
open (REFERENCEPOSITION, "<$ARGV[1].split.minus.$ARGV[3].CGMAP");
open (OUTPUT, ">$ARGV[2]");

my %h;

while(<REFERENCEPOSITION>){
        chomp;
        @array = split "\t", $_;
        $window = "$array[0]\t$array[1]\t$array[2]";
        $position = $array[3];
        $CGposition = $array[4];
        $position_to_window{"$array[0]:$position"} = $window;
        if($window_sum_CG{$window}<$CGposition){
            $window_sum_CG{$window} =$CGposition;
        }
        $position_absolute_site_in_window{"$array[0]:$position"} = $CGposition ;
}


print OUTPUT "chr\tstart\tend\tread\n";

while(<EPIREAD>){
        chomp;
        @array = split "\t", $_;
        $chr = $array[0];
        $start = $array[4];
        $end = $start + 1;
        $sequence = $array[5];
        if(!defined($position_to_window{"$chr:$start"})){next;}
        if($array[3] eq "+"){next;}
        $current_window = $position_to_window{"$chr:$start"};
        $heading_N = $position_absolute_site_in_window{"$chr:$start"};
        $sequence_length = length($sequence);
        $tailing_N = $window_sum_CG{$current_window} - $heading_N - $sequence_length;
        $text = $sequence;
        $text = ( '0,' x $heading_N ) . $text;
        $text =  $text . ( '0,' x $tailing_N );
        $text =~ s/C/1,/gi;
        $text =~ s/T/-1,/gi;
        $text =~ s/N/0,/gi;
        $text =~ s/,$//gi;
		$h{"$array[1]-$array[2]"} = $text;
        print OUTPUT "$current_window\t$text\n";
}
close EPIREAD;
close CGMAP;
close OUTPUT;

#system "perl /gpfs/output/Tequila/MCBS/bin/fixwindow.pl $ARGV[2] $ARGV[3]";

################################################################################################


open (INPUT_WINDOW, "<$ARGV[2]");
open (FIXED_WINDOW, ">$ARGV[2].fixed");

$maxnum = $ARGV[3];

$maxlinelength = 0;

if($maxlinelength<=$maxnum){
    $maxlinelength = $maxnum;
}

$count = 0;
while(<INPUT_WINDOW>){
    chomp;
    $count ++;
    $hashtable{$count} = $_;
}


foreach $count (sort keys %hashtable){
    $line = $hashtable{$count};
    @array = split ",", $line;
    while(@array < $maxlinelength){
        push(@array, 0);
    }
    while(@array > $maxlinelength){
        $voidplaceholder = pop(@array);
    }
    $printline = join ",", @array;
    print FIXED_WINDOW "$printline\n";
}



################################################################################################

open (SORTED, ">$ARGV[2].sorted");

print SORTED "chr\tstart\tend\tread\n";

close SORTED;

system "sort -k1,1V -k2,2n -k3,3n $ARGV[2].fixed | grep -v chr >> $ARGV[2].sorted;";

system "mv $ARGV[2].sorted $ARGV[2].readbed";

system "/gpfs/bin/htslib/bgzip -f $ARGV[2].readbed && /gpfs/bin/htslib/tabix -f -S1 -b2 -e3 -0 $ARGV[2].readbed.gz;";
system "rm -f $ARGV[2]  $ARGV[2].fixed";

#system "rm -f $ARGV[1].temp.CGMAP";
#system "rm -f $ARGV[0].temp.$randomnum.bam*";
#system "rm -f $ARGV[1].$randomnum.*";



#====annotate bam,add by yh===
#open F,"sambamba view -h $ARGV[0].target.bam|";
#open O,"|/gpfs/bin/samtools-1.3.1/samtools view -Sbh -F 1024 > $ARGV[0].minus.tag.bam";
#while(<F>){
#	if(/^\@/){
#		print O $_;
#	}else{
#		chomp;
#		my @a = split /\t/;
 #       my $k;
  #      if($a[1] & 0x40){
            #R1
   #         $k = "$a[0]-1";
    #    }elsif($a[1] & 0x80){
            #R2
     #       $k = "$a[0]-2";
      #  }
       # print O my $content = exists $h{$k} ? "$_\tBC:Z:$h{$k}\n" : "$_\n";
	#}
#}
#close O;
#`/gpfs/bin/samtools-1.3.1/samtools index $ARGV[0].minus.tag.bam`;
