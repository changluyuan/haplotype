#!/usr/bin/perl
######Usage: CGMAP_generator.minus.pl target.bed

system "bedtools makewindows -b $ARGV[0] -w 600 > $ARGV[0].split.600.bed";
system "bedtools getfasta -fi /gpfs/genomedb/b37/v37_decoy_plus_phage.fasta -bed $ARGV[0].split.600.bed -fo $ARGV[0].temp.get.fasta";
open (FASTA, "<$ARGV[0].temp.get.fasta");
open (CGMAP, ">$ARGV[0].split.minus.300.CGMAP");
while(<FASTA>){
    chomp;
    if($_=~">"){
        $current_window = $_;
        $current_window =~ s/>//g;
    }
    else{
        $sequence = $_;
        ($chr, $range) = split ":", $current_window;
        ($start, $end) = split "-", $range;
        @array = split //, $sequence;
        $n=0;
        $CGcounts = -1;
        while($n<length($sequence)){
                if($array[$n] eq "G" || $array[$n] eq "g"){
                    $CGcounts ++;
                    $absolute_position = $start + $n - 1;
                    print CGMAP "$chr\t$start\t$end\t$absolute_position\t$CGcounts\n";
                }
                $n++;
        }
    }
}
close FASTA;
close CGMAP;
system "rm -f $ARGV[1].temp.get.fasta";
