#!/usr/bin/perl
use strict;
use warnings;

my $hmmer_file     = shift @ARGV;
my %seq_2_pfam = ();

open my $fh, "<", $hmmer_file
    or die "Can't open $hmmer_file : $!\n";

while (<$fh>) {
    chomp;
    next if /^#/;
    my @fields = split /\s+/, $_;
    my ($seq,
    $pfam,
    $seq_start,
    $seq_end,
    $pfam_start,
    $pfam_end) = ($fields[0], $fields[4], $fields[17], $fields[18], $fields[15], $fields[16]);
    $seq_2_pfam{$seq} = [] unless exists $seq_2_pfam{$seq};
    push @{ $seq_2_pfam{$seq}}, [$pfam, $seq_start, $seq_end, $pfam_start, $pfam_end];
}


foreach my $seq (keys %seq_2_pfam) {
    my @pfams    = @{ $seq_2_pfam{$seq} };
    print "$seq\t";
    my @sorted = sort {
        $a->[1] <=> $b->[1]
    } @pfams;

    foreach my $pf (@sorted) {
        print $pf->[0], " ", "$pf->[1]-$pf->[2]", " $pf->[3]-$pf->[4]", "\t";
    }
    print "\n";
}
