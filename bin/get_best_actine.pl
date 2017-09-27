#!/usr/bin/perl
use warnings;
use strict;

my %best_pairs = ();
$best_pairs{'ACTB'} = ();
$best_pairs{'ACTG1'} = ();
$best_pairs{'ACTB'}{'eval'}  = 100;
$best_pairs{'ACTG1'}{'eval'} = 100;

while (<>) {
    chomp;
    my @cols  = split /\t/;
    if ($cols[2] eq "ACTB") {
        my $eval = 10;
        if ($cols[11] =~ m/(\-\d+)/) {
            $eval = $1;
        } else {
            next;
        }
        if ($eval < $best_pairs{'ACTB'}{'eval'}) {
            $best_pairs{'ACTB'}{'eval'} = $eval;
            $best_pairs{'ACTB'}{'homolog'} = $cols[0];
        }
    }
    if ($cols[3] eq "ACTB") {
        my $eval = 10;
        if ($cols[12] =~ m/(\-\d+)/) {
            $eval = $1;
        } else {
            next;
        }
        if ($eval < $best_pairs{'ACTB'}{'eval'}) {
            $best_pairs{'ACTB'}{'eval'}    = $eval;
            $best_pairs{'ACTB'}{'homolog'} = $cols[1];
        }
    }

    if ($cols[2] eq "ACTG1") {
        my $eval = 10;
        if ($cols[11] =~ m/(\-\d+)/) {
            $eval = $1;
        } else {
            next;
        }
        if ($eval < $best_pairs{'ACTG1'}{'eval'}) {
            $best_pairs{'ACTG1'}{'eval'} = $eval;
            $best_pairs{'ACTG1'}{'homolog'} = $cols[0];
        }
    }
    if ($cols[3] eq "ACTG1") {
        my $eval = 10;
        if ($cols[12] =~ m/(\-\d+)/) {
            $eval = $1;
        } else {
            next;
        }
        if ($eval < $best_pairs{'ACTG1'}{'eval'}) {
            $best_pairs{'ACTG1'}{'eval'}    = $eval;
            $best_pairs{'ACTG1'}{'homolog'} = $cols[1];
        }
    }
}

foreach my $gene (qw(ACTB ACGT1)){
    print "$gene\t$best_pairs{$gene}{homolog}\t$best_pairs{$gene}{eval}\n"
}
