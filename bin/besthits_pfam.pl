#!/usr/bin/perl

use warnings;
use strict;

my $file = shift @ARGV;

open my $fh, "<", $file
    or die "Can't open $file : $!\n";
my %best_hits = ();

while (<$fh>) {
    chomp;
    my @columns = split /\t/;
    if (not exists $best_hits{$columns[0]}->{$columns[1]}) {
        $best_hits{$columns[0]}->{$columns[1]} = [$columns[2], $columns[4]];
    } else {
        if ($columns[4] > $best_hits{$columns[0]}->{$columns[1]}[1]) {
            $best_hits{$columns[0]}->{$columns[1]} = [$columns[2], $columns[4]];
        }
    }
}

foreach my $code (keys %best_hits) {
    foreach my $hit (keys %{ $best_hits{$code} }) {
        print "$code\t$hit\t$best_hits{$code}->{$hit}->[0]\t$best_hits{$code}->{$hit}->[1]\n";
    }
}
