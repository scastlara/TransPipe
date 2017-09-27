#!/usr/bin/perl

use warnings;
use strict;

my $file = shift @ARGV;

open my $fh, "<", $file
    or die "Can't open $file : $!\n";
my %brhs = ();
my $first_lett = "";
print "TRANSCRIPT\tSUBJECT\tSCORE\n";

while(<$fh>) {
    chomp;
    my ($letter, $q, $t, $sc) = split /\t/;
    $first_lett = $letter if not $first_lett;
    if ($letter eq $first_lett) {
        $brhs{$q}->{$t} = $sc;
    } else {
        if (exists $brhs{$t}->{$q}) {
            if ($q =~ m/HGNC/) {
                print "$t\t$q\t$sc\n";
            } else {
                print "$q\t$t\t$sc\n";
            }
        }
    }

}
