#!/usr/bin/perl
use warnings;
use strict;

my $fw_file   = shift @ARGV;
my $rec_file  = shift @ARGV;
my %reciprocs = ();

open my $FW, "<", "$fw_file"
    or die "Can't open $fw_file : $!\n";

while (<$FW>) {
    chomp;
    my @columns = split /\t/, $_;
    last if $columns[0] eq "T";
    my $query = $columns[1];
    my ($target) = split /\s/, $columns[13];
    $target =~ s/\|.+//;
    $reciprocs{$query}->{$target} = $columns[12];
}

close $FW;

open my $RV, "<", "$rec_file"
    or die "Can't open $rec_file : $!\n";

print "TRANSCRIPT\tSUBJECT\tTRANS_COV\tSUBJ_COV\n";
while (<$RV>) {
    chomp;
    my @columns = split /\t/, $_;
    last if $columns[0] eq "T";
    my $query = $columns[1];
    my ($target) = split /\s/, $columns[13];
    $query =~ s/\|.+//;
    print "$target\t$query\t$reciprocs{$target}->{$query}\t$columns[12]\n" if exists $reciprocs{$target}->{$query};
}
close $RV;
