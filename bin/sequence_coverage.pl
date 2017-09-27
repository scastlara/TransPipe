#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

my $protein_name   = shift @ARGV;
my @files          = @ARGV;
my %align_coords   = ();
my $protein_length = 0;
my %align_counts   = ();

foreach my $transcriptome (@files) {
    my $tr_name = $transcriptome;
    $tr_name =~ (s/transpipe\_//g);
    open my $fh, "<", $transcriptome . "/BLAST/transcriptome-x-subject.blastx.out"
        or (print STDERR "File not found for $tr_name\n" && next);

    $align_coords{$tr_name} = ();
    while (<$fh>) {
        chomp;
        my @cols = split /\t/;
        next unless $cols[2] =~ m/$protein_name/;
        $protein_length = $cols[3] unless $protein_length;
        push @{ $align_coords{$tr_name} }, [$cols[6], $cols[7]];
    }
}

foreach my $transcriptome (keys %align_coords) {

    next if $transcriptome eq "planmine2";
    foreach my $coord (@{ $align_coords{$transcriptome} }) {
        foreach my $position ($coord->[0]..$coord->[1]) {
            $align_counts{$transcriptome}{$position} = 0 unless exists $align_counts{$transcriptome}{$position};
            $align_counts{$transcriptome}{$position}++;
        }
    }

    foreach my $coord (1..$protein_length) {
        # Add missing positions as 0 counts and print
        if (not $align_counts{$transcriptome}{$coord})  {
            $align_counts{$transcriptome}{$coord} = 0;
        }
        my $count = $align_counts{$transcriptome}{$coord};
        my $tr_name = $transcriptome;
        $tr_name =~ s/planmine/Dresden/g;
        $tr_name =~ s/cthulhu2/Cthulhu/g;
        print ucfirst($tr_name), "\t$coord\t$count\n";
    }

}
