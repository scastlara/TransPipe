#!/usr/bin/perl
use strict;
use warnings;
use JSON;
use Data::Dumper;


# MAIN
my $file = shift @ARGV;
if (not $file) {
    die "Gimme a file\n";
}

my $exp_data = read_json($file);

for my $transcript (@{ $exp_data->{rows} }) {
    print "$transcript\n";
}



# FUNCTIONS
sub read_json {
    my $file = shift;
    my $json = "";

    open my $fh, "<", $file
        or die "Can't open file : $!\n";
    while (<>) {
        chomp;
        $json .= $_;
    }
    my $exp_data = decode_json($json);
    return $exp_data;
}
