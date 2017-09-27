#!/usr/bin/perl

use warnings;
use strict;
use Graph::Directed;
use Data::Dumper;
# Prints an adjacency table of GO terms, like:
# GO [tab] Parent1 [tab] Parent2 ...
#
#

sub read_go {
    my $file = shift;
    my %graphs = (
     "biological_process" => Graph::Directed->new,
     "molecular_function" => Graph::Directed->new,
     "cellular_component" => Graph::Directed->new
    );

    my %roots = (
        "GO:0003674"  => "molecular function",
        "GO:0008150"  => "biological process",
        "GO:0005575 " => "cellular component"
    );

    local $/ = "[Term]";

    open my $fh, "<", $file
        or die "Can't open $file :$!\n";

    <$fh>; # Skip big header
    while (<$fh>) {
        chomp;
        my @lines = split /\n/;
        @lines = grep {/^id|^is_a|^is_obsolete|^namespace/} @lines;
        if (grep {/^is_obsolete/} @lines) { next; }
        my ($id, $namespace, @is_a) = (@lines);
        $namespace =~ s/namespace\:\s?//g;
        my $jnk;
        ($jnk, $id) = split /id\: /, $id;
        next unless $id =~ m/^GO\:/;
        @is_a = map {
            my ($go_parent) = split /\s?!\s?/;
            $go_parent =~ s/is_a: //;
            $go_parent if not exists $roots{$go_parent} and $go_parent =~ m/^GO/ # Remove roots
        } @is_a;

        foreach my $parent (@is_a) {
            next unless $id =~ m/^GO\:/ and $parent =~ m/^GO\:/;
            $graphs{$namespace}->add_edge($parent, $id);
        }

    }

    return \%graphs;


}

## MAIN
my $file  = shift @ARGV;
my $graphs = read_go($file);

foreach my $namespace ("biological_process","molecular_function","cellular_component") {
    open my $fh, ">", $namespace . "_adjacency.tbl"
        or die "Can't write file :$!\n";

    foreach my $v ($graphs->{$namespace}->vertices) {
        my @parent_gos = $graphs->{$namespace}->all_predecessors($v);
        print $fh "$v\t", join("\t", @parent_gos), "\n";
    }
    close $fh;
}
