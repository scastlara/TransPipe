#!/usr/bin/perl
use warnings;
use strict;
use REST::Neo4p;


my $go_file = shift @ARGV;
my $IP   = "http://192.168.0.2:7474";

sub connect_to_db {
    eval {
    	REST::Neo4p->connect("$IP", "neo4j", "5961");
    };

    if ($@) {
    	die "Can't connect to NEO4j database at $IP: \n$@\n";
    };
    return;
}

connect_to_db();

open my $fh, "<", $go_file
    or die "Can't open $go_file : $!\n";

my $stmt    = "MATCH (n:Human {symbol: {gene}}) RETURN n.symbol AS symbol";
my $node_q  = REST::Neo4p::Query->new($stmt);
my $go_stmt = "MATCH (n:Human {symbol: {gene}})
               MERGE (n)-[r:HAS_GO]->(go:Go {accession: {accession}, domain: {domain}})";
my $go_q    = REST::Neo4p::Query->new($go_stmt);

<$fh>; # skip first line
while (<$fh>) {
    chomp;
    my ($ens, $hgnc, $symbol, $domain, $go) = split /\t/;
    next unless ($symbol and $domain and $go);
    $node_q->execute(gene => $symbol);

    while (my $row = $node_q->fetch) {
        print STDERR "SYMBOL: $symbol\t";
        $go_q->execute(
            gene => $symbol,
            accession => $go,
            domain    => $domain
        );
        print STDERR "\n";
    }
}
