#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use REST::Neo4p;
use Data::Dumper;
use IO::Zlib;

my %OPTS = ();
my $IP   = "http://192.168.0.2:7474";


# PARSE COMMAND LINE ARGUMENTS
GetOptions (
    \%OPTS         ,
    'help|?'       ,
    'human=s'      ,
    'species=s'    ,
    'fasta=s'      ,
    'length=s'     ,
    'pfamdb=s'     ,
    'dir=s'
 );


my $ORF_FILE      = $OPTS{'dir'} . "/ORFS/longestorfs.fa";
my $INT_FILE      = $OPTS{'dir'} . "/INTERACTOME.tbl";
my $PFAM_FILE     = $OPTS{'dir'} . "/PFAM/transcriptome_pfam.fix.tbl";
my $HOMOLOGY_FILE = $OPTS{'dir'} . "homology_infojoiner.tbl";
my $PFAM_DBFILE   = $OPTS{'pfamdb'};

# START TIME
my $init_time  = localtime();
my $start_time = time();
print STDERR <<EOF
#
####### PROGRAM STARTED #######
# Start time: $init_time
# Connecting to $IP ...

EOF
;


die unless $OPTS{'species'};

connect_to_db();
if ($OPTS{'human'}) {
    upload_human($OPTS{'human'});
}

my %seq_info = ();
read_fasta(\%seq_info, $OPTS{fasta}, "sequence");
read_fasta(\%seq_info, $ORF_FILE,    "orf");
upload_transcriptome(\%seq_info);

# Right now we can't upload homology information for all the sequences because the file
# is too big. We should create a "simpler" file with info_joiner so that we can use it here.
#upload_homology($HOMOLOGY_FILE);

upload_interactome(\%seq_info, $INT_FILE);


if ($PFAM_DBFILE) {
    my %pfam_info = ();
    read_pfam(\%pfam_info, $PFAM_DBFILE);
    upload_PFAM(\%pfam_info, $OPTS{'species'}, $PFAM_FILE);
}

my $end_time = localtime();
print STDERR <<EOF

#
####### PROGRAM FINISHED #######
# Finish time: $end_time

EOF
;


#===============================================================================
# FUNCTIONS
#===============================================================================

sub connect_to_db {
    eval {
    	REST::Neo4p->connect("$IP", "neo4j", "5961");
    };

    if ($@) {
    	die "Can't connect to NEO4j database at $IP: \n$@\n";
    };
    return;
}

#--------------------------------------------------------------------------------
sub upload_human {
    my $file = shift;

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    my $node_q = node_query("Human");
    my $edge_q = edge_query("INTERACT_WITH");


    while (<$fh>) {
        chomp;
        my @n = split /\t/;
        add_nodes(\@n, $node_q);
        add_edges($n[0], $n[1], $edge_q);
    }
    return;
}

#--------------------------------------------------------------------------------
sub node_query {
    my $sp = shift;
    my $stmt  = "MERGE (n:$sp {symbol: {gene}})";
    my $node_q = REST::Neo4p::Query->new($stmt);

    return $node_q
}

#--------------------------------------------------------------------------------
sub edge_query {
    my $type = shift;
    my $stmt    = "MATCH (n{symbol: {source}}), (b{symbol: {target}}) " .
		 	      "MERGE (n)-[:$type]->(b)";
    my $edge_q = REST::Neo4p::Query->new($stmt);

    return $edge_q;
}

#--------------------------------------------------------------------------------
sub add_nodes {
    my $nodes = shift;
    my $query = shift;

    print STDERR "Adding ", scalar(@{$nodes}), " nodes...\n";
    foreach my $node (@{$nodes}) {
        $query->execute(gene => $node);
        print STDERR "\tnode: $node uploaded.\n";
    }
    return;
}

#--------------------------------------------------------------------------------
sub add_edges {
    my $source = shift;
    my $target = shift;
    my $query  = shift;

    print STDERR "Adding interaction...\n";
    $query->execute(
        source => $source,
        target => $target
    );

    return;
}

#--------------------------------------------------------------------------------
sub read_fasta {
    my $seq_info = shift;
    my $file     = shift;
    my $att_name = shift;

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    local $/ = ">";
    <$fh>; #skip 1st line
    while (<$fh>) {
        chomp;
        my ($id, @seq) = split /\n/;
        $seq_info->{$id}->{$att_name} = join("", @seq);

        if ($att_name eq "sequence") {
            $seq_info->{$id}->{length} = length(join("", @seq));
        }
    }

    return;
}

#--------------------------------------------------------------------------------
sub node_query_pred {
    my $sp = shift;
    my $stmt  = "MERGE (n:$sp {symbol: {gene}, sequence: {seq}, orf: {orf}, length: {length}})";
    my $node_q = REST::Neo4p::Query->new($stmt);

    return $node_q;
}

#--------------------------------------------------------------------------------
sub add_nodes_pred {
    my $nodes = shift;
    my $query = shift;
    my $seq_info = shift;

    print STDERR "Adding ", scalar(@{$nodes}), " nodes...\n";
    foreach my $node (@{$nodes}) {
        $query->execute(
            gene   => $node,
            seq    => $seq_info->{$node}->{sequence},
            orf    => $seq_info->{$node}->{orf},
            length => $seq_info->{$node}->{length}
        ) or die $query->errstr, "\n";
        print STDERR "\tnode: $node uploaded.\n";
    }
    return;
}

#--------------------------------------------------------------------------------
sub edge_query_pred {
    my $stmt    = "MATCH (n{symbol: {source}}), (b{symbol: {target}}) " .
		 	      "MERGE (n)-[r:INTERACT_WITH {
                      path_length: {plen},
                      dom_int_sc:  {domint},
                      molfun_nto:  {molfun},
                      bioproc_nto: {bioproc},
                      cellcom_nto: {cellcom},
                      int_prob   : {intprob}
                   }]->(b)";
    my $edge_q = REST::Neo4p::Query->new($stmt);

    return $edge_q;
}

#--------------------------------------------------------------------------------
sub add_edges_pred {
    my $source   = shift;
    my $target   = shift;
    my $query    = shift;
    my $features = shift;

    print STDERR "\n\tAdding Pred interaction...\n";
    $query->execute(
        source  => $source,
        target  => $target,
        plen    => $features->[0],
        domint  =>   $features->[15],
        molfun  =>   $features->[16],
        bioproc => $features->[17],
        cellcom => $features->[18],
        intprob => $features->[20]
    ) or die $query->errstr, "\n";;

    return;
}

#--------------------------------------------------------------------------------
sub edge_query_hom {
    my $type = shift;
    my $stmt    = "MATCH (n{symbol: {source}}), (b{symbol: {target}}) " .
		 	      "MERGE (n)-[r:$type {
                      blast_eval: {blast_eval},
                      blast_cov:  {blast_cov},
                      blast_brh:  {blast_brh},
                      nog_eval:   {nog_eval},
                      nog_brh:    {nog_brh},
                      pfam_sc:    {pfam_sc},
                      pfam_brh:   {pfam_brh}
                  }]->(b)";
    my $edge_q = REST::Neo4p::Query->new($stmt);

    return $edge_q;
}

#--------------------------------------------------------------------------------
sub add_edges_hom {
    my $source = shift;
    my $target = shift;
    my $query  = shift;
    my $blast_eval = shift;
    my $blast_cov = shift;
    my $blast_brh = shift;
    my $nog_eval = shift;
    my $nog_brh = shift;
    my $pfam_sc = shift;
    my $pfam_brh = shift;

    print STDERR "\tAdding homolog...\n";
    $query->execute(
        source => $source,
        target => $target,
        blast_eval => $blast_eval,
        blast_cov => $blast_cov,
        blast_brh => $blast_brh,
        nog_eval => $nog_eval,
        nog_brh => $nog_brh,
        pfam_sc => $pfam_sc,
        pfam_brh => $pfam_brh
    ) or die $query->errstr, "\n";

    return;
}

#--------------------------------------------------------------------------------
sub upload_transcriptome {
    my $seq_info = shift;
    my $query = node_query_pred($OPTS{species});
    add_nodes_pred([keys %{ $seq_info }], $query, $seq_info);
    return;
}

#--------------------------------------------------------------------------------
sub upload_homology {
    my $homology_file = shift;
    my $hom_query  = edge_query_hom("HOMOLOG_OF");
    open my $fh, "<", $homology_file
        or die "Can't open $homology_file : $!\n";

    my %added = ();
    my $first_line = <$fh>;
    my $n = 1;
    while (<$fh>) {
        chomp;
        $n++;
        my ($tr1, $tr2, $hom1, $hom2, @features) = split /\t/;
        add_edges_hom($tr1, $hom1, $hom_query, $features[1], $features[3], $features[5], $features[7], $features[9],  $features[11], $features[13]);
        add_edges_hom($tr2, $hom2, $hom_query, $features[2], $features[4], $features[6], $features[8], $features[10], $features[12], $features[14]);
    }
    return;
}

#--------------------------------------------------------------------------------
sub upload_interactome {
    my $seq_info   = shift;
    my $file       = shift;
    my $node_query = node_query_pred($OPTS{species});
    my $edge_query = edge_query_pred();
    my $hom_query  = edge_query_hom("HOMOLOG_OF");

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";
    my $header = <$fh>;
    my ($h_tr1, $h_tr2, $h_hom1, $h_hom2, @h_features) = split /\t/, $header;
    my $n = 1;
    while (<$fh>) {
        chomp;
        print STDERR "Interaction number $n...\n";
        $n++;
        my ($tr1, $tr2, $hom1, $hom2, @features) = split /\t/, $_;
        # All nodes should be already added!!!!
        #add_nodes_pred([$tr1, $tr2], $node_query, $seq_info);
        add_edges_hom($tr1, $hom1, $hom_query, $features[1], $features[3], $features[5], $features[7], $features[9],  $features[11], $features[13]);
        add_edges_hom($tr2, $hom2, $hom_query, $features[2], $features[4], $features[6], $features[8], $features[10], $features[12], $features[14]);
        add_edges_pred($tr1, $tr2, $edge_query, \@features);
    }

    return;

}

#--------------------------------------------------------------------------------
sub read_pfam {
    # Reads the Pfam-A.hmm.dat.gz file and saves the info about the domains.
    my $pfam_info = shift;
    my $file      = shift;


    local $/ = "//";
    open my $fh, "gzip -dc $file |"
        or die "Can't open $file :$!\n";

    DOMAIN:
    while (<$fh>) {
        chomp;
        my @fields = split /\n/;
        my %data = ();
        LINE:
        foreach my $field (@fields) {
            next LINE unless $field =~ m/^#=GF/;
            my ($junk, $entry_name, @value) = split /\s+/, $field; # Value is array because it may have spaces.
            next LINE unless $entry_name =~ m/ID|DE|AC|ML/;
            $data{$entry_name} = join(" ", @value);
        }
        if (%data) {
            $pfam_info->{ $data{"AC"} } = () unless exists $pfam_info->{$data{"AC"}};
            $pfam_info->{ $data{"AC"} }->{"ID"} = $data{"ID"};
            $pfam_info->{ $data{"AC"} }->{"DE"} = $data{"DE"};
            $pfam_info->{ $data{"AC"} }->{"ML"} = $data{"ML"};
        }
    }
    close $fh;

    return;
}


#--------------------------------------------------------------------------------
sub upload_PFAM {
    # Reads PFAM file and uploads the necessary PFAM domains and edges.
    my $pfam_info = shift;
    my $species   = shift;
    my $pfam_file = shift;

    my $dom_stmt = "MERGE (dom:Pfam { accession: {acc}, description: {desc}, identifier: {ident}, mlength: {mlen} })";
    my $al_stmt  = "MATCH (n:$species { symbol: {symbol} }), (dom:Pfam { accession: {acc} })
                    MERGE (n)-[r:HAS_DOMAIN {
                        pfam_start: {pstart},
                        pfam_end:   {pend},
                        s_start:    {sstart},
                        s_end:      {send},
                        perc:       {perc}
                    }]->(dom)
    ";

    my $dom_query = REST::Neo4p::Query->new($dom_stmt);;
    my $al_query  = REST::Neo4p::Query->new($al_stmt);

    open my $fh, "<", $pfam_file
        or die "Can't open $pfam_file : $!\n";

    print STDERR "\n\n#Adding PFAM Domains\n";
    my $n = 1;
    while (<$fh>) {
        chomp;
        my ($seq, @domains) = split /\t/;
        foreach my $domain (@domains) {
            print "\tUploading domain number $n\n";
            $n++;
            my ($acc, $seq_coord, $pfam_coord, $perc) = split /\s+/, $domain;
            my ($pfam_start, $pfam_end) = split /\-/, $pfam_coord;
            my ($s_start, $s_end)       = split /\-/, $seq_coord;

            my $ident = $pfam_info->{$acc}->{"ID"};
            my $desc  = $pfam_info->{$acc}->{"DE"};
            my $mlen  = $pfam_info->{$acc}->{"ML"};

            $dom_query->execute(
                acc   => $acc,
                desc  => $desc,
                ident => $ident,
                mlen  => $mlen
            ) or die $dom_query->errstr, "\n";;

            $al_query->execute(
                symbol => $seq,
                acc    => $acc,
                pstart => $pfam_start,
                pend   => $pfam_end,
                sstart => $s_start,
                send   => $s_end,
                perc   => $perc
            ) or die $al_query->errstr, "\n";
        }
    }
}
