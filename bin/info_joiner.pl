#!/usr/bin/perl
# CASES:
    # Col1: Trans/Gene 1
    # Col2: Trans/Gene 2
    # Col3: Homolog 1
    # Col4: Homolog 2
# FEATURES
    # Col5: H1 -> H2 distance in human graph
    # Col6: Nog evalue 1
    # Col6: Nog evalue 2
    # Col7: Nog BRH 1
    # Col8: Nog BRH2
    # Col9: Blast evalue 1
    # Col10: Blast evalue 2
    # Col11: Blast BRH 1
    # Col12: Blast BRH 2
    # Col 13: Blast Coverage 1
    # Col 14: Blast Coverage 2
    # Col 15: Pfam BRH 1
    # Col 16: Pfam BRH 2
    # Col 17: Pfam Coverage 1
    # Col 18: Pfam Coverage 2
    # Col 19: Domain interaction probability score
    # Col 20: GO molfunction distance
    # Col 21: GO bioprocess distance
    # Col 22: GO cellcomponent distance
    # Col 23: TRUE/FALSE


use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);;
use Data::Dumper;
use List::Util qw[min max];
use PerlIO::gzip;


#===============================================================================
# VARIABLES AND OPTIONS
#===============================================================================

my %OPTS;
GetOptions (
    \%OPTS        ,
    'blast=s'     ,
    'pfam=s'      ,
    'eggnogs=s'   ,
    'rh_blast=s'  ,
    'rh_pfam=s'   ,
    'rh_nogs=s'   ,
    'network=s'   ,
    'dict=s'      ,
    'train_dict=s',
    'plen=s'      ,
    'Pfam_ints=s' ,
    'domains=s'   ,
    'go_annot=s'  ,
    'go_parents=s'  ,
    'hgnc_2_symbol=s' ,
    'Train'     , # Can be T or F or empty
 );

#===============================================================================
# MAIN
#===============================================================================

print STDERR "\n\n# Reading alignments data and computing Best Hits...";
my %homo_data = ();
read_blast($OPTS{'blast'}, \%homo_data);
read_pfam( $OPTS{'pfam'},  \%homo_data);
my $nog_dict = read_nogs( $OPTS{'eggnogs'},  \%homo_data, $OPTS{'dict'});

foreach my $method ("blast", "pfam") {
    my $str = "rh_" . $method;
    read_brh($method, $OPTS{$str}, \%homo_data, $nog_dict);
}

my $bestest_hits = select_bestest_hits(\%homo_data);
fill_bestdata($bestest_hits);
print STDERR "done.\n\n";

my $train_dict;
if ($OPTS{Train}) {
    $train_dict     = read_train_dict($OPTS{train_dict});
}

my $hgnc_2_symbol  = read_hgnc_alias($OPTS{hgnc_2_symbol});

print STDERR "# Reading Shortest paths...\n";
my $plengths       = read_shortestpaths($OPTS{plen});
print STDERR "\ndone.\n\n";

print STDERR "# Reading GO data...";
my $go_annotations = read_go_annotation($OPTS{go_annot});
my $go_parents     = read_go_parents($OPTS{go_parents});
print STDERR "done.\n\n";

print STDERR "# Reading PFAM data...";
my $pfam_ints      = read_pfam_ints($OPTS{Pfam_ints});
my $domains        = read_domains($OPTS{domains});
print STDERR"done.\n\n";


print STDERR "# Iterating through all possible interactions...";
print_header();

my $done_hash = ();
my $num_int = 0;

Q1:
foreach my $query1 (keys %{ $bestest_hits }) {
    Q2:
    foreach my $query2 (keys %{ $bestest_hits }) {
        my $p1 = $OPTS{Train} ? $train_dict->{$query1} : $query1;
        my $p2 = $OPTS{Train} ? $train_dict->{$query2} : $query2;
        my ($hom1) = keys %{ $bestest_hits->{$query1} };
        my ($hom2) = keys %{ $bestest_hits->{$query2} };
        my $symbol1 = $hgnc_2_symbol->{$hom1};
        my $symbol2 = $hgnc_2_symbol->{$hom2};


        # Skip repeated interactions A -> B == B -> A
        next if exists $done_hash->{$p1}->{$p2} or exists $done_hash->{$p2}->{$p1};
        $done_hash->{$p1}->{$p2} = 1;

        # Get path length between hom1 and hom2 in human interactome
        my $path_length = "";
        # Both genes in interactome!
        if (exists $plengths->{$symbol1}->{$symbol2}) {
            # There is a path between them! Even if it's -1/Inf
            $path_length = $plengths->{ $symbol1 }->{ $symbol2 };
        } elsif (exists $plengths->{$symbol2}->{$symbol1}) {
            # There is a path between them! Even if it's -1/Inf
            $path_length = $plengths->{ $symbol2 }->{ $symbol1 };
        } else {
            # One of the nodes (or both) is not in the graph\n";
            $path_length = "NA";
            next Q2;
        }

        # Skip pairs too far away in human interactome
        next Q2 unless $path_length <= 2;
        next Q2 if     $path_length == -1;


        # Compute GO similarities
        my $molfun_nto  = compute_nto("molecular_function", $go_parents, $go_annotations, $hom1, $hom2);
        my $bioproc_nto = compute_nto("biological_process", $go_parents, $go_annotations, $hom1, $hom2);
        my $cellcom_nto = compute_nto("cellular_component", $go_parents, $go_annotations, $hom1, $hom2);
        #qprint STDERR "$molfun_nto and $bioproc_nto and $cellcom_nto\n";

        # Skip pairs with small cellcomp nto
        #next Q2 if $cellcom_nto == 0;

        # Compute pfam ints score
        my $pfam_int_score = compute_pfam_score($pfam_ints, $domains, $query1, $query2);

        # Print results
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            $p1,
            $p2,
            $symbol1,
            $symbol2,
            $path_length,
            $bestest_hits->{$query1}->{$hom1}->{blast}->{evalue},
            $bestest_hits->{$query2}->{$hom2}->{blast}->{evalue},
            $bestest_hits->{$query1}->{$hom1}->{blast}->{cov},
            $bestest_hits->{$query2}->{$hom2}->{blast}->{cov},
            $bestest_hits->{$query1}->{$hom1}->{blast}->{BRH},
            $bestest_hits->{$query2}->{$hom2}->{blast}->{BRH},
            $bestest_hits->{$query1}->{$hom1}->{nogs}->{evalue},
            $bestest_hits->{$query2}->{$hom2}->{nogs}->{evalue},
            $bestest_hits->{$query1}->{$hom1}->{nogs}->{BRH},
            $bestest_hits->{$query2}->{$hom2}->{nogs}->{BRH},
            $bestest_hits->{$query1}->{$hom1}->{pfam}->{cov},
            $bestest_hits->{$query2}->{$hom2}->{pfam}->{cov},
            $bestest_hits->{$query1}->{$hom1}->{pfam}->{BRH},
            $bestest_hits->{$query2}->{$hom2}->{pfam}->{BRH},
            $pfam_int_score,
            $molfun_nto,
            $bioproc_nto,
            $cellcom_nto
        );

    }
}

print STDERR "done.\n\n";


#===============================================================================
# FUNCTIONS
#===============================================================================
#--------------------------------------------------------------------------------
sub read_blast {
    my $file      = shift;
    my $homo_data = shift;
    my %q_len     = ();

    open my $fh, '<', $file
        or die "Can't open $file : $!\n";

    while (<$fh>) {
        chomp;
        next unless /^Q/;
        my @cols = split /\s+/;

        # Remove additional names/description from subject
        $cols[13] =~ s/\|.+//;
        my $query  = $cols[1];
        my $target = $cols[13];
        my $cov    = $cols[12];
        my $evalue = $cols[16];

        # Save data
        $homo_data->{$query}->{$target}->{'blast'}->{'evalue'} = $evalue;
        $homo_data->{$query}->{$target}->{'blast'}->{'cov'}    = $cov;
        $homo_data->{$query}->{$target}->{'blast'}->{'BRH'}    = 0;
    }

    return;
}


#--------------------------------------------------------------------------------
sub read_pfam {
    my $file      = shift;
    my $homo_data = shift;

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    while (<$fh>) {
        chomp;
        next unless /^C/;
        my ($j, $query, $target, $cov) = split /\t/;

        $homo_data->{$query}->{$target}->{'pfam'}->{'cov'} = $cov;
        $homo_data->{$query}->{$target}->{'pfam'}->{'BRH'} = 0;
    }

    return;
}

#--------------------------------------------------------------------------------
sub read_nogs {
    my $file        = shift;
    my $homo_data   = shift;
    my $dict_file   = shift;
    my $nog_to_hgnc = get_dict($dict_file);
    my $best_hits   = ();

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    <$fh>;
    while (<$fh>) {
        chomp;
        my @cols = split /\t/;
        next unless $cols[9]; # It's a best HIT
        my ($tr, $enog, $evalue) = ($cols[0], $cols[1], $cols[-1]);
        foreach my $hgnc (@{ $nog_to_hgnc->{$enog} }) {
            $homo_data->{$tr}->{$hgnc}->{'nogs'}->{'evalue'} = $evalue;
            if ($cols[-2]) {
                $homo_data->{$tr}->{$hgnc}->{'nogs'}->{'BRH'} = 1;
            } else {
                $homo_data->{$tr}->{$hgnc}->{'nogs'}->{'BRH'} = 0;
            }

        }
    }

    return $nog_to_hgnc;
}

#--------------------------------------------------------------------------------
sub get_dict {
    my $file         = shift;
    my %nogs_to_hgnc = ();

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    while (<$fh>) {
        chomp;
        my ($enog, $ensp, $hgncs) = split /\t/;
        next if $hgncs eq "NA";
        my @hgncs = split /\s/, $hgncs;

        foreach my $h (@hgncs) {
            push @{ $nogs_to_hgnc{$enog} }, $h  ;
        }

    }
    return \%nogs_to_hgnc;

}

#--------------------------------------------------------------------------------
sub read_brh {
    my $method      = shift;
    my $file        = shift;
    my $homo_data   = shift;
    my $nog_to_hgnc = shift;

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    <$fh>; # skip header
    while (<$fh>) {
        chomp;
        my ($q, $t) = split /\s+/;
        $homo_data->{$q}->{$t}->{$method}->{'BRH'} = 1;

    }

    return;
}

#--------------------------------------------------------------------------------
sub read_train_dict {
    my $file = shift;
    my %dict = ();

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    while (<$fh>) {
        chomp;
        my ($tr, $gene) = split /\t/;
        $dict{$tr} = $gene;
    }
    return \%dict;
}


#--------------------------------------------------------------------------------
sub select_bestest_hits {
    my $homo_data = shift;
    my $bestest   = ();

    # Rule:
    #   1- NOGS BRH
    #   2- More keys
    #   3- Is NOG best hit
    #   4- Is BLAST best hit
    #   5- warning

    foreach my $query (keys %{$homo_data}) {
        my %targets = ();

        foreach my $target (keys %{ $homo_data->{$query} }) {
            $targets{$target}->{"numk"} = keys %{ $homo_data->{$query}->{$target} };

            if (exists $homo_data->{$query}->{$target}->{"nogs"}) {
                $targets{$target}->{"has_nog"} = 1;
                $targets{$target}->{"n_brh"} = $homo_data->{$query}->{$target}->{"nogs"}->{"BRH"};
            } else {
                $targets{$target}->{"has_nog"} = 0;
                $targets{$target}->{"n_brh"}   = 0;
            }

            if (exists $homo_data->{$query}->{$target}->{'blast'}) {
                $targets{$target}->{"has_blast"} = 1;
            } else {
                $targets{$target}->{"has_blast"} = 0;
            }
        }

        # Select best target
        my $best_t = "";

        # ------
        # 1st NOG BRH
        my @with_nog_brh = grep {$targets{$_}->{"n_brh"}} keys %targets;

        if (not @with_nog_brh) {
            # No target has NOG BRH, all of them to next rule
            @with_nog_brh = keys %targets;
        } elsif (@with_nog_brh == 1) {
            # THE BEST ONE!
            $bestest->{$query}->{$with_nog_brh[0]} = $homo_data->{$query}->{$with_nog_brh[0]};
            next;
        }

        # ------
        # 2nd More keys
        @with_nog_brh = sort { $targets{$b}->{"numk"} <=> $targets{$a}->{"numk"} } @with_nog_brh;
        my $biggest_num = $targets{$with_nog_brh[0]}->{"numk"}; # this is the first (biggest number of methods)
        my @with_morek  = ();
        foreach my $with (@with_nog_brh) {
            if ($targets{$with}->{"numk"} >= $biggest_num) {
                push @with_morek, $with;
            }
        }

        if (@with_morek == 1) {
            # THE BEST ONE!
            $bestest->{$query}->{$with_morek[0]} = $homo_data->{$query}->{$with_morek[0]};
            next;
        }

        # ------
        # 3rd rule
        my @with_nog = ();
        foreach my $with (@with_morek) {
            if ($targets{$with}->{"has_nog"}) {
                push @with_nog, $with;
            }
        }

        if (not @with_nog) {
            # No NOG! all of them to next rule
            @with_nog = @with_morek;
        } elsif (@with_nog == 1) {
            # THE BEST ONE!
            $bestest->{$query}->{$with_nog[0]} = $homo_data->{$query}->{$with_nog[0]};
            next;
        }

        # ------
        # 4th rule
        my @with_blast = ();
        foreach my $with (@with_morek) {
            if ($targets{$with}->{"has_blast"}) {
                push @with_blast, $with;
            }
        }

        if (not @with_blast) {
            # None have BLAST nor NOG. Since there are no rules, get one randomly.
            # It should get the one with the higher PFAM coverage, right?
            $bestest->{$query}->{$with_nog[0]} = $homo_data->{$query}->{$with_nog[0]};
        } elsif (@with_blast == 1) {
            # THE BEST ONE!
            $bestest->{$query}->{$with_blast[0]} = $homo_data->{$query}->{$with_blast[0]};
        } elsif (@with_blast > 1) {
            # More than one have BLAST and NOG. Get the one with lowest BLAST evalue
            # This shouldn't happen, because BLAST only gives us ONE Best Hit.
            die "PROBLEM. NOT IMPLEMENTED\n";
        }

    } # for each query

    $homo_data = undef;
    return ($bestest);
}

#--------------------------------------------------------------------------------
sub read_shortestpaths {
    my $path_file = shift;
    my %path_length = ();

    if (defined open(my $fh, "<:gzip", $path_file)) {
        my $million_lines = 0;

        while (<$fh>) {
            chomp;
            if ($. % 1000000 == 0) {
                $million_lines++;
                print STDERR "\t# $million_lines Million Lines read\n";
            }
            my ($parent, $child, $length) = split /\t/;
            if ($length == "2147483647") { # This number means no interaction. I can't think of a better rule.
                $path_length{$parent}->{$child} = "-1"; # -1 means Infinite path between nodes
            } else {
                $path_length{$parent}->{$child} = $length;
            }
        }
    } else {
        die "Can't gzip $path_file :$!\n";
    }


    return \%path_length;
}


#--------------------------------------------------------------------------------
sub fill_bestdata {
    # This function fills the data for all the transcripts.
    # Instead of doing this I should initialize it BEFORE the actual data is stored,
    # but since I did this program a long time ago I prefer to do it a posteriori to not
    # destroy anything.

    my $data = shift;

    foreach my $trans (keys %{$data}) {
        my ($hgnc) = keys %{$data->{$trans}};
        foreach my $method ("blast", "nogs", "pfam") {
            next if exists $data->{$trans}->{$hgnc}->{$method};
            if ($method eq "blast") {
                $data->{$trans}->{$hgnc}->{$method}->{"evalue"} = 10;
                $data->{$trans}->{$hgnc}->{$method}->{"cov"}    = 0;
                $data->{$trans}->{$hgnc}->{$method}->{"BRH"}    = 0;
            } elsif ($method eq "pfam") {
                $data->{$trans}->{$hgnc}->{$method}->{"BRH"} = 0;
                $data->{$trans}->{$hgnc}->{$method}->{"cov"} = 0;
            } else {
                $data->{$trans}->{$hgnc}->{$method}->{"evalue"} = 10;
                $data->{$trans}->{$hgnc}->{$method}->{"BRH"}    = 0;
            }
        }
    }

    return;
}

#--------------------------------------------------------------------------------
sub read_hgnc_alias {
    my $file = shift;
    my %hgnc_2_symbol;

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    <$fh>; # skip header
    while (<$fh>) {
        chomp;
        my ($hgnc, $symbol) = split /\t/;
        $hgnc = "HGNC:". $hgnc;
        $hgnc_2_symbol{$hgnc} = $symbol;
    }

    return \%hgnc_2_symbol
}

#--------------------------------------------------------------------------------

sub read_go_annotation {
    my $file = shift;
    my %go_annotations = ();

    open my $fh, "<", $file
        or die "Can't open $file : $!\n";

    <$fh>; # skip header
    while (<$fh>) {
        chomp;
        my ($ens, $hgnc,  $symbol, $dom, $go) = split /\t/;

        # Using HGNC as gene id
        $go_annotations{$dom}->{$hgnc} = () unless exists $go_annotations{$dom}->{$hgnc};
        push @{ $go_annotations{$dom}->{$hgnc} }, $go;

        # Using ENSEMBL as gene id
        $go_annotations{$dom}->{$ens} = () unless exists $go_annotations{$dom}->{$ens};
        push @{ $go_annotations{$dom}->{$ens} }, $go;

    }

    return \%go_annotations;

}


#--------------------------------------------------------------------------------
sub read_go_parents {
    my $folder       = shift;
    my %go_parents = ();

    foreach my $domain ("biological_process", "cellular_component", "molecular_function") {
        my $file = $folder . "/" . $domain . "_adjacency.tbl";
        open my $fh, "<", $file
            or die "Can't open $file : $!\n";

        while (<$fh>) {
            chomp;
            my ($go, @parents) = split /\t/;
            $go_parents{$go} = () unless exists $go_parents{$go};
            push @{ $go_parents{$go} }, @parents;
        }
    }
    return \%go_parents;

}


#--------------------------------------------------------------------------------
sub compute_nto {
    my $domain         = shift;
    my $go_parents     = shift;
    my $go_annotations = shift;
    my $gene_1         = shift;
    my $gene_2         = shift;


    if (not exists $go_annotations->{$domain}->{$gene_1} or
        not exists $go_annotations->{$domain}->{$gene_2}) {
        return "-1";
    }

    my %g1_parents = ();
    my %g2_parents = ();
    foreach my $g1_go (@{ $go_annotations->{$domain}->{$gene_1} }) {
        $g1_parents{$g1_go } = 1;
    }

    foreach my $g1_go (@{$go_annotations->{$domain}->{$gene_2}}) {
        $g2_parents{$g1_go} = 1;
    }

    my $score = scalar(grep {exists $g1_parents{$_}; } keys %g2_parents );

    if ($score == 0) {
        return 0;
    } else {
        return $score / min(scalar(keys %g1_parents), scalar(keys %g2_parents));
    }
}

#--------------------------------------------------------------------------------
sub read_pfam_ints {
    my $file = shift;
    my %pfam_ints = ();

    open my $fh, "<", $file
        or die "Can't open $file :$!\n";

    while (<$fh>) {
        chomp;
        my ($p1, $p2) = split /\t/;
        $pfam_ints{$p1}->{$p2} = 1;
    }

    return \%pfam_ints;
}

#--------------------------------------------------------------------------------

sub read_domains {
    my $file = shift;
    my %pfam_domains = ();

    open my $fh, "<", $file
        or die "Can't open $file :$!\n";

    while (<$fh>) {
        chomp;
        my ($gene, @pfams) = split /\t/;
        @pfams = map { my ($pfam) = split /\s/; $pfam; } @pfams;
        $pfam_domains{$gene} = [@pfams];
    }

    return \%pfam_domains;

}


#--------------------------------------------------------------------------------
sub compute_pfam_score {
    my $pfam_ints = shift;
    my $domains   = shift;
    my $query1    = shift;
    my $query2    = shift;

    if (not exists $domains->{$query1} or not exists $domains->{$query2}) {
        return 0;
    }

    my $score = 0;
    foreach my $dom1 (@{ $domains->{$query1} }) {
        foreach my $dom2 (@{ $domains->{$query2} }) {
            $score++ if exists $pfam_ints->{$dom1}->{$dom2} or exists $pfam_ints->{$dom2}->{$dom1};
        }
    }

    return $score;

}



#--------------------------------------------------------------------------------
sub print_header {
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        "TRANS_1",
        "TRANS_2",
        "HOM_1",
        "HOM_2",
        "PATH_LENGTH",
        "BLAST_EVAL_1",
        "BLAST_EVAL_2",
        "BLAST_COV_1",
        "BLAST_COV_2",
        "BLAST_BRH_1",
        "BLAST_BRH_2",
        "NOG_EVAL_1",
        "NOG_EVAL_2",
        "NOG_BRH_1",
        "NOG_BRH_2",
        "PFAM_COV_1",
        "PFAM_COV_2",
        "PFAM_BRH_1",
        "PFAM_BRH_2",
        "DOMAIN_INT_SCORE",
        "MOLFUN_NTO",
        "BIOPROC_NTO",
        "CELLCOM_NTO"
    );

    return;
}
