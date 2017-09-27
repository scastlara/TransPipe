#!/usr/bin/perl
use strict;
use warnings;
use IO::Uncompress::Gunzip qw($GunzipError);
use Data::Dumper;

my $fasta          = shift @ARGV;
my $pfamdata       = shift @ARGV;
my %uniprot_2_hugo = ();

read_hugos($fasta);

my $hgnc_2_pfam = read_pfam(\%uniprot_2_hugo, $pfamdata);

# Sort pfams by coordinates

foreach my $hgnc (keys %{$hgnc_2_pfam}) {
    # THIS DOESN't WORK! NEEDS FIX. (stupid mistake)
    print "$hgnc\t";
    my @newfams;
    foreach my $pfam (@{ $hgnc_2_pfam->{$hgnc} }) {
        my ($id, $ac, $uni, $coords) = @{ $pfam };
        my ($start, $end) = split /\-/, $coords;
        push @newfams, [$id, $ac, $uni, $start, $end];
    }
    my @sorted = sort {
        $a->[3] <=> $b->[3]
    } @newfams;
    foreach my $pf (@sorted) {
        print $pf->[0], " ", $pf->[1], " ", $pf->[2], " ", "$pf->[3]-$pf->[4]", "\t";
    }
    print "\n";
}



# ---------------------------------
sub read_hugos {
    my $fasta = shift;

    open my $fh, "<", $fasta
        or die "Can't open $fasta : $!\n";

    while (<$fh>) {
        chomp;
        $_ =~ s/^>//g;
        my ($hgnc, $uniprot) = split /\|/;
        next unless ($hgnc && $uniprot);
        $uniprot_2_hugo{$uniprot} = $hgnc;
    }

}




#----------------------------------------------
sub read_pfam {
	my $uniprot_2_hugo = shift;
    my $pfamdata       = shift;
	my %hgnc_2_pfam    = ();
	my %repeated       = ();

    my $fh = IO::Uncompress::Gunzip->new( $pfamdata )
        or die "IO::Uncompress::Gunzip failed: $GunzipError\n";

	local $/ = "//";

	while (<$fh>) {
		chomp;
		my @lines = split "\n", $_;
		my ($pfam_id, $pfam_acc) = qw(NA NA);

		for (my $i = 0; $i < @lines; $i++) { # for line

			if ($lines[$i] =~ m/^#=GF ID/) { # if ID
				my @fields = split /\s+/, $lines[$i];
				$pfam_id = $fields[2];
			} elsif ($lines[$i] =~ m/^#=GF AC/) { # if Accession
				my @fields = split /\s+/, $lines[$i];
				$pfam_acc = $fields[2];

			} elsif ($lines[$i] =~ m/^#=GS/) {
				my ($sym, $protein) = split /\s+/, $lines[$i];
				my ($uniprot, $coords) = split /\//, $protein;

				next unless (exists $uniprot_2_hugo->{$uniprot});
				next if exists ($repeated{$uniprot . $coords . $pfam_acc});
				$repeated{$uniprot . $coords . $pfam_acc} = undef;

				$hgnc_2_pfam{ $uniprot_2_hugo->{$uniprot} } = []
					unless exists $hgnc_2_pfam{ $uniprot_2_hugo->{$uniprot} };
				push @{ $hgnc_2_pfam{ $uniprot_2_hugo->{$uniprot} } }, [$pfam_acc, $pfam_id, $uniprot, $coords];
			} # if

		} # for line

	} # while (for pfam)


	return \%hgnc_2_pfam;
}
