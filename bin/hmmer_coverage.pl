#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $domains_file    = shift @ARGV;
my %best_queries    = ();
my %best_targets    = ();
my %scores          = ();
my %lengths         = ();
my %reciprocals     = ();
my %results         = ();
my $hsp_count       = 0;

print "SEQUENCE", "\t",
	  "MODEL", "\t",
	  "SCORE", "\t",
	  "SEQLENGTH", "\t",
	  "MODELLENGTH", "\t",
	  "ALIGNLENGTH", "\t",
	  "ALHMMERLENGTH", "\t",
	  "SEQCOVERAGE", "\t",
	  "HMMERCOVERAGE", "\t",
	  "SEQBEST", "\t",
	  "MODELBEST", "\t",
	  "RECIPROCAL\t",
	  "EVALUE\n";



open my $dom_fh, "<", $domains_file
	or die "Can't open $domains_file : $!\n";


while (<$dom_fh>) {
	chomp;
	next if /^#/;
	$hsp_count++;
	my @fields = split /\s+/;
	my ($cthulhu, $cth_length, $nogg, $model_length, $eval, $h_from, $h_to, $al_from, $al_to) = (
		$fields[0],
		$fields[2],
		$fields[3],
		$fields[5],
		$fields[6],
		$fields[15],
		$fields[16],
		$fields[17],
		$fields[18]
		);
	$nogg =~ s/meNOG\.(.+)\..+/$1/;
	my $score = $fields[7];

	# Is it a best hit?

	if (not exists $best_queries{$cthulhu}) {
		$best_queries{$cthulhu}->{score} = $score;
		$best_queries{$cthulhu}->{evalue} = $eval;
		$best_queries{$cthulhu}->{target}->{$nogg} = undef;
	} else {
		if ($eval < $best_queries{$cthulhu}->{evalue}) {
			$best_queries{$cthulhu}->{score} = $score;
			$best_queries{$cthulhu}->{score} = $eval;
			$best_queries{$cthulhu}->{target} = undef;
			$best_queries{$cthulhu}->{target}->{$nogg} = undef;
		}
	}

	if (not exists $best_targets{$nogg}) {
		$best_targets{$nogg}->{score} = $score;
		$best_targets{$nogg}->{evalue} = $eval;
		$best_targets{$nogg}->{query}->{$cthulhu} = undef;
	} else {
		if ($eval < $best_targets{$nogg}->{evalue}) {
			$best_targets{$nogg}->{score} = $score;
			$best_targets{$nogg}->{evalue} = $eval;
			$best_targets{$nogg}->{query} = undef;
			$best_targets{$nogg}->{query}->{$cthulhu} = undef;
		}
	}



	# SAVE RESULTS
	$results{$cthulhu}->{$nogg} = undef unless exists $results{$cthulhu}->{$nogg};


	$results{$cthulhu}->{$nogg}->{SCORE}       = $score;
	$results{$cthulhu}->{$nogg}->{EVALUE}      = $eval;
	$results{$cthulhu}->{$nogg}->{SEQLENGTH}   = $cth_length;
	$results{$cthulhu}->{$nogg}->{MODELLENGTH} = $model_length;

	$results{$cthulhu}->{$nogg}->{ALINTERVALS} = []
		unless exists $results{$cthulhu}->{$nogg}->{ALINTERVALS};
	push @{ $results{$cthulhu}->{$nogg}->{ALINTERVALS} }, [$al_from, $al_to];

	$results{$cthulhu}->{$nogg}->{HMMINTERVALS} = []
		unless exists $results{$cthulhu}->{$nogg}->{HMMINTERVALS};
	push @{ $results{$cthulhu}->{$nogg}->{HMMINTERVALS} }, [$h_from, $h_to];


}


# Get reciprocals
foreach my $query (keys %best_queries) {
	my ($target) = keys %{$best_queries{$query}->{target} };

	if (exists $best_targets{$target}->{query}->{$query}) {
		$reciprocals{$query}->{$target} = undef;
	}

}


# MERGE INTERVALS AND PRINT RESULTS
foreach my $cthulhu (keys %results) {
	foreach my $nogg (keys %{ $results{$cthulhu} }) {

		my @al_starts  =  ();
		my @al_ends    =  ();
		my @hmm_starts =  ();
		my @hmm_ends   =  ();

		my @al_intervals  = @{ $results{$cthulhu}->{$nogg}->{ALINTERVALS} };
		my @hmm_intervals = @{ $results{$cthulhu}->{$nogg}->{HMMINTERVALS} };

		@al_intervals = sort {
			$a->[0] <=> $b->[0]
		} @al_intervals;

		@hmm_intervals = sort {
			$a->[0] <=> $b->[0]
		} @hmm_intervals;




		# ALIGNMENT LENGTH (ON CTHULHU SEQUENCE)
		foreach my $interval (@al_intervals) {
			push @al_starts, $interval->[0];
			push @al_ends, $interval->[1];
		}

		for (1..$#al_starts) {
    		# extra check on array bounds, since we edit in-place
    		last unless $_ < @al_starts;
    		# don't need to collapse if no overlap with previous end
    		next unless $al_starts[$_] <= $al_ends[$_-1];
    		# delete this start and the previous end
    		splice(@al_starts,$_,1);
    		splice(@al_ends,$_-1,1);
    		# rerun this loop for the same value of $_ since it was deleted
    		redo;
		}

		my $al_length = 0;
		for my $i (0..$#al_starts) {
			$al_length += $al_ends[$i] - $al_starts[$i];
		}


		# ALIGNMENT LENGTH (ON HMMER MODEL)
		foreach my $interval (@hmm_intervals) {
			push @hmm_starts, $interval->[0];
			push @hmm_ends, $interval->[1];
		}

		for (1..$#hmm_starts) {
    		# extra check on array bounds, since we edit in-place
    		last unless $_ < @hmm_starts;
    		# don't need to collapse if no overlap with previous end
    		next unless $hmm_starts[$_] <= $hmm_ends[$_-1];
    		# delete this start and the previous end
    		splice(@hmm_starts,$_,1);
    		splice(@hmm_ends,$_-1,1);
    		# rerun this loop for the same value of $_ since it was deleted
    		redo;
		}

		my $hmm_length = 0;
		for my $i (0..$#hmm_starts) {
			$hmm_length += $hmm_ends[$i] - $hmm_starts[$i];
		}

		my $seqbest       = 0;
		my $modelbest     = 0;
		my $reciprocal    = 0;
		my $seqcoverage   = $al_length / $results{$cthulhu}->{$nogg}->{SEQLENGTH};
		my $hmmercoverage = $hmm_length / $results{$cthulhu}->{$nogg}->{MODELLENGTH};

		if (exists $best_queries{$cthulhu}->{target}->{$nogg}) {
			$seqbest = 1;
		}

		if (exists $best_targets{$nogg}->{query}->{$cthulhu}) {
			$modelbest = 1;
		}

		if (exists $reciprocals{$cthulhu} and exists $reciprocals{$cthulhu}->{$nogg}) {
			$reciprocal = 1;
		}

		print "$cthulhu", "\t",
		      "$nogg", "\t",
		      $results{$cthulhu}->{$nogg}->{SCORE}, "\t",
		      $results{$cthulhu}->{$nogg}->{SEQLENGTH}, "\t",
		      $results{$cthulhu}->{$nogg}->{MODELLENGTH}, "\t",
		      "$al_length", "\t",
		      "$hmm_length", "\t",
		      "$seqcoverage", "\t",
		      "$hmmercoverage", "\t",
		      "$seqbest", "\t",
		      "$modelbest", "\t",
		      "$reciprocal", "\t",
			  $results{$cthulhu}->{$nogg}->{EVALUE}, "\n";

	}
}





# Print some numbers!

print STDERR "\n\n",
			 "# TOTAL HSPs   : ", $hsp_count, "\n",
             "# TOTAL QUERIES: ", scalar(keys %best_queries), "\n",
             "# TOTAL TARGETS: ", scalar(keys %best_targets), "\n",
             "# RECIPROCALS  : ", scalar(keys %reciprocals), "\n",
			 "\n\n";
