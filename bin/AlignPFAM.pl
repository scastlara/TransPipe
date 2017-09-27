#!/usr/bin/perl
#################################################################################
#                               AlignPFAM.pl									#
#################################################################################

#================================================================================
#        Copyright (C) 2014 - Sergio CASTILLO
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#================================================================================


use strict;
use warnings;
use Data::Dumper;

#================================================================================
# VARIABLES
#================================================================================

my $pfam_length_file   = shift @ARGV;
my $transcriptome_pfam_file  = shift @ARGV;
my $hgnc_pfam_file     = shift @ARGV;
my %domains            = ();

#================================================================================
# MAIN LOOP
#================================================================================
my $pfam_length = get_pfamlength($pfam_length_file);

my $c_flag = 1;
my $transcriptome_profiles = read_file(
	$transcriptome_pfam_file,
	\%domains,
	$c_flag
);

$c_flag = 0;
my $hgnc_profiles = read_file(
	$hgnc_pfam_file,
	\%domains,
	$c_flag
);

align_pfams_hgnc(
	$hgnc_profiles,
	$transcriptome_profiles,
	\%domains
);

align_pfams_transcriptome(
	$transcriptome_profiles,
	$hgnc_profiles,
	\%domains
);


#================================================================================
# FUNCTIONS
#================================================================================
#--------------------------------------------------------------------------------
sub get_pfamlength {
	my $filename = shift;
	my %length   = ();

	open my $l_fh, '<', $filename
		or die "Can't open $filename : $!\n";

	while (<$l_fh>) {
		chomp;
		my ($pfam, $l) = split /\t/, $_;
		next unless $pfam;
		$length{$pfam} = $l;
	}

	return \%length;
}

#--------------------------------------------------------------------------------
sub read_file {
	my $filename         = shift;
	my $domains          = shift;
	my $c_flag           = shift;
	my %profiles = ();
	my $letter = "";

	open my $cth_fh, '<', $filename
		or die "Can't open $filename : $!\n";

	while (<$cth_fh>) {
		chomp;
		my ($seq_id, @pfams) = split /\t/, $_;

		foreach my $pfam_record (@pfams) {
			my ($p_id,     $s_coords,
				$p_coords, $perc,
				$perc_1,
				$s_start,  $s_end,
				$p_start,  $p_end
				);

			if ($c_flag) {
				$letter = "C";
				($p_id, $s_coords, $p_coords, $perc) = split /\s/, $pfam_record;
				($s_start, $s_end) = split /\-/, $s_coords;
				($p_start, $p_end) = split /\-/, $p_coords;
				$perc_1 = $perc / 100;

			} else {
				$letter = "H";
				my ($foo, $bar);
				($p_id, $foo, $bar, $s_coords) = split /\s/, $pfam_record;
				($s_start, $s_end) = split /\-/, $s_coords;
				$perc = 80; # LOOOOOOOL

			}

			# Add domain to all_domains lib
			$domains->{$p_id}->{$letter}->{$seq_id} = undef
				unless exists $domains->{$p_id}->{$letter}->{$seq_id};


			# Add transcriptome pfam domain info
			$profiles{$seq_id} = []
				unless exists $profiles{$seq_id};

			if ($c_flag) {
				push @{ $profiles{$seq_id} }, [
					$p_id,
					$s_start,
					$s_end,
					$p_start,
					$p_end,
					$perc_1
				];
			} else {
				push @{ $profiles{$seq_id} }, [
					$p_id,
					$s_start,
					$s_end,
					$perc
				];
			}

		}
	}

	return \%profiles;

}

#--------------------------------------------------------------------------------
sub align_pfams_hgnc {
	my $hgnc_profiles    = shift;
	my $transcriptome_profiles = shift;
	my $domain_info      = shift;

	foreach my $seq_id (keys %{$hgnc_profiles}) {
		my @shared_transcriptomes   = ();
		my @seq_profile       = ();
		my %shared_c_profiles = ();

		# Get HGNC profile and shared transcriptomes
		foreach my $pfam (@{ $hgnc_profiles->{$seq_id} }) {
			my @transcriptomes_with_prof = keys %{ $domain_info->{$pfam->[0]}->{C} };

			push @shared_transcriptomes, @transcriptomes_with_prof;
			push @seq_profile, $pfam->[0];

		}

		next unless @shared_transcriptomes;

		# Get shared transcriptome profiles
		foreach my $transcriptome (@shared_transcriptomes) {

			if (! exists $shared_c_profiles{$transcriptome}) {
				$shared_c_profiles{$transcriptome} = [];

				foreach my $c_pfam (@{ $transcriptome_profiles->{$transcriptome} }) {
					push @{ $shared_c_profiles{$transcriptome} }, [$c_pfam->[0], $c_pfam->[5]];
				}

			}

		} # foreach cth

		make_alignment(
			$seq_id,
			\@seq_profile,
			\%shared_c_profiles
		);

	} # foreach HGNC

	return;
}

#--------------------------------------------------------------------------------
sub align_pfams_transcriptome {
	my $transcriptome_profiles = shift;
	my $hgnc_profiles    = shift;
	my $domain_info      = shift;

	foreach my $seq_id (keys %{ $transcriptome_profiles }) {
		my @shared_hgnc       = ();
		my @seq_profile       = ();
		my %shared_h_profiles = ();

		foreach my $pfam (@{ $transcriptome_profiles->{$seq_id} }) {
			my @hgnc_with_prof = keys %{ $domain_info->{$pfam->[0]}->{H} };

			push @shared_hgnc, @hgnc_with_prof;
			push @seq_profile, [$pfam->[0], $pfam->[5]];
		}

		next unless @shared_hgnc;

		foreach my $hgnc (@shared_hgnc) {

			if (! exists $shared_h_profiles{$hgnc}) {
				$shared_h_profiles{$hgnc} = [];

				foreach my $h_pfam (@{ $hgnc_profiles->{$hgnc} }) {
					push @{ $shared_h_profiles{$hgnc} }, $h_pfam->[0];
				}

			}

		} # foreach hgnc

		make_alignment_transcriptome(
			$seq_id,
			\@seq_profile,
			\%shared_h_profiles
		);


	} # For Cthulhu

	return;

}


#--------------------------------------------------------------------------------
sub make_alignment {
	my $seq_id              = shift;
	my $seq_profile         = shift;
	my $other_seqs_profiles = shift;
	my %results = ();

	foreach my $other_seq ( keys %{$other_seqs_profiles} ) {
		my @other_profile_entries = @{ $other_seqs_profiles->{$other_seq}};

		my @other_profile = ();
		my @profile_perc = ();

		foreach my $prof_entry (@other_profile_entries) {
			push @other_profile, $prof_entry->[0];
			push @profile_perc, $prof_entry->[1];;
		}


		# Change pfam names for numbers
		my %alias_dict = ();
		my $counter    = 0;
		# HGNC...
		my @hgnc_alias_profile = ();

		my $c_flag = 0;
		my $hgnc_alias_profile = assign_alias(
			\$counter,
			$seq_profile,
			\%alias_dict,
		);

		$c_flag = 1;
		my $transcriptome_alias_profile = assign_alias(
			\$counter,
			\@other_profile,
			\%alias_dict,
		);


		my ($score, $identity, $al1, $al2) = needleman_wunsch(
			$hgnc_alias_profile,
			$transcriptome_alias_profile,
			\@profile_perc
		);

		my $hgnc_align    = translate_to_pfams($al1, \%alias_dict);
		my $transcriptome_align = translate_to_pfams($al2, \%alias_dict);


		print "H", "\t", $seq_id, "\t", $other_seq,
		      "\t", $identity, "\t", $score, "\t",
		      join("_",@{ $hgnc_align }), "\t",
		      join("_",@{ $transcriptome_align }), "\n";
	}

	return;

}


sub make_alignment_transcriptome {
	my $seq_id              = shift;
	my $seq_profile         = shift;
	my $other_seqs_profiles = shift;
	my %results   = ();
	my @profile_perc = ();
	my @transcriptome_profile = ();


	# Get transcriptome profile and create perc dictionary
	foreach my $prof_entry (@{$seq_profile}) {
		push @transcriptome_profile, $prof_entry->[0];
		push @profile_perc, $prof_entry->[1];
	}

	foreach my $other_seq (keys %{$other_seqs_profiles}) {
		my @other_profile = @{ $other_seqs_profiles->{$other_seq}};

		my %alias_dict = ();
		$c_flag = 0;
		my $counter = 0;
		my $hgnc_alias_profile = assign_alias(
			\$counter,
			\@other_profile,
			\%alias_dict,
		);

		$c_flag = 1;
		my ($transcriptome_alias_profile, $alias_perc) = assign_alias(
			\$counter,
			\@transcriptome_profile,
			\%alias_dict,
		);


		my ($score, $identity, $al1, $al2) = needleman_wunsch(
			$hgnc_alias_profile,
			$transcriptome_alias_profile,
			\@profile_perc
		);

		my $hgnc_align    = translate_to_pfams($al1, \%alias_dict);
		my $transcriptome_align = translate_to_pfams($al2, \%alias_dict);

		print "C", "\t", $seq_id, "\t", $other_seq,
		      "\t", $identity, "\t", $score, "\t",
		      join("_",@{ $transcriptome_align }), "\t",
		      join("_",@{ $hgnc_align }), "\n";

	} # foreach HGNC

	return;
}

#--------------------------------------------------------------------------------
sub assign_alias {
	my $counter_ref  = shift;
	my $in_profile   = shift;
	my $alias_dict   = shift;
	my @out_profile  = ();


	foreach my $pf ( @{ $in_profile } ) {

		if (not exists $alias_dict->{$pf}) {
			$$counter_ref++;
			$alias_dict->{$pf} = $$counter_ref;
			push @out_profile, $$counter_ref;

		} else {
			push @out_profile, $alias_dict->{$pf};
		}

	} # foreach pfam


	return \@out_profile;

}

#--------------------------------------------------------------------------------
sub translate_to_pfams {
	my $profile    = shift;
	my $dictionary = shift;
	my %reversed = reverse %{$dictionary};

	my @new_profile = map {
		$_ eq '-' ? $_ : $reversed{$_};
	} @{ $profile };

	return \@new_profile;
}

#--------------------------------------------------------------------------------
sub needleman_wunsch {
	my $first_profile  = shift;
	my $second_profile = shift;
	my $profile_perc     = shift;
	my @results        = ();
	my $score          = 0;

	# SETTINGS
	my $MATCH     = 30;
	my $MISSMATCH = -30;
	my $GAP       = -5;


	$results[0]->[0]->{score}   = 0;
	$results[0]->[0]->{pointer} = "none";


	for (my $i = 0; $i < @{ $first_profile }; $i++) {
		$score += $GAP;
		$results[$i+1]->[0]->{score}   = $score;
		$results[$i+1]->[0]->{pointer} = "up";
	}

	my $scnd_score = 0;
	for (my $j = 0; $j < @{ $second_profile }; $j++) {
		$scnd_score += $GAP;
		$results[0]->[$j+1]->{score}   = $scnd_score;
		$results[0]->[$j+1]->{pointer} = "left";
	}


	for (my $i = 0; $i < @{ $first_profile }; $i++) {

		for (my $j = 0; $j < @{ $second_profile }; $j++) {
			# Domains
			my $first_domain  = $first_profile->[$i];
			my $second_domain = $second_profile->[$j];

			# SCORES
			my ($upscore, $leftscore, $diagscore, $match_flag);

			# DIAGONAL SCORE
			if ($first_domain == $second_domain) {
				$diagscore  = $results[$i]->[$j]{score} + $MATCH * $profile_perc->[$j];  # * % transcriptome per domain!!!
				$match_flag = 1;
			} else {
				$diagscore = $results[$i]->[$j]{score} + $MISSMATCH;

			}

			# GAP SCORES
			$upscore   = $results[$i]->[$j+1]->{score} + $GAP;
			$leftscore = $results[$i+1]->[$j]->{score}  + $GAP;

			if ($diagscore >= $upscore) {

				if ($diagscore >= $leftscore) {
					$results[$i+1]->[$j+1]->{score}   = $diagscore;
					$results[$i+1]->[$j+1]->{pointer} = "diagonal";
					$results[$i+1]->[$j+1]->{match}   = 1 if $match_flag;
				} else {
					$results[$i+1]->[$j+1]->{score}   = $leftscore;
					$results[$i+1]->[$j+1]->{pointer} = "left";
				}

			} else {

				if ($upscore >= $leftscore) {
					$results[$i+1]->[$j+1]->{score}   = $upscore;
					$results[$i+1]->[$j+1]->{pointer} = "up";
				} else {
					$results[$i+1]->[$j+1]->{score}   = $leftscore;
					$results[$i+1]->[$j+1]->{pointer} = "left";
				}

			}
		}

	}

	my @align1 = ();
	my @align2 = ();

	# start at last cell of matrix
	my $j = scalar(@{ $second_profile });
	my $i = scalar(@{ $first_profile });

	# STATS
	my $total_score = $results[$i]->[$j]->{score};
	my $identity = 0;
	while (1) {
		last if $results[$i]->[$j]->{pointer} eq "none";

		if ($results[$i]->[$j]->{pointer} eq "diagonal") {
			push @align1, $first_profile->[$i-1];
			push @align2, $second_profile->[$j-1];
			$identity++ if exists $results[$i]->[$j]->{match};
			$i--;
			$j--;
		}
		elsif ($results[$i]->[$j]->{pointer} eq "left") {
			push @align1, "-";
			push @align2, $second_profile->[$j-1];
			$j--;
		} elsif ($results[$i]->[$j]->{pointer} eq "up") {
			push @align1, $first_profile->[$i-1];
			push @align2, "-";
			$i--;
		}
	}

	# Alignment length (num of domains)
	my $al_length = scalar(@align1);


	my @first_align  = reverse(@align1);
	my @second_align = reverse(@align2);

	my $rounded_id = sprintf("%.2f", ($identity / $al_length) * 100);

	return(
		$total_score,
		$rounded_id,
		\@first_align,
		\@second_align
	);

}
