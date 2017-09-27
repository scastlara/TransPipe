#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my %length = ();
my $length_filename = shift @ARGV;
my $pfam_file       = shift @ARGV;

open my $fh, '<', $length_filename
	or die "Can't open $length_filename : $!\n";

while (<$fh>) {
	chomp;
	my ($pfam, $l) = split /\t/, $_;
	next unless $pfam;
	$length{$pfam} = $l;
}

close $fh;

open my $p_fh, "<", $pfam_file
	or die "Can't open $pfam_file : $!\n";


while (<$p_fh>) {
	chomp;
	my ($cthulhu, @pfams) = split /\t/, $_;
	my @coords = ();

	foreach my $pfam (@pfams) {
		my ($id,
			$c_coords,
			$pfam_coords)     = split /\s/, $pfam;
		my ($c_start, $c_end) = split /\-/, $c_coords;
		my ($p_start, $p_end) = split /\-/, $pfam_coords;
		my $length            = $p_end - $p_start + 1;
		my $perc              = ($length / $length{$id}) * 100;
		my $round_p           =  sprintf("%.2f", $perc);

		push @coords, [
			$id,
			$c_start,
			$c_end,
			$round_p,
			$p_start,
			$p_end
		];

	}

	#print "BEFORE:\n";
	#print Dumper(\@coords);
	my $flag = 1;
	my @merged = ();
	if (@coords > 1) {
		merge_coords($flag, \@coords, \@merged);
	} else {
		@merged = @coords;
	}

	#print "AFTER:\n";
	#print Dumper(\@merged);

	print $cthulhu, "\t";

	foreach my $coord (@merged) {
		print $coord->[0], " ",
		      $coord->[1] . "-" . $coord->[2], " ",
		      $coord->[4] . "-" . $coord->[5], " ",
		      $coord->[3];

		print "\t";
	}
	print "\n";

}





sub merge_coords {
	my $flag   = shift;
	my $coords = shift;
	my $merged = shift;

	my @merged_coords = ();
	return unless ($flag);
	$flag = 0;

	if (@{ $coords } == 1) {
		push @merged_coords, [@{$coords->[0]}];
		# print Dumper(\@merged_coords);
		return;
	}

	for (my $i = 0; $i < @{ $coords }; $i++) {
		my $current_id      = $coords->[$i]->[0];
		my $current_c_start = $coords->[$i]->[1];
		my $current_c_end   = $coords->[$i]->[2];
		my $current_perc    = $coords->[$i]->[3];
		my $current_p_start = $coords->[$i]->[4];
		my $current_p_end   = $coords->[$i]->[5];

		if ($i == $#{ $coords }) {
			push @merged_coords, [
				$current_id,
				$current_c_start,
				$current_c_end,
				$current_perc,
				$current_p_start,
				$current_p_end
			];
			next;
		}

		my $next_id      = $coords->[$i + 1]->[0];
		my $next_c_start = $coords->[$i + 1]->[1];
		my $next_c_end   = $coords->[$i + 1]->[2];
		my $next_perc    = $coords->[$i + 1]->[3];
		my $next_p_start = $coords->[$i + 1]->[4];
		my $next_p_end   = $coords->[$i + 1]->[5];

		my $overlap        = $current_p_end - $next_p_start;
		my $overlap_factor = $length{$current_id} * 0.25; # 25 % of pfam domain length

		my $cthulhu_distance = $next_c_start - $current_c_end;
		my $distance_max     = $next_p_start - $current_p_end + $overlap_factor;

		# print STDERR "Df: $current_p_end\tDi: $next_p_start\n",
		#              "Ofactor: $overlap\tOmax: $overlap_factor\n" if ($current_id eq $next_id);

		#print "Overlap: $overlap - Factor: $overlap_factor\n",
		#       "Distance: $cthulhu_distance - Max: $distance_max\n";

		# We will merge domains with a distance (on chtulhu) between them
		# that is equal or less than the "missing" pfam domain length + 25%  (arbitrary).

		if ($current_id eq $next_id  and $overlap < $overlap_factor and $cthulhu_distance <= $distance_max) {
		# same domain with overlap < 25% of pfam length and "short" distance on cthulhu
			my $new_c_start = $current_c_start;
			my $new_c_end   = $next_c_end;

			my ($new_p_start) = sort { # minimum value!
				$a > $b;
			} ($current_p_start, $next_p_start);

			my ($new_p_end) = sort { # max value!
				$a < $b;
			} ($current_p_end, $next_p_end);

			my $new_length = $new_p_end - $new_p_start + 1;
			my $new_perc   = ($new_length / $length{$current_id}) * 100;
			my $rounded_perc = sprintf("%.2f", $new_perc);

			push @merged_coords, [
				$current_id,
				$new_c_start,
				$new_c_end,
				$rounded_perc,
				$new_p_start,
				$new_p_end
			];
			$flag = 1;
			$i++;

		} else {
			push @merged_coords, [
				$current_id,
				$current_c_start,
				$current_c_end,
				$current_perc,
				$current_p_start,
				$current_p_end
			];

		}# if same domain
	}

	@{$merged} = @merged_coords;
	merge_coords($flag, \@merged_coords, $merged);
	return;

}
