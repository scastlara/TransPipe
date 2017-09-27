#!/usr/bin/perl
=pod

=cut



use warnings;
use strict;
use Graphs::Compare;
use Getopt::Long;
use Cwd 'abs_path';
use IO::Zlib;


#===============================================================================
# VARIABLES AND OPTIONS
#===============================================================================


our $PROGRAM       = "TransPipe.pl";
our $VERSION       = 'v0.0.1';
our $USER          = $ENV{ USER };
our $HOST          = $ENV{ HOSTNAME } ? $ENV{ HOSTNAME } : "host";
our $W_DIRECTORY   = $ENV{ PWD };
our $INSTALL_PATH  = get_installpath();

my %OPTS;

# DEFAULTS
$OPTS{blastevalue} = "1e-20";
$OPTS{nogevalue}   = "1e-20";
$OPTS{pfamevalue}  = "1e-20";

# PARSE COMMAND LINE ARGUMENTS
GetOptions (
    \%OPTS            ,
    'all'             ,
    'help|?'          ,
    'transcriptome=s' ,
    'subject=s'       ,
    'nogs=s'          ,
    'pfam=s'          ,
    'blastevalue=s'   ,
    'nogevalue=s'     ,
    'pfamevalue=s'    ,
 );

 # CREATE DIRECTORY STRUCTURE
my %DIRS = (
    "BLAST"      => "BLAST/",
    "PFAM"       => "PFAM/",
    "NOGS"       => "NOGS/",
    "CLUSTERS"   => "CLUSTERS/",
    "NOG_HMM"    => "NOGS/HMM/",
    "ORFS"       => "ORFS/",
    "PLOTS"      => "PLOTS/",
    "LOGS"       => "LOGs/",
    "INFOSEQ"    => "INFOSEQ/",
    "BLAST_DB"   => "BLAST/DB/",
    "GO"         => "GO/",
    "COMPARISON" => "COMPARISON/"
);


#===============================================================================
# STARTING
#===============================================================================

# OPEN LOG FILEHADLE
open my $LOG, ">>", "pipeline.log"
    or die "Can't create logfile : $!\n";

# PRINT START LOG
print_start();
my $done_hash = check_already_done();

# DO YOU WANT TO DO EVERYTHING?
if ($OPTS{all}) {
    $done_hash = undef;
}


#===============================================================================
# FUNCTIONS
#===============================================================================

sub get_installpath {
    my $path = abs_path($0);
    $path =~ s/(.+)\/.*?$/$1\//;
    return($path);
}

# ------------------------------------------------------------------------------
sub print_start {
    my $starttime = localtime();
    print $LOG <<EOF

#------------------------------------------------------------#
#                                                            #
#                       FINDING HOMOLOGS                     #
#                         TransPipe.pl                       #
#------------------------------------------------------------#

                         USER: $USER
                         HOST: $HOST
                      VERSION: $VERSION
            WORKING DIRECTORY: $W_DIRECTORY
           PROGRAM STARTED AT: $starttime

           OPTIONS:
                  transcriptome: $OPTS{transcriptome}
                  subject      : $OPTS{subject}
                  nogs         : $OPTS{nogs}
                  pfam         : $OPTS{pfam}
                  blastevalue  : $OPTS{blastevalue}
                  nogevalue    : $OPTS{nogevalue}
                  pfamevalue   : $OPTS{pfamevalue}
                  all          : $OPTS{all}



EOF
;


    return;
}

# ------------------------------------------------------------------------------
sub print_end {
    my $endtime = localtime();
    print $LOG <<EOF


#------------------------------------------------------------#
#                                                            #
#                       PROGRAM FINISHED                     #
#                         TransPipe.pl                       #
#------------------------------------------------------------#
                         USER: $USER
                         HOST: $HOST
                      VERSION: $VERSION
            WORKING DIRECTORY: $W_DIRECTORY
          PROGRAM FINISHED AT: $endtime
EOF
;
    return;
}

# ------------------------------------------------------------------------------
sub print_job {
    my $string       = shift;
    my $current_time = localtime();
    my $title_l      = length($string);
    my $remaining    = 52 - $title_l;
        # 52 is the number of columns in the title
    my $job_title = " " x int($remaining/2) . $string . " " x int($remaining/2);
    my $header    = "#" . "-" x (int($remaining/2) * 2 + $title_l) . "#";

    my $s_toprint = <<"EOF"
    $header
    #$job_title#
    $header
    #   Job started at: $current_time \
EOF
;

print $LOG $s_toprint, "\n";

}

# ------------------------------------------------------------------------------
sub check_already_done {
    my $logfile   = "pipeline.log";
    my $done_hash = ();

    open my $fh, "<", $logfile
        or do {
            print $LOG "    #    No pipeline.log file. Will do all the jobs.\n";
            return;
        };

    while (<$fh>) {
        chomp;
        next unless /\.\.\. ok/;
        if (/([\w\d_:\-]+)\.\.\./i) {
            $done_hash->{$1} = 1;

        }
    }
    return($done_hash);

}

# ------------------------------------------------------------------------------
sub check_status {
    my $job_string   = shift;
    my $status       = shift;
    my $out_files    = shift;
    my $spaces       = " " x 4;
    my $current_time = localtime();

    my @status_arr = &check_syscall_exit($status);

    if ($status_arr[0] == 1) {
        print $LOG "    #    $job_string... ok\n";
        if ($out_files) {
            print $LOG "    #    Outfiles: ", join("\n$spaces#" . " " x 14 , @{$out_files}), "\n";
        }
        print $LOG "    #    Status: $status_arr[1]\n";
        print $LOG "    #    Job finished at: $current_time\n\n";

    } else {
        print $LOG "    #    $job_string... not ok\n";
        print $LOG "    #    [FATAL ERROR]\n$spaces#${spaces}$status_arr[1]\n\n";
        exit(1);
    }
    return;
}

# ------------------------------------------------------------------------------
sub check_syscall_exit {
    my $prog_exit = 0xffff & shift;
    my ($T, $F) = (1, 0);
    my ($exitflg,$exitstr) = ($F,'');
    $exitstr = sprintf("Command returned %#04x : ", $prog_exit);
    if ($prog_exit == 0) {
        $exitflg = $T;
        $exitstr .= "ran with normal exit ...";
    }
    elsif ($prog_exit == 0xff00) {
        $exitstr .= "command failed: $! ...";
    }
    elsif (($prog_exit & 0xff) == 00) {
        $prog_exit >>= 8;
        $exitstr .= "ran with non-zero exit status $prog_exit ...";
    }
    else {
        $exitstr .= "ran with ";
        if ($prog_exit &   0x80) {
            $prog_exit &= ~0x80;
            $exitstr .= "coredump from ";
            };
        $exitstr .= "signal $prog_exit ...";
    };
    return ($exitflg,$exitstr,$prog_exit);
} # check_syscall_exit

# ------------------------------------------------------------------------------
sub check_program {
    # This checks if a program is installed.
    my $program = shift;
    my $cmd     = `which $program`;
    return $cmd ? 1 : 0;
}

# ------------------------------------------------------------------------------
sub try_load {
    # Checks if a particular module is installed
    my $mod = shift;

    eval("use $mod");
    if ($@) {
        return(0);
    } else {
        return(1);
    }
}

#===============================================================================
# SUB-PROCESSES
#===============================================================================

sub check_all_programs {
    # This checks if all the programs are installed.
    my $job_title = shift;
    print_job($job_title);
    my @programs = qw(
        perl
        zcat
        gawk
        infoseq
        Rscript
        blastx
        tblastn
        makeblastdb
        hmmsearch
        /usr/local/molbio/bin/coverage_blastshorttbl.pl
        /usr/local/molbio/bin/cdna2orfs.pl
        graphcompare
    );

    foreach my $p (@programs) {
        if (check_program($p)) {
            print $LOG "    #    $p installed.\n";
        } else {
            print $LOG "    #    [FATAL ERROR] $p not installed ('which' not working).\n";
            exit(1);
        }
    }
    print $LOG "\n";
    return;
}

# ------------------------------------------------------------------------------
sub check_modules {
    # Check if perl modules are installed
    my $job_title = shift;
    print_job($job_title);
    my @modules = qw(
        global
        List::Util
    );

    foreach my $m (@modules) {
        if (try_load($m)) {
            print $LOG "    #    $m available.\n";
        } else {
            print $LOG "    #    [FATAL ERROR] $m not available ('use module' not working).\n";
            exit(1);
        }
    }

    print $LOG "\n";

    return;
}

# ------------------------------------------------------------------------------
sub create_dirs {
    # Creates the necessary directories.
    my $job_title = shift;

    return if exists $done_hash->{$job_title};
    print_job("CREATING_DIRECTORIES");

    mkdir($DIRS{$_}) && print $LOG "    #    Creating $DIRS{$_} directory\n"
        for sort keys %DIRS;

    print $LOG "    #    $job_title... ok\n";
    return;
}

# ------------------------------------------------------------------------------
sub infoseq {
    # This function runs the program infoseq and creates a file with the output.

    my $job_title = shift;
    my $infile    = shift;
    my $seqtitle  = shift;

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "infoseq $infile -only -name -length \\
        > $DIRS{INFOSEQ}/$seqtitle.infoseq 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$DIRS{INFOSEQ}/$seqtitle.infoseq"]);
    return;
}

# ------------------------------------------------------------------------------
sub plot_seqlen {
    # This function takes two infoseq files and creates a plot with the sequence
    # lengths. It uses the program infoseq_length.R, bundled with TransPipe.

    my $job_title = shift;
    my $outfile   = shift;

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/infoseq_length.R 2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall,[$outfile]);

    return;

}

# ------------------------------------------------------------------------------
sub makeblastdb {
    # This uses the program makeblastdb to create a BLAST database out of a FASTA file.

    my $job_title = shift;
    my $outfile   = shift;
    my $dbtype    = shift;
    my $seqtitle  = shift;

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "makeblastdb  \\
        -in $OPTS{$seqtitle}  \\
        -dbtype \"$dbtype\"  \\
        -out $DIRS{BLAST_DB}$seqtitle.db \\
        >> $DIRS{LOGS}/$seqtitle.makeblastdb.log \\
        2>> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, [$outfile]);

}

# ------------------------------------------------------------------------------
sub get_blast_format {
    # This is the format that will be used for the BLASTS

    my @blast_format = qw (
    6 qseqid
    qlen sseqid
    slen qstart
    qend sstart
    send length
    score evalue
    bitscore pident
    nident mismatch
    positive gapopen
    gaps ppos
    qframe sframe
    qseq sseq
    );
    my $format = join(" ", @blast_format);
    return $format;
}

# ------------------------------------------------------------------------------
sub run_blast {
    # Runs blast!

    my $job_title = shift;
    my $program   = shift;
    my $format    = shift;
    my ($q, $t, $out)     = ("", "", "");

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    if ($program eq "blastx") {
        $q   = "$OPTS{transcriptome}";
        $t   = "$DIRS{BLAST_DB}subject.db";
        $out = "transcriptome-x-subject.blastx";
    } elsif($program eq "tblastn") {
        $q   = "$OPTS{subject}";
        $t   = "$DIRS{BLAST_DB}transcriptome.db";
        $out = "subject-x-transcriptome.tblastn";
    }

    my $code =
        "$program \\
            -query \"$q\" \\
            -db \"$t\" \\
            -num_threads 12 \\
            -evalue $OPTS{blastevalue} \\
            -outfmt '$format' \\
            -out $DIRS{BLAST}/$out.out \\
            2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);

    check_status($job_title, $syscall, ["$out.out"]);
    return;
}

# ------------------------------------------------------------------------------
sub get_best_hits {
    # Uses coverage_blastshorttbl to get BLAST best hits

    my $job_title = shift;
    my $program   = shift;
    my $in  = "";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    if ($program eq "blastx") {
        $in  = "$DIRS{BLAST}/transcriptome-x-subject.blastx.out";

    } elsif($program eq "tblastn") {
        $in = "$DIRS{BLAST}/subject-x-transcriptome.tblastn.out";
    }
    my $out = $in . ".besthits";

    my $code =
        "perl /usr/local/molbio/bin/coverage_blastshorttbl.pl \\
        -prog " . uc($program) . " \\
         \"$in\" \\
         > \"$out\" \\
         2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub get_brh {
    # Reads two blast besthits files and creates a new file with BRH

    my $job_title = shift;
    my $out       = "$DIRS{BLAST}/BLAST_BRHs.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "$INSTALL_PATH/bin/getBRH.pl \\
        $DIRS{BLAST}/transcriptome-x-subject.blastx.out.besthits \\
        $DIRS{BLAST}/subject-x-transcriptome.tblastn.out.besthits \\
        > $out 2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub plot_brh {
    # Function that reads output of getBRH.pl and plots a scatterplot of coverage

    my $job_title = shift;
    my $out = "PLOTS/BLAST_BRH_COVERAGE.pdf";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/BRH_coverage.R \\
            $DIRS{BLAST}/BLAST_BRHs.tbl \\
            $out \\
            BLAST Transcript Subject\\
            2> $DIRS{LOGS}/$job_title.log;";
    my $syscall = system($code);

    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub predict_orfs {
    # Function that predicts longest orfs out of fasta file

    my $job_title = shift;
    my $out = "$DIRS{ORFS}/longestorfs.fa";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "/usr/local/molbio/bin/cdna2orfs.pl -VlBlu \\
        $OPTS{transcriptome} 2> $DIRS{LOGS}/$job_title.log | \\
        perl -ne \'if (/^>/) {s/\\.\\d[+-]\\..*?\$//;print} else {print}\' \\
        > $out";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub plot_longorf {
    # This plots a scatterplot of transcript length vs orflength

    my $job_title = shift;
    my $out = "PLOTS/LONGEST_ORFS.pdf";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/longestorf_length.R 2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub hmmsearch {
    # Aligns ORFs with NOGs or PFAM HMM database

    my $job_title = shift;
    my ($out_tbl, $out_dom, $out_raw, $hmm, $eval);

    if ($job_title eq "HMMSEARCH_NOGS") {
        $out_tbl = "$DIRS{NOGS}/transcript-x-nogs.hmmsearch.tbl";
        $out_dom = "$DIRS{NOGS}/transcript-x-nogs.hmmsearch.domains";
        $out_raw = "$DIRS{NOGS}/transcript-x-nogs.hmmsearch.gz";
        $hmm     = "$OPTS{nogs}";
        $eval    = "$OPTS{nogevalue}";
    } else {
        $out_tbl = "$DIRS{PFAM}/transcript-x-pfam.hmmsearch.tbl";
        $out_dom = "$DIRS{PFAM}/transcript-x-pfam.hmmsearch.domains";
        $out_raw = "$DIRS{PFAM}/transcript-x-pfam.hmmsearch.gz";
        $hmm     = "$OPTS{pfam}.hmm";
        $eval    = "$OPTS{pfamevalue}";
    }

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "bash -c '(set -o pipefail && hmmsearch \\
            -E $eval \\
            --domE $eval \\
            --domtblout $out_dom \\
            --tblout $out_tbl \\
            $hmm $DIRS{ORFS}/longestorfs.fa \\
            2> $DIRS{LOGS}/$job_title.log \\
            | gzip > $out_raw;)'";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out_tbl", "$out_dom", "$out_raw"]);
    return;
}

# ------------------------------------------------------------------------------
sub nogs_coverage {
    # Reads hmmsearch output and computes coverage for each hit

    my $job_title = shift;
    my $out       = "$DIRS{NOGS}/NOGS_coverage.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "perl $INSTALL_PATH/bin/hmmer_coverage.pl \\
        $DIRS{NOGS}/transcript-x-nogs.hmmsearch.domains \\
        > $out 2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub nogs_brh {
    # Gets Nogs Best Reciprocal hits

    my $job_title = shift;
    my $out       = "$DIRS{NOGS}/NOGS_BRHs.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        'gawk "BEGIN{OFS = \"\t\"; print \"TRANSCRIPT\",\"SUBJECT\", \"TRANS_COV\", \"SUBJ_COV\"}
               { if (\$12 == 1) { print \$1, \$2, \$8*100, \$9*100 } }" ' .
        "$DIRS{NOGS}/NOGS_coverage.tbl > $out \\
        2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub plot_nogcoverage {
    # Compares coverage of ORF and HMM in NOG alignment

    my $job_title = shift;
    my $out       = "$DIRS{PLOTS}/NOGS_BRH_COVERAGE.pdf";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/BRH_coverage.R \\
            $DIRS{NOGS}/NOGS_BRHs.tbl \\
            $out \\
            NOGS ORF NOG\\
            2> $DIRS{LOGS}/$job_title.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub pfam_lengths {
    # Reads PFAM hmm database and gets all the lengths of the PFAM domains

    my $job_title = shift;
    my $out = "$DIRS{PFAM}/pfam_lengths.tbl";
    my $in  = "$OPTS{pfam}.hmm.dat.gz";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    open(my $fh, "gzip -dc $in |")
        or ( print $LOG "Can't open to $in $!\n" &&
             check_status($job_title, 1, ["$out"]) );

    open my $ofh, ">", $out
        or ( print $LOG "Can't write to $out $!\n" &&
             check_status($job_title, 1, ["$out"]) );

    local $/ = "//";
    while (<$fh>) {
        chomp;
        my ($domain, $length);
        if (/#=GF\s+AC\s+(.+)\n/g) {
            $domain = $1;
        }
        if (m/#=GF\s+ML\s+(.+)\n/g) {
            $length = $1;
        }
        next unless ($domain && $length);
        print $ofh "$domain", "\t", "$length\n";
    }
    close $fh;
    close $ofh;

    my $syscall = 0; # If made it to here, everything is fine
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub subject_pfamtable {
    # Gets all the PFAM domains of the human subject proteins

    my $job_title = shift;
    my $out       = "$DIRS{PFAM}/subject_pfam.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "$INSTALL_PATH/bin/subject_pfam_table.pl \\
              $OPTS{subject} $OPTS{pfam}.full.gz \\
              > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub transcriptome_pfamtable {
    # Creates table with pfam domains of transcriptome sequences

    my $job_title = shift;
    my $out       = "$DIRS{PFAM}/transcriptome_pfam.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
    "$INSTALL_PATH/bin/transcriptome_pfam_table.pl \\
          $DIRS{PFAM}/transcript-x-pfam.hmmsearch.domains \\
          > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub pfam_overlap {
    # Fixes pfam domain overlap in query sequences

    my $job_title = shift;
    my $out       = "$DIRS{PFAM}/transcriptome_pfam.fix.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "$INSTALL_PATH/bin/pfam_overlap.pl \\
              $DIRS{PFAM}/pfam_lengths.tbl $DIRS{PFAM}/transcriptome_pfam.tbl \\
              > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub meta_align {
    # Meta-alignment of two sequences using pfam-domains

    my $job_title = shift;
    my $out       =  "$DIRS{PFAM}pfam_alignment.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "$INSTALL_PATH/bin/AlignPFAM.pl \\
            $DIRS{PFAM}/pfam_lengths.tbl \\
            $DIRS{PFAM}/transcriptome_pfam.fix.tbl \\
            $DIRS{PFAM}/subject_pfam.tbl \\
            > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub filter_pfamscore {
    # Removes pfam alignments with score < 0
    my $job_title = shift;
    my $out       = "$DIRS{PFAM}pfam_alignment.filtered.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        'gawk "{if (\$5 > 0) {print \$0}}" ' .
        "$DIRS{PFAM}pfam_alignment.tbl \\
        > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub pfam_besthits {
    # Gets meta-alignment best hits

    my $job_title = shift;
    my $out       =  "$DIRS{PFAM}pfam_alignment.filtered.besthits.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "$INSTALL_PATH/bin/besthits_pfam.pl \\
        $DIRS{PFAM}pfam_alignment.filtered.tbl \\
        > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub pfam_brh {
    # Gets meta-alignment BRHs

    my $job_title = shift;
    my $out       = "$DIRS{PFAM}/PFAM_BRHs.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "$INSTALL_PATH/bin/getBRH_pfam.pl \\
         $DIRS{PFAM}pfam_alignment.filtered.besthits.tbl \\
         > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub plot_pfamscore {
    # Plots pfam scores of meta-alignment

    my $job_title = shift;
    my $out       = "$DIRS{PLOTS}/PFAM_SCORES.pdf";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/pfam_score.R \\
         $DIRS{PFAM}/PFAM_BRHs.tbl $out \\
         2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub count_brh {
    # Counts BRHs from NOG, BLAST and PFAM (transcript and subject)

    my $job_title = shift;
    my @out = (
        "BLAST_transcripts.tbl",
        "NOGS_transcripts.tbl",
        "PFAM_transcripts.tbl",
        "BLAST_subject.tbl",
        "NOGS_subject.tbl",
        "PFAM_subject.tbl"
    );

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "gawk -f $INSTALL_PATH/bin/get_transcripts.awk \\
         outdir='COMPARISON/' type='transcripts' col=1 $DIRS{BLAST}BLAST_BRHs.tbl $DIRS{NOGS}NOGS_BRHs.tbl $DIRS{PFAM}PFAM_BRHs.tbl; \\
         gawk -f $INSTALL_PATH/bin/get_transcripts.awk \\
          outdir='COMPARISON/' type='subject' col=2 $DIRS{BLAST}BLAST_BRHs.tbl $DIRS{NOGS}NOGS_BRHs.tbl $DIRS{PFAM}PFAM_BRHs.tbl;";

    my $syscall = system($code);
    check_status($job_title, $syscall, \@out);
    return;
}

# ------------------------------------------------------------------------------
sub plot_brhcount {
    # Compares BRH

    my $job_title = shift;
    my @out = ("TRANSCRIPTOME_BRH_VENN.svg", "$DIRS{COMPARISON}transcripts_brhs.tbl");

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "graphcompare -in $DIRS{COMPARISON}BLAST_transcripts.tbl,$DIRS{COMPARISON}NOGS_transcripts.tbl,$DIRS{COMPARISON}PFAM_transcripts.tbl \\
                      -v $DIRS{PLOTS}TRANSCRIPTOME_BRH_VENN.svg \\
                      -n $DIRS{COMPARISON}transcriptome_brhs.tbl \\
                      > /dev/null 2> $DIRS{LOGS}/$job_title.transcripts.log; \\
         graphcompare -in $DIRS{COMPARISON}BLAST_subject.tbl,$DIRS{COMPARISON}NOGS_subject.tbl,$DIRS{COMPARISON}PFAM_subject.tbl \\
                      -v $DIRS{PLOTS}SUBJECT_BRH_VENN.svg \\
                      -n $DIRS{COMPARISON}subject_brhs.tbl \\
                      > /dev/null 2>> $DIRS{LOGS}/$job_title.subject.log;";

    my $syscall = system($code);
    check_status($job_title, $syscall, \@out);
    return;
}

# ------------------------------------------------------------------------------
sub join_info {
    # Joins all the homology information

    my $job_title = shift;
    my $out       = "homology_infojoiner.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "perl $INSTALL_PATH/bin/info_joiner.pl \\
            -b $DIRS{BLAST}/transcriptome-x-subject.blastx.out.besthits \\
            -pf $DIRS{PFAM}/pfam_alignment.filtered.besthits.tbl \\
            -e $DIRS{NOGS}NOGS_coverage.tbl \\
            -dict $INSTALL_PATH/data/meNOG.members.ensp.hgnc \\
            -rh_blast $DIRS{BLAST}/BLAST_BRHs.tbl -rh_pfam PFAM/PFAM_BRHs.tbl \\
            -rh_nogs $DIRS{NOGS}/NOGS_BRHs.tbl \\
            -Pfam_ints $INSTALL_PATH//data/3did.filtered.tbl \\
            -domains $DIRS{PFAM}/transcriptome_pfam.fix.tbl \\
            -go_annot $INSTALL_PATH/data/HGNC_2_GO.txt \\
            -go_parents $INSTALL_PATH/data/ \\
            -train_dict $INSTALL_PATH/data/dmel_tr_2_gn.tbl \\
            -plen $INSTALL_PATH/data/Pnet.short.tbl.gz \\
            -hgnc_2_symbol $INSTALL_PATH/data/ALIAS_HUGO_IDS.tbl.2 \\
            > $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub print_homologs {
    # Creates a file with all the homology relationships
    my $job_title = shift;
    my $out       = "homologs.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);
    print "HOLA";
    my $code =
    'gawk \'FNR > 1{print $1"\t"$3"\n"$2"\t"$4}\' ' .
    "homology_infojoiner.tbl | sort | uniq > $out 2> $DIRS{LOGS}/$job_title.log";
    my $syscall = system($code);
    check_status($job_title, $syscall, [$out])
}

# ------------------------------------------------------------------------------
sub predict_interactions {
    # Uses random forest to predict interactions

    my $job_title = shift;
    my $out       = "RESULTS.tbl";

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/rf_classifier.R \\
        $INSTALL_PATH/data/randomforest.downsampled.RData \\
        homology_infojoiner.tbl \\
        $out 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, ["$out"]);
    return;
}

# ------------------------------------------------------------------------------
sub plot_ppairs {
    my $job_title = shift;
    my @out = qw(
        PLOTS/PREDICTION_MOLFUN_NTO.png
        PLOTS/PREDICTION_BIOPROC_NTO.png
        PLOTS/PREDICTION_CELLCOM_NTO.png
        PLOTS/PREDICTION_MOLFUN_NTO.png
        PLOTS/PREDICTION_DOMAIN_INT_SCORE.png
        PLOTS/PREDICTION_PATH_LENGTH.png
    );

    return if exists $done_hash->{$job_title};
    print_job($job_title);

    my $code =
        "Rscript $INSTALL_PATH/bin/prediction_plots.R \\
         RESULTS.tbl 2> $DIRS{LOGS}/$job_title.log";

    my $syscall = system($code);
    check_status($job_title, $syscall, \@out);
    return;
}


# ------------------------------------------------------------------------------
sub save_predictions {
    my $job_title = shift;
    my $out       = "INTERACTOME.tbl";

    open my $fh, "<", "RESULTS.tbl"
        or die "Can't open RESULTS.tbl : $!\n";

    open my $outfh, ">", $out
        or die "Can't write to $out : $!\n";

    my $header = <$fh>; # Skip header
    print $outfh "$header";

    while (<$fh>) {
        chomp;
        my @cols = split /\t/;
        if ($cols[-1] >= 0.6) {
            print $outfh "$_\n";
        }
    }

    close $fh;
    close $outfh;

    my $syscall = 0; # If made it to here, everything is fine
    check_status($job_title, $syscall, ["$out"]);
    return;
}



#===============================================================================
# PIPELINE
#===============================================================================

# Starting (These are always executed)
check_all_programs("CHECKING_PROGRAMS");
check_modules("CHECKING_MODULES");
create_dirs("CREATING_DIRECTORIES");

# Sequence length
infoseq("INFOSEQ_TRANSCRIPTOME", "$OPTS{transcriptome}", "transcriptome");
infoseq("INFOSEQ_SUBJECT", "$OPTS{subject}", "subject");
plot_seqlen("PLOT:SEQUENCE_LENGTH", "PLOTS/SEQUENCE_LENGTH.svg");

# BLAST alignment
makeblastdb("MAKEBLASTDB_SUBJECT", "$DIRS{BLAST_DB}subject.db", "prot", "subject");
makeblastdb("MAKEBLASTDB_TRANSCRIPTOME", "$DIRS{BLAST_DB}transcriptome.db", "nucl", "transcriptome");
my $blast_format = get_blast_format();
run_blast("BLASTX_TRANS-x-SUBJ",  "blastx",  $blast_format);
run_blast("TBLASTN_SUBJ-x-TRANS", "tblastn", $blast_format);
get_best_hits("BLASTX_BESTHITS", "blastx");
get_best_hits("TBLASTN_BESTHITS", "tblastn");
get_brh("BLAST_BRHs");
plot_brh("PLOT:BLAST_BRHs");

# ORF and protein prediction
predict_orfs("LONGEST_ORFS");
infoseq("INFOSEQ_LONGESTORF", "ORFS/longestorfs.fa", "longestorf");
plot_longorf("PLOT:ORFS");

# NOGs alignment
hmmsearch("HMMSEARCH_NOGS");
nogs_coverage("NOGS_COVERAGE");
nogs_brh("NOGS_BRHs");
plot_nogcoverage("PLOT:NOGS_BRHs");

# PFAM meta-alignment
hmmsearch("HMMSEARCH_PFAM");
pfam_lengths("PFAM_LENGTHS");
subject_pfamtable("SUBJECT_PFAM_TABLE");
transcriptome_pfamtable("TRANSCRIPTOME_PFAM_TABLE");
pfam_overlap("FIX_OVERLAP");
meta_align("PFAM_ALIGNMENT");
filter_pfamscore("FILTER_PFAM_SCORES");
pfam_besthits("PFAM_BESTHITS");
pfam_brh("PFAM_BRHs");
plot_pfamscore("PLOT:PFAM_SCORE");

# BRH comparison
count_brh("COUNT_BRHs");
plot_brhcount("PLOT:BRH_count");

# Join information and prepare for classification
join_info("JOINING_INFORMATION");
print_homologs("PRINTING_HOMOLOGS");

# Predict interactions
predict_interactions("PREDICTING_INTERACTIONS");
save_predictions("SAVE_PREDICTIONS");
plot_ppairs("PLOT_PREDICTIONS");

# Ending
print_end();
