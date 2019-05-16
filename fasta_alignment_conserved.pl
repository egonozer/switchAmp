#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
fasta_alignment_conserved.pl

Identifies conserved regions in a fasta-formatted alignment. Outputs
coordinates of non-conserved regions relative to one of the
aligned sequences (i.e. a reference) and sequences of the genomes in each of
the non-conserved regions.

Required:
  -f    fasta-formatted alignment
  -r    ID of the reference sequence in the alignment
  
Optional:
  -o    prefix of output sequences
        (default: 'output')

";

# command line processing
use Getopt::Std;
our ($opt_f, $opt_r, $opt_o);
getopts('f:r:o:');

die $usage unless $opt_f and $opt_r;

my $afile   = $opt_f;
my $refid   = $opt_r;
my $pref    = $opt_o ? $opt_o : "output";

open (my $in, "<$afile") or die "ERROR: Can't open $afile: $!\n";
my @seqs;
my $rseq;
my ($id, $seq);
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            if ($id eq $refid){
                $rseq = $seq;
            } else {
                push @seqs, ([$id, $seq]);
            }
        }
        $seq = "";
        $id = substr($line, 1);
        next;
    }
    $line =~ s/\s//g;
    $line = uc($line);
    $seq .= $line;
}
close ($in);
if ($id){
    if ($id eq $refid){
        $rseq = $seq;
    } else {
        push @seqs, ([$id, $seq]);
    }
}
die "ERROR: reference sequence $refid not found\n" unless $rseq;
die "ERROR: Only one sequence present in alignment file\n" unless @seqs;

my $aleng = length($rseq);
my $rpos = 0;
my ($astart, $astop, $rstart, $rstop);
my @conserved;
for my $i (0 .. $aleng - 1){
    my $rbase = substr($rseq, $i, 1);
    $rpos++ unless $rbase !~ m/[ACGT]/;
    my %bhash;
    $bhash{$rbase} = 1;
    foreach (@seqs){
        my ($qid, $qseq) = @{$_};
        $bhash{substr($qseq, $i, 1)} = 1;
    }
    if (keys %bhash == 1){
        #conserved
        $astart = $i + 1 unless ($astart);
        $rstart = $rpos unless ($rstart);
        ($astop, $rstop) = ($i+1, $rpos);
    } else {
        if ($astart){
            push @conserved, ([$astart, $astop, $rstart, $rstop]);
        }
        ($astart, $astop, $rstart, $rstop) = ("") x 4;
    }
}
if ($astart){
    push @conserved, ([$astart, $astop, $rstart, $rstop]);
}

open (my $sout, ">$pref.non-conserved.seqs.txt");
my ($nc_astart, $nc_rstart) = (0) x 2;
my $seg = 0;
foreach (@conserved){
    ($astart, $astop, $rstart, $rstop) = @{$_};
    my $nc_astop = $astart - 1;
    my $nc_rstop = $rstart;
    my $dist = $nc_astop - $nc_astart;
    if ($dist == 0){
        ($nc_astart, $nc_rstart) = ($astop, $rstop);
        next;
    }
    $seg++;
    print $sout "=$seg,$refid,$nc_rstart,$nc_rstop,$dist\n";
    my $rsub = substr($rseq, $nc_astart, $dist);
    $rsub =~ s/-//g;
    print $sout "$refid\t$rsub\n";
    foreach (@seqs){
        my ($qid, $qseq) = @{$_};
        my $qsub = substr($qseq, $nc_astart, $dist);
        $qsub =~ s/-//g;
        print $sout "$qid\t$qsub\n";
    }
    ($nc_astart, $nc_rstart) = ($astop, $rstop);
}
unless ($nc_astart == $aleng){
    $seg++;
    my $dist = ($aleng - 1) - $nc_astart;
    my $nc_rstop = $aleng + 1;
    print $sout "=$seg,$refid,$nc_rstart,$nc_rstop,$dist\n";
    my $rsub = substr($rseq, $nc_astart, $dist);
    $rsub =~ s/-//g;
    print $sout "$refid\t$rsub\n";
    foreach (@seqs){
        my ($qid, $qseq) = @{$_};
        my $qsub = substr($qseq, $nc_astart, $dist);
        $qsub =~ s/-//g;
        print $sout "$qid\t$qsub\n";
    }
}