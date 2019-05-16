#!/usr/bin/perl

my $license = "
    SwitchAmp
    Copyright (C) 2019 Egon A. Ozer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see [http://www.gnu.org/licenses/].
";

use strict;
use warnings;
use Text::Levenshtein::XS qw/distance/;

$|++;

my $version = 0.2.1;

## Added calculation of number of unique sequences that have only a set number of differences from the parental/wt
## Added ability to count unique sequences with one or more particular variant types (i.e. variants incompatible with growth)
## Calculate total number of reads in each output type including mean, median, min, and max

my $usage = "
switchAmp.pl

Required:
  -q    Input sequences in fastq format (can be gzipped). Multiple files can
        be separated by commas
  -n    file of non-conserved coordinates and sequences
        (as output by fasta_alignment_conserved.pl)

Options:
  -s    fasta file containing full parental sequence. This should be first
        sequence in the file.
        If there is sequence upstream or downstream of the parental sequence
        expected to also be found in the amplicon sequences, then they should
        be included in the file after the parental sequence. The upstream and
        downstream sequences must be given the IDs '>up' and '>dn' respectively
        If no upstream or downstream sequences are included, it will be assumed
        that there is no upstream and/or downstream sequence in the reads
        (default: The parental sequence of GC strain 1-81-S2 and default
        upstream and downstream sequences will be used. Use the -S option to
        print these sequences)
  -S    Print the default parental, upstream, and downstream sequences and quit

  -r    minimum number of reads in each output group
        (default: 3)
  -t    maximum number of differences (mismatches / insertions / deletions)
        in upstream and downstream sequences
        (default: 2)
  -d    do not ignore variants in conserved regions
        (default: variants in conserved regions will be ignored)
  -p    maximum number of differences relative to the parental sequence to be
        counted as parental. This will only apply if none of the variants are
        perfect matches for one or more of the non-conserved sequences. Reads
        fulfilling this criterion will be grouped with the parental sequence.
        (default: 1)
        
  -x    maximum variant sequence difference.
        If a read contains one or more variable regions that differ by more
        than this number of bases from either the parental sequence or any of
        the non-conserved / silent sequences, the read will be omitted
        (default: 2)
        
  -v    file with list of variant types associated with a phenotype of interest.
            Format:
                variant_number1:id1,id2,id3
                variant_number2:id1,id3
                etc.
            'variant number' should match the numbering in the file given to -n.
            'id' should match sequence IDs in the file given to -n. 

  -o    output file prefix
        (default 'output')        
  -h    linkage method. Choices are 'single' or 'complete'
        (default: 'complete')

";

use Getopt::Std;
our ($opt_f, $opt_q, $opt_k, $opt_m, $opt_b, $opt_o, $opt_r, $opt_t, $opt_h, $opt_d, $opt_n, $opt_p, $opt_v, $opt_x, $opt_s, $opt_S);
getopts('f:q:k:m:b:o:r:h:t:dm:n:p:v:x:s:S');
die $usage unless ($opt_f or $opt_q) and $opt_n;

my $keyfile         = $opt_k;
my $ncfile          = $opt_n;
my $max_u_dist      = defined $opt_m ? $opt_m : 2;
my $bintarget       = $opt_b if $opt_b;
my $pref            = $opt_o ? $opt_o : "output";
my $min_rep         = $opt_r ? $opt_r : 3;
my $cl_type         = $opt_h ? $opt_h : "complete";
my $max_flank_dist  = defined $opt_t ? $opt_t : 2;
my $max_wt_diff     = defined $opt_p ? $opt_p : 1;
my $varfile         = $opt_v if $opt_v;
my $max_var_dist    = defined $opt_x ? $opt_x : 2;
my $seqfile         = $opt_s if $opt_s;

my ($wt, $us, $ds);
## Hard-coded parental, upstream, and downstream sequences
## These can be substituted by the user
$wt =
"atgaatacccttcaaaaaggctttacccttatcgagctgatgattgtgatcgctatcgtcggcattttggcgg
cagtcgcccttcccgcctaccaagactacaccgcccgcgcgcaagtttccgaagccatccttttggccgaagg
tcaaaaatcagccgtcaccgagtattacctgaatcacggcatatggccgaaagacaacacttctgccggcgtg
gcatccgcttcaacaatcaaaggcaaatatgttcagaaagttgaagtcgcaaaaggcgtcgttaccgcccaaa
tggcttcaaccggcgtaaacaaagaaatccaagacaaaaaactctccctgtgggccaagcgtcaagacggttc
ggtaaaatggttctgcggacagccggttacgcgcaccggcgacaacgacgacaccgttgccgacgccaacaac
gccatcgacaccaagcacctgccgtcaacctgccgcgatgaatcatctgccacctaaggcaaattaggcctta
aa";
$us =
"TTTCCCCTTTCAATTAGGAGTAATTTT";
$ds =
"TTTTAAATAAATCAAGCGGTAAGTGATTTCCCACGGCCGCCCGGATCAACCCGGGCGGCTTGTCTTTTAAGGG
TTTGCAAGGCGGGCGGGGTCGTCCGTTCCGGTGGAAATAATATATCGATTGCGCTTCAAGGCCCTGCATGTAC
CTCATTGCCACCCGTTTAAACACGGTTTTTATCTGACAGGCGCGCAATCCGCCCCCTCATTTGTTAATCCGCC
ATATTGTATTGAAACACCGCCCGGAACCC";

if ($opt_S){
    print "\nparental:\n$wt\n";
    print "upstream:\n$us\n";
    print "downstream:\n$ds\n\n";
    die ("Exiting\n");
}

if ($seqfile){
    ($wt, $us, $ds) = ("") x 3;
    my ($id, $seq);
    open (my $in, "<$seqfile") or die "ERROR: Can't open $seqfile: $!\n";
    while (my $line = <$in>){
        chomp $line;
        $line =~ s/\s*$//;
        if ($line =~ m/^>/){
            if ($id){
                if (!$wt){
                    $wt = $seq;
                } else {
                    if ($id eq "up"){
                        $us = $seq;
                    } elsif ($id eq "dn"){
                        $ds = $seq;
                    } else {
                        die "ERROR: in file $seqfile, the sequence ID $id is not recognized. Should be 'up' or 'dn' only\n";
                    }
                }
            }
            $id = substr($line, 1);
            $seq = "";
            next;
        }
        $line =~ s/\s//g;
        $seq .= $line;
    }
    close ($in);
    if ($id){
        if (!$wt){
            $wt = $seq;
        } else {
            if ($id eq "up"){
                $us = $seq;
            } elsif ($id eq "dn"){
                $ds = $seq;
            } else {
                die "ERROR: in file $seqfile, the sequence ID $id is not recognized. Should be 'up' or 'dn' only\n";
            }
        }
    }
    die "ERROR: No parental sequence found in file $seqfile\n" unless $wt;
}

$wt =~ s/\s//g;
$wt = uc($wt);
my ($lengus, $lengds) = (0) x 2;
if ($us){
    $us =~ s/\s//g;
    $us = uc($us);
    $lengus = length($us);
}
if ($ds){
    $ds =~ s/\s//g;
    $ds = uc($ds);
    $lengds = length($ds);
}

open (my $statout, ">$pref.stats.txt");
my $refseqs = $seqfile ? $seqfile : "default";
print $statout "version: $version
input sequence: $opt_q
non-conserved sequence file: $opt_n
parental/us/ds sequence file: $refseqs
maximum number of variants in upstream or downstream flanking sequences: $max_flank_dist
minimum representatives: $min_rep
maximum variant distance: $max_var_dist
hierarchical clustering linkage method: $cl_type
";
print $statout "ignore variants in conserved regions: ", $opt_d ? "no\n" : "yes\n";
print $statout "variants of interest file: $varfile\n" if $varfile;
print $statout "\n";

my @strings;

## first, we'll read in the coordinates of the non-conserved regions
open (my $nin, "<$ncfile") or die "ERROR: Can't open $ncfile: $!\n";
my @nc_array;
my @nc_seqs;
my @tmp_nc_seqs;
my $wtid;
my $lstop = 1;

my @seg_names;
my @var_seg_names;
my @varseqs_wt;
my $wt_vseq_stg;
#;
while (my $line = <$nin>){
    chomp $line;
    if ($line =~ m/^=/){
        if (@tmp_nc_seqs){
            push @nc_seqs, [@tmp_nc_seqs];
        }
        @tmp_nc_seqs = ();
        my $tmp = substr($line, 1);
        my ($num, $wt, $start, $stop, $leng, $partial) = split(",", $tmp);
        
        if ($start >= $lstop){
            push @seg_names, "cons:$lstop-$start";
            #print "\tcons:$lstop-$start";
        }
        push @seg_names, "var$num:$start-$stop";
        push @var_seg_names, "var$num:$start-$stop";
        #print "\tvar$num";
        $lstop = $stop;
        
        $wtid = $wt;
        push @nc_array, ([$num, $start, $stop, $leng, $partial]);
        next;
    }
    my ($id, $seq) = split("\t", $line);
    push @tmp_nc_seqs, ([$id, $seq]);
    if ($id eq $wtid){
        @{$varseqs_wt[$#var_seg_names]} = ($id, $seq);
        $wt_vseq_stg .= uc($seq);
    }
}
close ($nin);
#print "\n";
if (@tmp_nc_seqs){
    push @nc_seqs, [@tmp_nc_seqs];
}
undef @tmp_nc_seqs;

## next, if there is a variant-of-interest file, read it in
my %vois; # "vars of interest"
if ($varfile){
    open (my $vin, "<$varfile") or die "ERROR: Can't open $varfile: $!\n";
    while (my $line = <$vin>){
        chomp $line;
        next if $line =~ m/^\s*$/;
        my ($id, $varstring) = split(":", $line);
        $id =~ s/\s//g;
        my @vars = split(",", $varstring);
        foreach my $var (@vars){
            $var =~ s/\s//g;
            $vois{$id}{$var} = 1;
        }
    }
}

## now, read in the amplicon sequences
my $read_in_count = 0;
if ($opt_q){
    @strings = ();
    my @files = split(",", $opt_q);
    foreach my $file (@files){
        my $in;
        if ($file =~ m/\.gz$/){
            open ($in, "gzip -cd $file | ");
        } else {
            open ($in, "<$file") or die "ERROR: Can't open $file: $!\n";
        }
        print STDERR "\rRead in $read_in_count sequences";
        while (<$in>){
            chomp(my $id = $_);
            chomp(my $seq = <$in>);
            chomp(my $id2 = <$in>);
            chomp(my $qual = <$in>);
    
            ## ccs debugging
            push @strings, ([uc($seq), $id]);
            
            $read_in_count++;
            print STDERR "\rRead in $read_in_count sequences" if $read_in_count % 1000 == 0;
        }
        close $in;
    }
    print STDERR "\rRead in $read_in_count total sequences\n";
    print $statout "Read in $read_in_count total sequences\n";
}


#first, group identical sequences
my %hash1;
my $stg_proc_count = 0;
my $filtered = 0;
my $filtered_no_us = 0;
my $filtered_no_ds = 0;
my $filtered_no_us_or_ds = 0;
my ($var_us, $var_ds) = (0) x 2;
print STDERR "\rProcessed $stg_proc_count sequences";

#foreach my $stg (@strings){
#    $stg = uc($stg);
## ccs debugging
foreach my $slice (@strings){
    my ($stg, $id) = @{$slice};
    $stg_proc_count++;
    print STDERR "\rProcessed $stg_proc_count sequences";
    
    if ($us or $ds){
        my ($no_us, $no_ds) = (1) x 2;
        $no_us = 0 unless $us;
        $no_ds = 0 unless $ds;
        my ($close_us, $close_ds) = (0) x 2;
        my $us_leng = length($us) if $us;
        my $ds_leng = length($ds) if $ds;
        my $rcd;
        if ($us){
            # going to assume for now that upstream/downstream sequences will be the same length (no indels or balanced indels)
            # this is obviously not likely to be 100% true, but hopefully will capture some of the errors
            # could potentially use needleman-wunsch to align us and ds sequences for more sensitivity
            my $stg_us = substr($stg, 0, $us_leng);
            if ($stg_us eq $us){
                $no_us = 0;
            } elsif ($max_flank_dist > 0){
                my $dist = distance($stg_us, $us);
                if ($dist <= $max_flank_dist){
                    $no_us = 0;
                    $var_us++;
                }
            }
            if ($no_us){
                ## try reverse complement
                $stg = reverse($stg);
                $stg =~ tr/ACGT/TGCA/;
                $rcd = 1;
                $stg_us = substr($stg, 0, $us_leng);
                if ($stg_us eq $us){
                    $no_us = 0;
                } elsif ($max_flank_dist > 0){
                    my $dist = distance($stg_us, $us);
                    if ($dist <= $max_flank_dist){
                        $no_us = 0;
                        $var_us++;
                    }
                }
            }
            $filtered_no_us++ if $no_us == 1;
            $stg = substr($stg, $us_leng) unless $no_us; #trim the upstream sequence
        }
        if ($ds){
            my $stg_ds = substr($stg, -$ds_leng, $ds_leng);
            if ($stg_ds eq $ds){
                $no_ds = 0;
            } elsif ($max_flank_dist > 0){
                my $dist = distance($stg_ds, $ds);
                if ($dist <= $max_flank_dist){
                    $no_ds = 0;
                    $var_ds++;
                }
            }
            if ($no_ds){
                $stg = reverse($stg);
                $stg =~ tr/ACGT/TGCA/;
                $stg_ds = substr($stg, -$ds_leng, $ds_leng);
                if ($stg_ds eq $ds){
                    $no_ds = 0;
                } elsif ($max_flank_dist > 0){
                    my $dist = distance($stg_ds, $ds);
                    if ($dist <= $max_flank_dist){
                        $no_ds = 0;
                        $var_ds++;
                    }
                }
            }
            $filtered_no_ds++ if $no_ds == 1;
            $stg = substr($stg, 0, -$ds_leng) unless $no_ds; #trim the downstream sequence
        }
        if ($no_us or $no_ds){
            $filtered_no_us_or_ds++ if $no_us and $no_ds;
            next;
        }
        ##ccs debugging
        $hash1{$stg}{'count'}++;
        $hash1{$stg}{'id'} = $id;
    } else {
        if (exists $hash1{$stg}){
            
            #$hash1{$stg}++;
            ##ccs debugging
            $hash1{$stg}{'count'}++;
            $hash1{$stg}{'id'} = $id;
            
            
            #print STDERR "stg:", substr($stg, 0, 5), "\n";
        } else {
            my $stgrc = reverse($stg);
            $stgrc =~ tr/ACGT/TGCA/;
            if (exists $hash1{$stgrc}) {
                
                #$hash1{$stgrc}++;
                ##ccs debugging
                $hash1{$stgrc}{'count'}++;
                $hash1{$stgrc}{'id'} = $id;
                
                #print STDERR "stgrc:", substr($stgrc, 0, 5), "\n";
            } else {
                ## check which orientation is closest to the wild-type sequence
                #my ($for_ls, $rev_ls) = (levenshtein($stg, $wt), levenshtein($stgrc, $wt));
                my ($for_ls, $rev_ls) = (distance($stg, $wt), distance($stgrc, $wt)); 
                if ($for_ls < $rev_ls){
                    
                    #$hash1{$stg}++;
                    ##ccs debugging
                    $hash1{$stg}{'count'}++;
                    $hash1{$stg}{'id'} = $id;
                    
                    #print STDERR "NEWstg:", substr($stg, 0, 5), "\n";
                } else {
                    
                    #$hash1{$stgrc}++;
                    ##ccs debugging
                    $hash1{$stgrc}{'count'}++;
                    $hash1{$stgrc}{'id'} = $id;
                    
                    #rint STDERR "NEWstgrc:", substr($stgrc, 0, 5), "\n";
                }
            }
        }
    }
}
print STDERR "\rProcessed $stg_proc_count sequences\n";
if ($us){
    statprint("\tFound $var_us sequences with imperfect upstream matches\n");
    statprint("\tFiltered $filtered_no_us sequences lacking upstream sequences only\n");
}
if ($ds){
    statprint("\tFound $var_ds sequences with imperfect downstream matches\n");
    statprint("\tFiltered $filtered_no_ds sequences lacking downstream sequences only\n");
}
if ($us and $ds){
    statprint ("\tFiltered $filtered_no_us_or_ds sequences lacking upstream and downstream sequences\n");
}

@strings = ();
my @str_counts;
my @tmp_names;
my @act_names;
my $tmpcount = 0;
my $maxdist = 0; #maximum levenshtein distance is the length of the longest string
my $kout;
if ($opt_d){
    open ($kout, ">$pref.sequence_key.txt");
    print $kout "id\tlength\tnum\tseq\n";
}

#foreach my $key (sort {$hash1{$b} <=> $hash1{$a}} keys %hash1){
#    my $rep = $hash1{$key};
##ccs debugging
foreach my $key (sort {$hash1{$b}{'count'} <=> $hash1{$a}{'count'}} keys %hash1){
    my $rep = $hash1{$key}{'count'};
    
    #next if $rep < $min_rep;
    push @strings, ($key);
    push @str_counts, ($rep);
    $tmpcount++;
    my $tmp_name = "seq" . sprintf("%05d", $tmpcount);
    push @tmp_names, $tmp_name;
    
    ##ccs debugging
    push @act_names, $hash1{$key}{'id'};
    
    my $key_leng = length($key);
    $maxdist = $key_leng if $key_leng > $maxdist;
    print $kout "$tmp_name\t$key_leng\t$rep\t$key\n" if $opt_d;
    #print STDERR "$tmp_name (#$rep) [$key_leng bp]\n";
}
close ($kout) if $opt_d;
print STDERR "\n";
my ($avg, $min, $max, $med, $num, $sum) = stats(\@str_counts);
statprint("unique sequences: ", scalar @strings, "\n");
statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
if ($min_rep > 1){
    ($avg, $min, $max, $med, $num, $sum) = stats(\@str_counts, $min_rep);
    statprint("unique sequences with at least $min_rep reads: $num\n");
    statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
}
print STDERR "\n";


## align the reads to parental sequence
my $rout;
if ($opt_d){
    open ($rout, ">$pref.results.tsv");
    print $rout "id\tcount\t";
    print $rout "has_var\t" if $varfile;
    print $rout join("\t", @seg_names), "\n";
}

my %seg_seqs;
my @varseqs;
my @vartypes;
my @has_voi_j;
my @wt_adjacent_unique;
my $wt_adjacent_total = 0;
my %wt_adjacent_j;
my %exceeds_max_var_d;
my %top_vars_j;
foreach my $j (0 .. $#strings){
    print STDERR "\rAnalyzing ", $j + 1 if ($j+1)%10 == 0;
    my $amp = $strings[$j];
    my $rep = $str_counts[$j];
    my $name = $tmp_names[$j];
    
    ##ccs debugging
    my $actname = $act_names[$j];
    print $rout "$name\t$rep" if $rep >= $min_rep and $opt_d;
    my $toprinttorout;
    my ($aln_w, $aln_a);
    nw_c(uc($wt), uc($amp), $aln_w, $aln_a); ## C implementation of Needleman-Wunsch.
    #($aln_w, $aln_a) = nw($wt, $amp); ## Pure Perl implemnation of Needleman-Wunsch. Sloooooow.

    ## Let's run through the alignment and check out the non-conserved regions
    my $aleng = length($aln_w);
    my @sub_nc = @nc_array; #num, start, stop, leng
    my @sub_sq = @nc_seqs;
    my ($unum, $ustart, $ustop, $uleng, $upart) = @{shift @sub_nc};
    my $last_ustop = $ustop;
    my @useqs = @{shift @sub_sq};
    my ($subseq_w, $subseq_a);
    my $pos_w = 0;
    my $seq_i_start = 0;
    my $ustarted = 0;
    my $seg_num = 0;
    my $var_seg_num = 0;
    my $dist_from_wt = 0;
    my $has_voi;
    my $exceeds_max_var_dist;
    my %top_vars;
    for my $i (0 .. $aleng - 1){
        
        my $base_w = substr($aln_w, $i, 1);
        unless ($base_w eq "-"){
            $pos_w++;
        }
        
        if ($pos_w == $ustop){
            next if $base_w eq "-";
            $ustarted = 0;
            my $dist = ($i - $seq_i_start); # don't need to add 1 since we're one i beyond
            $subseq_w = substr($aln_w, $seq_i_start, $dist);
            $subseq_a = substr($aln_a, $seq_i_start, $dist);
            ## screen the sequences in the non-conserved region against the variant sequences
            $upart = 0 unless $upart;
            my ($degapped, $type, $list) = nc_screen($subseq_a, $upart, \@useqs);
            $varseqs[$j][$var_seg_num] = $degapped;
            my @larray = @{$list};
            my $outtype;
            if ($type == 0 and $larray[0] eq $wtid){
                $toprinttorout .= "\t" if $rep >= $min_rep and $opt_d;
                #print $rout "\t" if $rep >= $min_rep and $opt_d;
            } else {
                $outtype = "$type:" . join(",", @larray);
                my $sub_has_voi;
                my $is_wt;
                foreach (@larray){
                    $top_vars{$_}++;
                    $sub_has_voi = 1 if %vois and $vois{$unum}{$_};
                    $is_wt = 1 if ($_ eq $wtid); 
                }
                if ($type == 0){ #matches a non-wt sequence perfectly
                    $dist_from_wt = -1;
                }
                $exceeds_max_var_dist = 1 if $type > $max_var_dist;
                if ($dist_from_wt >= 0){ #doesn't match any sequence perfectly
                    if ($is_wt){ #if the closest match includes wt, add the distance to the total
                        $dist_from_wt += $type;
                        $sub_has_voi = 0; ## if the closest variant is the same distance as parental, will not count it towards variant of interest
                    } else { #otherwise, remove from consideration
                        $dist_from_wt = -1;
                    }
                }
                $has_voi = 1 if $sub_has_voi;
                if ($rep >= $min_rep){
                    $toprinttorout .= "\t$outtype" if $opt_d;
                    unless ($seg_seqs{$seg_num}){
                        $subseq_w =~ s/-//g;
                        push @{$seg_seqs{$seg_num}}, ([$wtid, $subseq_w]);
                    }
                    push @{$seg_seqs{$seg_num}}, ([$name, $degapped]);
                }
                $vartypes[$j][$var_seg_num] = "$outtype";
            }
            $var_seg_num++;
            $seq_i_start = $i;
            if (@sub_nc){
                $last_ustop = $ustop;
                ($unum, $ustart, $ustop, $uleng, $upart) = @{shift @sub_nc};
                @useqs = @{shift @sub_sq};
                $seg_num++;
            }
        }
        if ($pos_w == $ustart){
            next if $ustarted;
            my $dist = ($i - $seq_i_start) + 1;
            $subseq_w = substr($aln_w, $seq_i_start, $dist);
            $subseq_a = substr($aln_a, $seq_i_start, $dist);
            if ($subseq_w eq $subseq_a){
                $toprinttorout .= "\t" if $rep >= $min_rep and $opt_d;
            } else {
                if ($rep >= $min_rep){
                    $toprinttorout .= "\t*" if $opt_d;
                    unless ($seg_seqs{$seg_num}){
                        $subseq_w =~ s/-//g;
                        push @{$seg_seqs{$seg_num}}, ([$wtid, $subseq_w]);
                    }
                    $subseq_a =~ s/-//g;
                    push @{$seg_seqs{$seg_num}}, ([$name, $subseq_a]);
                }
            }
            $seq_i_start = $i + 1;
            $ustarted = 1;
            $seg_num++;
        }
    }
    if ($ustarted){
        my $dist = (($aleng - 1) - $seq_i_start);
        $subseq_w = substr($aln_w, $seq_i_start, $dist);
        $subseq_a = substr($aln_a, $seq_i_start, $dist);
        ## screen the sequences in the non-conserved region against the variant sequences
        $upart = 0 unless $upart;
        my ($degapped, $type, $list) = nc_screen($subseq_a, $upart, \@useqs);
        $varseqs[$j][$var_seg_num] = $degapped;
        my @larray = @{$list};
        if ($type == 0 and $larray[0] eq $wtid){
            $toprinttorout .= "\t" if $rep >= $min_rep and $opt_d;
        } else {
            my $outtype = "$type:" . join(",", @larray);
            my $sub_has_voi;
            foreach (@larray){
                $top_vars{$_}++;
                $sub_has_voi = 1 if %vois and $vois{$unum}{$_};
            }
            if ($type == 0){ #matches a non-wt sequence perfectly
                $dist_from_wt = -1;
            }
            $exceeds_max_var_dist = 1 if $type > $max_var_dist;
            if ($dist_from_wt >= 0){ #doesn't match any sequence perfectly
                my $is_wt;
                foreach (@larray){
                    $is_wt = 1 if ($_ eq $wtid); 
                }
                if ($is_wt){ #if the closest match includes wt, add the distance to the total
                    $dist_from_wt += $type;
                    $sub_has_voi = 0; ## if the closest variant is the same distance as parental, will not count it towards variant of interest
                    
                } else { #otherwise, remove from consideration
                    $dist_from_wt = -1;
                }
            }
            $has_voi = 1 if $sub_has_voi;
            if ($rep >= $min_rep){
                $toprinttorout .= "\t$outtype" if $opt_d;
                print $rout "\t$outtype" if $opt_d;
                unless ($seg_seqs{$seg_num}){
                    $subseq_w =~ s/-//g;
                    push @{$seg_seqs{$seg_num}}, ([$wtid, $subseq_w]);
                }
                push @{$seg_seqs{$seg_num}}, ([$name, $degapped]);
            }
            $vartypes[$j][$var_seg_num] = "$outtype";
        }
    } else {
        my $dist = (($aleng - 1) - $seq_i_start) + 1;
        $subseq_w = substr($aln_w, $seq_i_start, $dist);
        $subseq_a = substr($aln_a, $seq_i_start, $dist);
        #print STDERR "conserved:\n$subseq_w\n$subseq_a\n";
        if ($subseq_w eq $subseq_a){
            $toprinttorout .= "\t" if $rep >= $min_rep and $opt_d;
            #print $rout "\t" if $rep >= $min_rep and $opt_d;
        } else {
            if ($rep >= $min_rep){
                $toprinttorout .= "\t*" if $opt_d;
                print $rout "\t*" if $opt_d;
                #print STDERR "!!conserved region $last_ustop-$ustart does not match\n";
                unless ($seg_seqs{$seg_num}){
                    $subseq_w =~ s/-//g;
                    push @{$seg_seqs{$seg_num}}, ([$wtid, $subseq_w]);
                }
                $subseq_a =~ s/-//g;
                push @{$seg_seqs{$seg_num}}, ([$name, $subseq_a]);
            }
        }
    }
    if ($rep >= $min_rep and $opt_d){
        if (%vois){
            if ($has_voi){
                print $rout "*\t";
            } else {
                print $rout "\t";
            }
        }
        print $rout "$toprinttorout\n"
    }
    $has_voi_j[$j] = 1 if $has_voi;
    if ($dist_from_wt >= 1 and $dist_from_wt <= $max_wt_diff and !$exceeds_max_var_dist){
        ## There is no filtering for read count here.
        ## All reads deemed to be wt_adjacent will be counted regardless of
        ## the setting of $min_rep. This is probably the correct approach as
        ## they are all counted as non-unique parental sequences, the total
        ## of which is likely higher than the $min_rep cutoff. If we don't
        ## want to count these, however, we should uncomment the block
        ## surrounding these two lines:
        
        #if ($rep >= $min_rep){
        push @wt_adjacent_unique, $rep;
        $wt_adjacent_total += $rep;
        #}
        
        $wt_adjacent_j{$j} = 1;
    } else {
        if (%top_vars){
            my @order = sort {$top_vars{$b} <=> $top_vars{$a}} keys %top_vars;
            my $max = $top_vars{$order[0]};
            my @tmp;
            foreach my $var (@order){
                my $count = $top_vars{$var};
                last if $count < $max;
                push @tmp, $var;
            }
            @tmp = sort{$a cmp $b}@tmp;
            my $string = join(",", @tmp);
            @{$top_vars_j{$j}} = ($string, $max);
        } 
    }
    $exceeds_max_var_d{$j} = 1 if $exceeds_max_var_dist;
}
print STDERR "\rAnalyzing ", scalar(@strings), "\n\n";
close ($rout) if $opt_d;

($avg, $min, $max, $med, $num, $sum) = stats(\@wt_adjacent_unique);
statprint("unique sequences with no more than $max_wt_diff difference(s) from the parental sequence in the variable regions: $num\n");
statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
if ($min_rep > 1){
    ($avg, $min, $max, $med, $num, $sum) = stats(\@wt_adjacent_unique, $min_rep);
    statprint("unique sequences with no more than $max_wt_diff difference(s) from the parental sequence in the variable regions and with at least $min_rep read(s): $num\n");
    statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n"); 
}
print STDERR "\n";

if ($opt_d){
    for my $i (0 .. $#seg_names){
        next unless exists $seg_seqs{$i};
        my $name;
        unless ($name = $seg_names[$i]){
            $name = "cons:last";
        }
        $name =~ s/:/_/;
        open (my $sout, ">$pref.variants.$name.fasta");
        foreach my $slice (@{$seg_seqs{$i}}){
            my ($id, $seq) = @{$slice};
            print $sout ">$id\n$seq\n";
        }
        close ($sout);
    }
}

unless ($opt_d){
    #combine variant sequences 
    my @var_groups;
    foreach my $j (0 .. $#strings){
        my $rep = $str_counts[$j];
        my @vseq = @{$varseqs[$j]};
        my $vseq_stg = join("", @vseq);
        next if defined $wt_adjacent_j{$j};
        
        unless (@var_groups){
            push @var_groups, ([$rep, $vseq_stg, $j]);
            next;
        }
        my $found;
        for my $i (0 .. $#var_groups){
            if ($vseq_stg eq $var_groups[$i][1]){
                $var_groups[$i][0] += $rep;
                $found = 1;
            }
        }
        unless ($found){
            push @var_groups, ([$rep, $vseq_stg, $j]);
        }
    }
    #add wt_adjacent sequence counts to the count of pure parental sequences
    for my $i (0 .. $#var_groups){
        my ($rep, $vseq_stg, $j) = @{$var_groups[$i]};
        if ($vseq_stg eq $wt_vseq_stg){
            $var_groups[$i][0] += $wt_adjacent_total;
            last;
        }
    }
    my @var_group_counts = map $_->[0], @var_groups;
    ($avg, $min, $max, $med, $num, $sum) = stats(\@var_group_counts);
    statprint("After ignoring conserved region variants, $num unique sequences remain\n");
    statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
    if ($min_rep > 1){
        ($avg, $min, $max, $med, $num, $sum) = stats(\@var_group_counts, $min_rep);
        statprint("After ignoring conserved region variants, $num unique sequences with at least $min_rep read(s) remain\n");
        statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n"); 
    }
    print STDERR "\n";
    
    #my $clustered = scalar @var_groups;
    #print STDERR "After ignoring conserved region variants, $clustered unique sequences remain\n";
    #print $statout "After ignoring conserved region variants, $clustered unique sequences remain\n";
    
    
    #resort the list by read counts
    @var_groups = sort{$b->[0] <=> $a->[0]} @var_groups;
    
    open (my $rvout, ">$pref.results.only_variable.tsv");
    print $rvout "id\tcount\t";
    print $rvout "has_var\t" if $varfile;
    print $rvout join("\t", @var_seg_names), "\ttop_vars\ttop_var_count\n";
    open (my $kvout, ">$pref.sequence_key.txt");
    print $kvout "id\tlength\tnum\tseq\n";
    my @vseqs_for_fasta;
    my @vseq_stgs;
    my @vseq_ids;
    @var_group_counts = ();
    my @filt_counts;
    foreach my $i (0 .. $#var_groups){
        my ($rep, $vseq_stg, $j) = @{$var_groups[$i]};
        if ($exceeds_max_var_d{$j}){
            push @filt_counts, $rep;
            next;
        } else {
            push @var_group_counts, $rep;
        }
        next if $rep < $min_rep;
        my $tmp_name = "seq" . sprintf("%05d", $i + 1);
        push @vseq_stgs, $vseq_stg;
        push @vseq_ids, "$tmp_name#$rep";
        my $amp = $strings[$j];
        my $aleng = length($amp);
        print $kvout "$tmp_name\t$aleng\t$rep\t$amp\n";
        my @vseq = @{$varseqs[$j]};
        my @vtype;
        @vtype = @{$vartypes[$j]} if $vartypes[$j];
        print $rvout "$tmp_name\t$rep";
        if ($varfile){
            if ($has_voi_j[$j]){
                print $rvout "\t*";
            } else {
                print $rvout "\t";
            }
        }
        for my $k (0 .. $#var_seg_names){
            my $os = "\t";
            if ($vtype[$k]){
                $os .= $vtype[$k];
                push @{$vseqs_for_fasta[$k]}, ([$tmp_name, $vseq[$k]]);
            }
            print $rvout "$os";
        }
        my ($tv, $tvc) = ("-") x 2;
        if ($top_vars_j{$j}){
            ($tv, $tvc) = @{$top_vars_j{$j}};
            
        }
        print $rvout "\t$tv\t$tvc\n";
    }
    close ($rvout);
    close ($kvout);
    for my $i (0 .. $#var_seg_names){
        if ($vseqs_for_fasta[$i]){
            my $name = $var_seg_names[$i];
            $name =~ s/:/_/;
            open (my $vfasta, ">$pref.variants.only_variable.$name.fasta");
            my @wt = @{$varseqs_wt[$i]};
            print $vfasta ">", join("\n", @wt),"\n";
            foreach my $slice (@{$vseqs_for_fasta[$i]}){
                print $vfasta ">", join("\n", @{$slice}), "\n";
            }
            close ($vfasta);
        }
    }
    
    ($avg, $min, $max, $med, $num, $sum) = stats(\@var_group_counts);
    statprint("After ignoring reads with > $max_var_dist differences from parental or silent sequences, $num unique sequences remain\n");
    statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
    if ($min_rep > 1){
        ($avg, $min, $max, $med, $num, $sum) = stats(\@var_group_counts, $min_rep);
        statprint("After ignoring reads with > $max_var_dist differences from parental or silent sequences, $num unique sequences with at least $min_rep reads remain\n");
        statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
    }
    print STDERR "\n";
    ($avg, $min, $max, $med, $num, $sum) = stats(\@filt_counts);
    statprint("Number of removed unique sequences with > $max_var_dist differences from parental or silent sequences: $num\n");
    statprint("\tCounts: total:$sum avg:",sprintf("%.2f", $avg), "($min-$max), median:$med\n");
    
    ## hierarchical clustering of reads based on variant sequences
    print STDERR "Clustering...\n";
    my @v_ls_array;
    for my $i (0 .. ($#vseq_stgs-1)){
        my $s1 = $vseq_stgs[$i];
        for my $j ($i+1 .. $#vseq_stgs){
            my $s2 = $vseq_stgs[$j];
            #my $lev = levenshtein($s1, $s2);
            my $lev = distance($s1, $s2);
            my $norm_lev = 1 - ($lev / $maxdist); #normalizes the levenshtein distance to where 1 = no differences and 0 = maximum differences
            $v_ls_array[$i][$j] = $lev;
        }
    }
    my @v_outgroups = hclust(\@v_ls_array, \@vseq_ids);
}
close ($statout);
die "Done\n";

#----------------------------------------------------------------------------
sub statprint{
    my $string = join("", @_);
    print $statout "$string";
    print STDERR "$string";
    return();
}

sub levenshtein
{
    ## adapted from http://www.perlmonks.org/?node_id=245428
    
    # $s1 and $s2 are the two strings
    # $len1 and $len2 are their respective lengths
    #
    my ($s1, $s2) = @_;
    my ($len1, $len2) = (length $s1, length $s2);

    # If one of the strings is empty, the distance is the length
    # of the other string
    #
    return $len2 if ($len1 == 0);
    return $len1 if ($len2 == 0);

    # Init the distance matrix
    #
    # The first row to 0..$len1
    # The first column to 0..$len2
    # The rest to 0
    #
    # The first row and column are initialized so to denote distance
    # from the empty string
    #
    my @mat;
    @{$mat[0]} = (0 .. $len2);
    for my $i (1 .. $len1){
        @{$mat[$i]} = (0) x ($len2 + 1); #this line is not strictly necessary, but in Benchmarks it timed essentialy the same with or without, so I'll just leave it in
        $mat[$i][0] = $i;
    }

    # Some char-by-char processing is ahead, so prepare
    # array of chars from the strings
    #
    my @ar1 = split(//, $s1);
    my @ar2 = split(//, $s2);


    for my $i (1 .. $len1){
        for my $j (1 .. $len2){
            # Set the cost to 1 if the ith char of $s1
            # equals the jth of $s2
            # 
            # Denotes a substitution cost. When the char are equal
            # there is no need to substitute, so the cost is 0
            #
            my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;

            # Cell $mat{$i}{$j} equals the minimum of:
            #
            # - The cell immediately above plus 1
            # - The cell immediately to the left plus 1
            # - The cell diagonally above and to the left + the cost
            #
            # We can either insert a new char, delete a char or
            # substitute an existing char (with an associated cost)
            #
            $mat[$i][$j] = (sort{$a <=> $b}($mat[$i-1][$j] + 1,
                                            $mat[$i][$j-1] + 1,
                                            $mat[$i-1][$j-1] + $cost))[0];
        }
    }

    # Finally, the distance equals the rightmost bottom cell
    # of the matrix
    #
    # Note that $mat{$x}{$y} denotes the distance between the 
    # substrings 1..$x and 1..$y
    #
    return $mat[$len1][$len2];
}

my %nodes;
my %index;
sub hclust{
    my @array = @{$_[0]};
    my @names = @{$_[1]};
    
    %index = ();
    for my $i (0 .. $#names){
        $index{$names[$i]} = $i;
    }
    
    my $secount = scalar @names;
    return unless $secount >= 2;
    
    %nodes = ();
    my $node_num = 0;
    my @outgroups;
    
    my $dim = scalar @{$array[0]};
    my $cluster_number = 0;
    while ($dim > 1){
        my $max = 1E6;
        my @max_i;
        ## find maximum correlation(s), i.e. the lowest levenshtein score
        for my $i (0 .. $dim - 2){
            for my $j ($i + 1 .. $dim - 1){
                my $val = $array[$i][$j];
                if ($val < $max){
                    $max = $val;
                    @max_i = ("$i,$j");
                } elsif ($val == $max){
                    push @max_i, "$i,$j";
                }
            }
        }
        
        #return if $max < $min_sim;
        
        #if more than one pair had the maximum correlation, group samples in pairs into subclusters
        my %cluster_gens;
        my $subclustnum = 0;
        foreach (@max_i){
            my ($one, $two) = split(",", $_);
            my ($onesc, $twosc);
            $onesc = $cluster_gens{$one} if $cluster_gens{$one};
            $twosc = $cluster_gens{$two} if $cluster_gens{$two};
            if ($onesc and $twosc){
                ## reset subcluster on non-previously-overlapping cluster members to lowest subcluster value
                if ($onesc < $twosc){
                    foreach my $key (keys %cluster_gens){
                        if ($cluster_gens{$key} == $twosc){
                            $cluster_gens{$key} = $onesc;
                        }
                    }
                } else {
                    foreach my $key (keys %cluster_gens){
                        if ($cluster_gens{$key} == $onesc){
                            $cluster_gens{$key} = $twosc;
                        }
                    }
                }
            } elsif ($onesc){
                $cluster_gens{$two} = $onesc;
            } elsif ($twosc){
                $cluster_gens{$one} = $twosc;
            } else {
                $subclustnum++;
                $cluster_gens{$one} = $subclustnum;
                $cluster_gens{$two} = $subclustnum
            }
        }
        my %subclust;
        foreach my $key (sort keys %cluster_gens){
            my $sc = $cluster_gens{$key};
            push @{$subclust{$sc}}, $key;
        }
        
        ## create new distance array
        my @newgroups;
        my @newnames;
        foreach my $sc (sort keys %subclust){
            my @index = @{$subclust{$sc}};
            $node_num++;
            my $nodename = "NODE$node_num";
            $cluster_number++;
            $nodes{$nodename}{'max'} = $max;
            my @nstring;
            foreach (@index){
                my $mem_name = $names[$_];
                my $dist = 0;
                if ($nodes{$mem_name}){
                    $dist = $nodes{$mem_name}{'max'};
                    $mem_name = $nodes{$mem_name}{'string'};
                }
                my $outdist = $max - $dist;
                push @nstring, "$mem_name:$outdist";
                push @{$nodes{$nodename}{'members'}}, $names[$_];
            }
            $nodes{$nodename}{'string'} = "(" . join(",", @nstring) . ")";
            push @newgroups, [@index];
            push @newnames, $nodename;
            
        }
        for my $i (0 .. $dim - 1){
            next if exists $cluster_gens{$i};
            push @newgroups, ([$i]);
            push @newnames, $names[$i];
        }
        my @newarray;
        for my $i (0 .. $#newgroups - 1){
            my @cmem1 = @{$newgroups[$i]};
            for my $j ($i + 1 .. $#newgroups){
                my @cmem2 = @{$newgroups[$j]};
                my $newval;
                foreach my $k (@cmem1){
                    foreach my $l (@cmem2){
                        my $val;
                        $val = $array[$k][$l] if defined $array[$k][$l];
                        $val = $array[$l][$k] if defined $array[$l][$k];
                        if (defined $newval){
                            if ($cl_type eq "single"){
                                $newval = $val if $val >= $newval; #single linkage
                            } else {
                                $newval = $val if $val <= $newval; #complete linkage
                            }
                        } else {
                            $newval = $val;
                        }
                    }
                }
                $newarray[$i][$j] = $newval;
            }
        }
        @array = @newarray;
        @names = @newnames;
        
        last unless @array;
        
        $dim = scalar @{$array[0]};
    }
    
    my $tree = $nodes{"NODE$node_num"}{'string'};
    open (my $out, ">$pref.distance_dendrogram.tre") or die "ERROR: Can't open $pref.distance_dendrogram.tre for writing: $!\n";
    print $out "$tree;\n";
    close $out;
    
    ## if all of the elements were joined into a single group:
    unless (@outgroups){
        my @tmp = (0 .. $secount - 1);
        push @outgroups, ([@tmp]);
    }
    
    return (@outgroups);
}

sub nodeparse{
    my $nodename = $_[0];
    my @names;
    if (exists $nodes{$nodename}{'members'}){
        my @members = @{$nodes{$nodename}{'members'}};
        foreach my $mem (@members){
            next if exists $nodes{$mem}{'out'};
            if ($mem =~ m/^NODE\d+$/){
                my @tmp = nodeparse($mem);                
                push @names, @tmp;
            } else {
                push @names, $index{$mem};
            }
            $nodes{$mem}{'out'} = 1;
        }
    } else {
        push @names, $index{$nodename};
    }
    return(@names);
}

sub nw{
    #needleman-wunsch alignment
    #adapted from http://www.tcoffee.org/Courses/Exercises/Aln/Perl_examples/aln_algorithms.html
    my $seq1 = uc(shift);
    my $seq2 = uc(shift);
    
    #Parameters:
    my $match=10;
    my $mismatch=-10;
    my $gep=-10;
    
    $seq1 =~ s/\s//g;
    $seq1 =~ s/-//g;
    $seq2 =~ s/\s//g;
    $seq2 =~ s/-//g;
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    
    my @smat;
    my @tb;
    for my $i (0 .. $len1){
        $smat[$i][0] = $i * $gep; 
        $tb[$i][0] = 1;
    }
    for my $i (0 .. $len2){
        $smat[0][$i] = $i * $gep;
        $tb[0][$i] = -1;
    }
    
    for my $i (1 .. $len1){
        for my $j (1 .. $len2){
            #calculate the score
            my $score = $mismatch;
            $score = $match if (substr($seq1, $i-1, 1) eq substr($seq2, $j-1, 1));
            my $sub = $smat[$i-1][$j-1]+$score;
            my $del = $smat[$i][$j-1]+$gep;
            my $ins = $smat[$i-1][$j]+$gep;
            if ($sub > $del and $sub > $ins){
                $smat[$i][$j]=$sub;
                $tb[$i][$j]=0;
            } elsif ($del >$ins){
                $smat[$i][$j]=$del;
                $tb[$i][$j]=-1;
            } else {
                $smat[$i][$j]=$ins;
                $tb[$i][$j]=1;
            }
        }
    }
    my $i = $len1;
    my $j = $len2;
    my ($aln1, $aln2);
    my $dist = 0;
    while (!($i == 0 and $j == 0)){
        if ($tb[$i][$j] == 0){
            $aln1 .= substr($seq1, --$i, 1);
            $aln2 .= substr($seq2, --$j, 1);
            $dist++ if substr($aln1, -1, 1) ne substr($aln2, -1, 1);
        } elsif ($tb[$i][$j] == -1){
            $aln1 .= "-";
            $aln2 .= substr($seq2, --$j, 1);
            $dist++;
        } elsif ($tb[$i][$j] == 1){
            $aln1 .= substr($seq1, --$i, 1);
            $aln2 .= "-";
            $dist++;
        }
    }
    $aln1 = reverse($aln1);
    $aln2 = reverse($aln2);
    return($aln1, $aln2, $dist);
}

sub nc_screen {
    my $subseq_a = shift;
    my $upart = shift;
    my @useqs = @{$_[0]};
     
    unless ($subseq_a =~ m/^-+$/){
        $subseq_a =~ s/-//g;
    }
    ## first, look for perfect string or substring matches. This should be the case for most of the sequences
    my @perfect;
    foreach (@useqs){
       my ($uid, $useq) = @{$_};
        next unless $useq =~ m/\S/;
        if ($useq eq $subseq_a){
            push @perfect, $uid;
        }
    }
    
    my $outdist;
    my @outlist;
    if (@perfect){
        $outdist = 0;
        ## if the list includes the wild type sequence, will only output wild-type
        foreach (@perfect){
            if ($_ eq $wtid){
                push @outlist, $wtid;
                last;
            }
        }
        @outlist = @perfect unless @outlist;
    } else {
        ## Closest Levenshtein distance
        my @lds;
        foreach (@useqs){
            my ($uid, $useq) = @{$_};
            next unless $useq =~ m/\S/;
            my $ld = distance($subseq_a, $useq);
            push @lds, ([$uid, $ld]);
        }
        @lds = sort{$a->[1] <=> $b->[1]} @lds;
        my $minld = $lds[0][1];
        my @top;
        foreach my $slice (@lds){
            my ($uid, $ld) = @{$slice};
            last if $ld > $minld;
            push @top, $uid;
        }
        $outdist = $minld;
        @outlist = @top;
    }
    $subseq_a =~ s/-//g;
    return($subseq_a, $outdist, \@outlist);
}

sub stats {
    my @inlist = @{$_[0]};
    my $cutoff = $_[1] if defined $_[1];
    my $sum = 0;
    my @list;
    foreach my $val (@inlist){
        if (defined $cutoff){
            next if $val < $cutoff;
        }
        $sum += $val;
        push @list, $val;
    }
    my $num = scalar @list;
    return(0,0,0,0,0) unless $num;
    @list = sort{$a <=> $b}@list;
    my $average = $sum / $num;
    my $min = $list[0];
    my $max = $list[$#list];
    my $median;
    my $mid = int($num / 2);
    if ($num % 2 == 0){
        $median = ($list[$mid] + $list[$mid - 1]) / 2;
    } else {
        $median = $list[$mid];
    }
    return($average, $min, $max, $median, $num, $sum);
}

use Inline C => <<'END_OF_C_CODE';

void nw_c(char* s1, char* s2, SV* aln1, SV* aln2) {
    int match=10;
    int mismatch=-10;
    int gep=-10;

    int len1 = strlen(s1) + 1;
    int len2 = strlen(s2) + 1;
    
    int **smat = (int**)malloc(len1*sizeof(int*));
    for (int i = 0; i < len1; i++){
        smat[i] = (int*)malloc(len2*sizeof(int));
    }
    int **tb = (int**)malloc(len1*sizeof(int*));
    for (int i = 0; i < len1; i++){
        tb[i] = (int*)malloc(len2*sizeof(int));
    }
    
    for (int i=0; i<len1; i++){
        smat[i][0] = i * gep; 
        tb[i][0] = 1;
    }
    for (int i=0; i<len2; i++){
        smat[0][i] = i * gep;
        tb[0][i] = -1;
    }
    
    len1--;
    len2--;
    
    for (int i=1; i<=len1; i++){
        for (int j=1; j<=len2; j++){
            int score = mismatch;
            if (s1[i-1] == s2[j-1]){
                score = match;
            }
            int sub = smat[i-1][j-1]+score;
            int del = smat[i][j-1]+gep;
            int ins = smat[i-1][j]+gep;
            if (sub > del && sub > ins){
                smat[i][j]=sub;
                tb[i][j]=0;
            } else if (del > ins){
                smat[i][j]=del;
                tb[i][j]=-1;
            } else {
                smat[i][j]=ins;
                tb[i][j]=1;
            }
        }
    }
    
    for (int i = 0; i < (len1+1); i++){
        int* currentIntPtr = smat[i];
        free(currentIntPtr);
    }
    free(smat);
    
    int i = len1;
    int j = len2;
    
    char tmp1[len1 + len2];
    char tmp2[len1 + len2];
    
    int acount = 0;
    while (!(i == 0 && j == 0)){
        if (tb[i][j] == 0){
            tmp1[acount] = s1[--i];
            tmp2[acount] = s2[--j];
        } else if (tb[i][j] == -1){
            tmp1[acount] = '-';
            tmp2[acount] = s2[--j];
        } else if (tb[i][j] == 1){
            tmp1[acount] = s1[--i];
            tmp2[acount] = '-';
        }
        acount++;
    }
    tmp1[acount] = '\0';
    tmp2[acount] = '\0';
    
    for (int i = 0; i < (len1+1); i++){
        int* currentIntPtr = tb[i];
        free(currentIntPtr);
    }
    free(tb);
    
    int lens1 = strlen(tmp1);
    int lens2 = strlen(tmp2);
    
    char *p1 = tmp1;
    char *p2 = tmp1 + lens1 - 1;
    while (p1 < p2){
        char tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
    
    p1 = tmp2;
    p2 = tmp2 + lens2 - 1;
    while (p1 < p2){
        char tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
    
    sv_setpvn(aln1, tmp1, lens1);
    sv_setpvn(aln2, tmp2, lens2);
  
}

char *str_reverse_in_place(char *str, int len)
{
    char *p1 = str;
    char *p2 = str + len - 1;

    while (p1 < p2) {
        char tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }

    return str;
}

 
END_OF_C_CODE
