#! /usr/bin/perl

=head1 NAME

	tetramer_freqs_esom.pl 
		calculates Tetramer frequencies for the given fasta file and produces 4 output files

=head1 USAGE
	
	perl tetramer_freqs_esom.pl -f fastaFile -a annotationFile [OPTIONS]

=head1 OPTIONS
	
	-f		Required	Fasta File, may include X:s and N:s
	-a		Required	Annotation File (3 columns; 1. full contig name, 2.annotation, 3.Class Number (Your metagenome has class#0, everything else 1+)) 
	-min	Optional	default=2500; Minimal length (in nt) of input contig to be included in output
	-max	Optional	default=5000
	Note:	The script will split sequence after each 'max' nt; join last part, if remaining seq shorter than 'max', with second-last part
			eg: in default settings, a sequence of 14 kb will be split into a 5 kb and a 9 kb fragment if window_size = 5 kb.

=head1 AUTHORS

	#Anders Andersson, 2007 (anders.andersson@scilifelab.se)
	
=head2 Modified by

	#Sunit Jain, 2010 (sunitj [AT] umich [DOT] edu)
	
=head1 CITATION

	If you use this script please cite:
	Dick, G.J., A. Andersson, B.J. Baker, S.S. Simmons, B.C. Thomas, A.P. Yelton, and J.F. Banfield (2009).
	Community-wide analysis of microbial genome sequence signatures. Genome Biology, 10: R85.
	
=head1 LICENSE & COPYRIGHT

	Copyright (C) 2007 Anders Andersson (anders.andersson@scilifelab.se)
	This is free software; see the COPYRIGHT file accompanying this script
	for details. There is NO warranty; not even for MERCHANTABILITY or 
	FITNESS FOR A PARTICULAR PURPOSE.

=cut


use Getopt::Long;
use File::Basename;

my $version=$0." v1.0.3";
my $sfile; #fasta file, may include X:s and N:s
my $annotationfile; #full contig name in left, annotation in right, column. headers (whatever) on first line 
my $min_length = 2500; #Minimal length (in nt) of input contig to be included in output
my $window_size = 5000; #split sequence after each window_size nt, 
                     #join last part, if shorter than window_size, 
                     #with second-last part (a sequence of 
                     #14 kb will be split into a 5 kb and a
                     #9 kb fragment if window_size = 5 kb)

GetOptions(
	'f|fasta=s'=> \$sfile,
	'a|ann=s'=>\$annotationfile,
	'min:i'=>\$min_length,
	'max:i'=>\$window_size,
	'v|version'=> sub{&licensing; exit},
	'h|help'=> sub{system('perldoc', $0); exit;},
);

&licensing;

if ((! $sfile) || (! $annotationfile)){print "[ERROR $0] Missing required input.\nFor help using the script, type 'perl $0 -help'\n"; exit;}

print "\n############### TetraNucleotide Frequencies ##################\n";
print "Minimum length (in bases) of input contig to be included in output:\n";
print "$min_length\n";
print "Window Size:\n$window_size\n";
my $seqfile = basename($sfile, ".fasta");
print $seqfile;
#!!! Program will automatically create outfiles called (whatever) infile.lrn and infile.names (and overwrite existing) !!!

$lrnfile = "Tetra_".$seqfile."_$min_length\.lrn";
$namesfile = "Tetra_".$seqfile."_$min_length\.names";
$classfile = "Tetra_".$seqfile."_$min_length\.cls";
$reffile= "Tetra_".$seqfile."_$min_length_$window_size\.fasta";

#$window_size = 5000;

$mer_length = 4; #can't change!
####### main #############
$n=0;
&make_list_of_possible_tetramers;
&calc_tetra_freqs;
&make_lrn_file;
&make_names_file;
&make_class_file;
&make_seq_file;
&getRowColESOM;

##### sub routines #######
sub make_list_of_possible_tetramers {
   #only works if mer_length == 4 ...
   @bases = ("A", "T", "C", "G");
   foreach $base1 (@bases) {
       foreach $base2 (@bases) {        
           foreach $base3 (@bases) {        
               foreach $base4 (@bases) {
                   $mer = $base1.$base2.$base3.$base4;
                   $rc_mer = &make_revcomp($mer);
                   if (!defined $allowed{$rc_mer}) {
                       push (@mers, $mer);
                       $allowed{$mer}++;
                   }
               }
           }
       }
   }
}

sub make_lrn_file {
    print "printing lrn file:   $lrnfile\n";
    open (OUT, ">$lrnfile") || "can't create outfile $lrnfile";
    $number_rows = @names;
    $number_cols = @mers + 1;
    print OUT "% $number_rows\n";
    print OUT "% $number_cols\n";
    print OUT "% 9";
    foreach $mer (@mers) {
        print OUT "\t1";
    }
    print OUT "\n";
    print OUT "% Key";
    foreach $mer (@mers) {
        print OUT "\t$mer";
    }
    print OUT "\n";
    $key = 0;
    foreach $tetra (@tetras) {
        $key++;
        print OUT "$key$tetra\n";
    }
    close (OUT);
}

sub make_names_file {
    print "printing names file: $namesfile\n";
    $number_rows = @names;
    open (OUT, ">$namesfile");
    print OUT "% $number_rows\n";
    foreach $name (@names) {
        print OUT "$name\n";
    }
    close (OUT);
}

sub make_class_file {
    open (INFILE, $annotationfile) || die ("can't open infile!");
    while (<INFILE>) {
        next if $_=~ /^#/;
        chomp;
        @fields = split(/\t/, $_);
        $fields[0] =~ s/\s+$//;
        $class{$fields[0]} = $fields[2];
    }
    close (INFILE);
    print "printing class file: $classfile\n";
    open (OUT, ">$classfile");
    $number_rows = @names;
    print OUT "% $number_rows\n";
    foreach $item (@names) {
        #print "$item\n";
        @fields = split(/\t/, $item);
        print OUT "$fields[0]\t$class{$fields[2]}\n";
    }
    close (OUT);
}
sub make_seq_file {
    print "printing seq file: $reffile\n";
    $number_rows = @names;
	$getseq= @seq_list;
	open (OUT, ">$reffile");
    print OUT "% $number_rows\n";

	foreach $item (@names) {
        @fields = split(/\t/, $item);
        print OUT ">$fields[1]\n";
		print OUT "$seq_list[$n]\n";
		$n++;
	}
    close (OUT);
}
sub calc_tetra_freqs {
    print "\ncalculating tetranucleotide frequencies ...\n";
    @names = ();
    @tetras = ();
    $total_index = 0;
    $read = 0;
    $seq = "";
    open (INFILE, $sfile) || die ("can't open infile!");
    while (<INFILE>) {
        chomp;
        if (substr($_, 0 ,1) eq ">") {
            if ($read == 0) {
                $read++;
                $id = $_;
                substr($id, 0, 1) = "";
                $id =~ s/\s+$//;
                next;
            }
            if (length($seq) >= $min_length) {
                &get_tetra_freqs;
            }
            $id = $_;
            substr($id, 0, 1) = "";
            $id =~ s/\s+$//;
            $seq = "";
            $read++;
        } else {
            $seq = $seq.$_;
        }
    }
    if (length($seq) >= $min_length) {
    &get_tetra_freqs;
    }
    close (INFILE);
}

sub get_tetra_freqs {
    #filter out short sequences between N's and X's 
    #as well as between these and beginning and end of sequence
    @lowqual = ();
    push(@lowqual, 0); #to get start position
    for ($i = 0; $i < (length($seq) - $mer_length); $i++) {
        $base = substr($seq, $i, 1);
        if ($base eq "N" || $base eq "X") {
            push(@lowqual, $i);
        }
		
	}
    push(@lowqual, length($seq)); #to get end position
    $filtered_seq = $seq;
    for ($i = 1; $i < @lowqual; $i++) {
        $length = $lowqual[$i] - $lowqual[$i - 1] - 1;
        if ($length < 50) {
            for ($j = $lowqual[$i - 1]; $j < $lowqual[$i]; $j++) {
                substr($filtered_seq, $j, 1) = "Z";
            }
        }
    }

    $seq = $filtered_seq;
    $seq = uc($seq);
    
    @sub_seq = ();
	#@seq_list=();
    if (length($seq) < 2*$window_size) {
        @sub_seq = ($seq);
		#push(@seq_list, $seq);
    } 
	else {
        for ($i = 0; $i < length($seq); $i = $i + $window_size) {
            $subseq = substr($seq, $i, $window_size);
            push(@sub_seq, $subseq);
		#	push(@seq_list, $subseq);
		}
		if (length($sub_seq[-1]) < $window_size) {
            $sub_seq[-2] = $sub_seq[-2].$sub_seq[-1];
            pop (@sub_seq);
		}
	}
	#calculate and print freqs for each subsequence
    $sub_index = 0;
    foreach $seq (@sub_seq) {
		
        $sub_index++;
        %this_mers = ();
        $sum = 0;
        for ($i = 0; $i < length($seq); $i++) {
            $mer = substr($seq, $i, $mer_length);
            if (defined $allowed{$mer}) {
                $this_mers{ $mer }++;
                $sum++;
            } else {
                $rc_mer = &make_revcomp($mer);
                if (defined $allowed{$rc_mer}) {
                    $this_mers{ $rc_mer }++;
                    $sum++;
                }
            }
        }
		
        #if ($sum > 5000) {
            $tetra = "";
            $total_index++;
            #$name = "$total_index\t$id"."_"."$sub_index\t$id"."_"."$sub_index";
            $name = "$total_index\t$id"."_"."$sub_index\t$id";
            push(@names, $name);
			push(@seq_list, $seq);
            foreach $mer (@mers) {
                if (defined $this_mers{$mer}) {
                    $counts = $this_mers{$mer}/$sum;
                    #print"\t$counts";
                    $tetra = $tetra."\t".$counts;
                } else {
                    #print"\t0";
                    $tetra = $tetra."\t0";
                }
            }            
            #print"\n";
            push(@tetras, $tetra);
        #}
    }
}

sub make_revcomp {
    local($seq) = $_[0];
    local($modseq) = "";
    local($i);
    local($nt);
    for ($i = 0; $i < length($seq); $i++) {
        $nt = substr($seq, $i, 1);
        if ($nt eq "A") {
             $modseq = T.$modseq;
        } elsif ($nt eq "T") {
             $modseq = A.$modseq;
        } elsif ($nt eq "C") {
             $modseq = G.$modseq;
        } elsif ($nt eq "G") {
             $modseq = C.$modseq;
        } else {
             $modseq = $nt.$modseq;
        }
    }
    return ($modseq);
}

sub getRowColESOM{
	my $acceptedSeq=@names;

	my $mapSpace= $acceptedSeq * 5.5;
	my $sqrt= sqrt($mapSpace);
	my $half= $sqrt/2;
	my $cols= (int(($half + $sqrt) + 0.5))*2;
	my $rows= (int(($mapSpace/$cols) + 0.5))*2;

	print "\nTry the following values for ESOM Training:\n\tRows:\t$rows\n\tCols:\t$cols\n";
	print "These values are just meant as suggestions, feel free to try your own\n";
}

sub licensing{
	print $version."\n\n";
	
	print "Copyright (C) 2007 Anders Andersson (anders.andersson@scilifelab.se)\nThis is free software; see the COPYRIGHT file accompanying this script for details. There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n";

	print "Please cite:\nDick, G.J., A. Andersson, B.J. Baker, S.S. Simmons, B.C. Thomas, A.P. Yelton, and J.F. Banfield (2009).  Community-wide analysis of microbial genome sequence signatures. Genome Biology, 10: R85.\n\n";
}
