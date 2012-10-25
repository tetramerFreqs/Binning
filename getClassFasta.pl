#! /usr/bin/perl

=head1 NAME

	getClassFasta
		This program takes a chopped fasta file, a names file, a class file and extracts the seqs for the contigs presents in your desired class"; 

=head1 USAGE
	
	perl getClassFasta.pl -cls <CLASS File> -names <NAMES File> -fasta <Split FASTA File> -c <CLASS NUMBER>;

=head1 OPTIONS
	
	-cls		Required	cls file produced by the tetramer script.
	-fasta		Required	Chopped Fasta File produced by the tetramer script.
	-names	Required	names file produced by the tetramer script.
	-c	Required	Class number you wish to extract.
	
=head1 AUTHORS

	Sunit Jain, 2010 (sunitj [AT] umich [DOT] edu)
	
=head3 About the author

	www [DOT] sunitjain [DOT] com
	
=head1 LICENSE & COPYRIGHT

	This software is released into the public domain. To the extent 
	possible under law, all copyright and related or neighboring
	rights are waived and permission is explicitly and irrevocably
	granted to copy, modify, adapt and distribute this software
	in any way you choose. This program is distributed in the hope
	that it will be useful, but WITHOUT ANY WARRANTY; without even 
	the implied warranty of MERCHANTABILITY or FITNESS FOR A 
	PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my ($classFile, $namesFile,$fastaFile, $classNum);
my $version=$0." v1.0.0";

GetOptions(
	'fasta=s'=> \$fastaFile,
	'cls|class=s'=>\$classFile,
	'name|names=s'=>\$namesFile,
	'c=i'=>\$classNum,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=> sub{system('perldoc', $0); exit;},
);

if ((!$classFile) || (!$namesFile) || (!$fastaFile) || (!$classNum)){system('perldoc', $0); exit;}


# Parse *.cls file to get SeqID for all Seqs in the desired class
my %clsHash;
open ( CLS, $classFile) || die "ERROR: $classFile.\n".$!;
	while (my $line=<CLS>){
		chomp($line);
		unless ($line=~ /^%/){
			my ($seqNum, $cls)=split(/\t/,$line);
			if ($cls==$classNum){
				$clsHash{$seqNum}=$cls;	# %clsHash {Sequence Number  => Class Number}
			}
		}
	}
close CLS;

# Parse the *.names file to id the seq names in the fasta file using the classHash from above.
my %seqNames;
open (NAMES, $namesFile) || die "ERROR: $namesFile.\n".$!;
	while (my $line=<NAMES>){
		chomp($line);
		unless ($line =~ /^%/){
			my ($seqNum, $seqFastaName, $seqName)=split(/\t/, $line);
			if ($clsHash{$seqNum}){
				$seqNames{$seqFastaName}=$seqName; # %seqNames {Name of Seq part of the contig => Name of the whole contig}
			}
		}
	}
close NAMES;
undef %clsHash;

# Parse the split fasta file to get the desired sequences, using the seqNames hash from above.
my $outFile=$classNum.".fasta";
open (OUT,">".$outFile )|| die "ERROR: $outFile.\n".$!;
open (FASTA, $fastaFile) || die "ERROR: $fastaFile.\n".$!;
	$/= ">";
	while (my $line = <FASTA>) {
		chomp $line;
		next unless $line;
		my ($name, @sequence) = split (/\n/, $line);
		my $seq = uc(join ("", @sequence));
		if ($seqNames{$name}){
			print OUT ">$name\n$seq\n";
		}
	}
	$/= "\n";
close FASTA;
close OUT;
