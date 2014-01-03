#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
my $path = getcwd;

# User is prompted to enter the name of the sequences against which the vtx subtyping scheme is to be applied
print "Welcome to the CFIA V-Typer.\n";
print "Please input the name of the sequence(s) for which you would like to determine the verotoxin subtype.\n";

my $input = <>;
chomp $input;

print "Thank you.\n";

open (INPUT, "<", $input) or die "Sorry. There was an error opening file $input: $!.\n";
my @data = (<INPUT>);
close INPUT;

# The file extension is removed to keep filenames pretty
my ($name, $a) = split('\.', $input);

# In order from making the same files every time, the unless statement is used
unless (-e "$name.famap") {
	# The famap program bundled with ePCR is used to create a mmapped-file with sequence data in random-accessible format
	system("$path/ePCR/famap -b $name.famap $input");
	# The fahash program precalculates positions of all words of the sequence database
	system("$path/ePCR/fahash -b $name.hash $name.famap");
}

print "Finding subtypes.\n";

my $out = "$name.blast";

my ($subtype, $strain, $b, $start, $finish, $c, $d, $length, $e, $f, $g, $h, $i, $accession, $ident, $sub, @subject_id, @unique_names);
my %seen = ();

system("$path/ePCR/re-PCR -S $name.hash -r + -m 10000 -n 0 -g 0 vxt_subtyping_primers.bak.txt -o $out");

open(BLAST, "<", $out);

while (<BLAST>) {
	chomp;
	if (/^#/) {
		next;
	} else {
		($subtype, $strain, $b, $start, $finish, $c, $d, $length) = split("\t", $_);
		#($e, $f, $g, $accession) = split('\|', $strain);
		($ident, $h) = split('\_', $strain);
		($sub, $i) = split('\_', $subtype);
		#print "$subtype, $strain, $b, $start, $finish, $c, $d, $length.\n";
		print "$ident has subtype $sub.\n";
		
	}
}














