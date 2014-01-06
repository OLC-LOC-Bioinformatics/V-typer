#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Time::Piece;

my $path = getcwd;

# Initialise variables
my ($subtype, $strain, $strand, $start, $finish, $c, $d, $length, $e, $f, $g, $i, $accession, $sub, @subject_id, @unique_names, @blast, @subtype, @subtypes);

# Print a welcome message
print "Welcome to the CFIA V-Typer.\n";

# The files must be in the "sequences" subfolder. This subfolder must only have sequences that you wish to examine, or the program won't be able to find them
chdir ("$path/sequences");
my @files = glob("*");

print "Finding subtypes.\n";

foreach my $line (@files) {
	# The file extension is removed to keep filenames pretty. Make sure that the files have good names, as these names will be used in generating the report
	my ($name, $a) = split('\.', $line);
	# The famap program bundled with ePCR is used to create a mmapped-file with sequence data in random-accessible format
	system("$path/ePCR/famap -b $name.famap $line 2>/dev/null");
	# The fahash program precalculates positions of all words of the sequence database
	system("$path/ePCR/fahash -b $name.hash $name.famap 2>/dev/null");
	# Uses re-PCR from NCBI to quickly search the hashes created above for the presence of sequences present between the primer pairs supplied in the
	# vxt_subtyping_primers.txt file. If changes to the primer pairs are required, please edit this file.
	system ("$path/ePCR/re-PCR -S $name.hash -r + -m 10000 -n 0 -g 0 -G -q -o $name.blast $path/vtx_subtyping_primers.txt 2>/dev/null");
	
	# In order to keep re-PCR from printing some text to the console, I used "2>/dev/null" in the above command. This precluded the use of backticks, so a temporary file was created.
	open(INPUT, "<", "$name.blast");
	while (<INPUT>) {
		chomp;
		# This loop reads through the ".blast" file and extracts the required data.
		# Format of the output:
=com		
		#- sts	seq	strand	from	to	mism	gaps	act_len/exp_len
		vtx1a	OLC797_446_length_572088_cov_34.382633	+	467890	468367	0	0	478/10000-10000
		#   STS                                  vtx1a   CCTTTCCAGGTACAACAGCGGTT...430...CCAGAATGGCATCTGATGAGTTTCC  
		#                                                |||||||||||||||||||||||   430   |||||||||||||||||||||||||  
		#   Seq OLC797_446_length_572088_cov_34.382633 taCCTTTCCAGGTACAACAGCGGTT...430...CCAGAATGGCATCTGATGAGTTTCCtt
		##########################################################################
		#- Done
=cut		
		# If statement skips any lines that start with #'s 
		if (/^#/) {
			next;
		} else {
			# Split the text on "tabs" and assign the split text to appropriate variables
			($subtype, $strain, $strand, $start, $finish, $c, $d, $length) = split("\t", $_);
			# In the supplied primer file, the different primer pairs, which yield a similar subtype are differentiated by numbers (e.g. primer pairs vtx2a_1_1_1 and vtx2a_1_4_4)
			# may differ by a few degenerate bases, but are both vtx2a) - the underscores and numbers are discarded.
			($sub, $i) = split('\_', $subtype);
			push(@subtype, $sub);			
		}
	}
		# Because a file can have either zero, one, or multiple verotoxin subtypes, this loop counts the length of @subtype, and acts accordingly based on the length
		# Foreach type, the results are printed to the terminal, and pushed to a summary array, and @subtype is cleared
		# No subtypes
		if (scalar @subtype == 0) {
			print "$name has no subtype.\n";
			push (@subtypes, "$name\n");
			undef @subtype;
		# Greater than one subtype - the @subtype array is printed as a comma-separated list using the join command
		} elsif (scalar @subtype > 1) {
			chomp @subtype;
			print "$name has subtypes ", join( ', ', @subtype),".\n";
			push (@subtypes, "$name\t", join( "\t", @subtype),"\n");
			undef @subtype;
		# One subtype
		} else {
			print "$name has subtype $sub.\n";
			push (@subtypes, "$name\t$sub\n");
			undef @subtype;
		}
	close INPUT;
}

# Remove all the files created during the analysis
unlink glob ("*.blast");
unlink glob ("*.hash");
unlink glob ("*.famap");

# Create and open the output file
chdir ("$path");
# Get the data for naming the output file
my $date = Time::Piece->new->strftime('%Y%m%d');
open (OUTPUT, ">", "vt_subtyping_results_$date.csv");
print OUTPUT "Strain\t", "Subtype\n";
# Print the values from the array to the file
foreach (@subtypes) {
	print OUTPUT "$_";
}
close OUTPUT;

print "Verocytoxin subtyping completed.\n";

exit;











