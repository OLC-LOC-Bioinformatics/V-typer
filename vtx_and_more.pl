#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Time::Piece;
use List::MoreUtils qw(uniq);
use Getopt::ArgParse;

# Adds an argument parser to allow the user to easily specify certain variables
my $ap = Getopt::ArgParse->new_parser(
        prog        => $0,
        description => 'Performs an automated ePCR using user supplied primers. The primers must be in the format: 
        				<name>TAB<forward primer>TAB<reverse primer>TAB<maximum size allowed between primers>NEWLINE.
        				Sequence files must be stored in <path>/sequences',
 );

# Set the accepted arguments
$ap->add_arg('--path', '-p', required => 1, help => 'The path of the folder that contains the run data');
$ap->add_arg('--primerFile', '-f', required => 0, default => "vtx_subtyping_primers.txt", 
			 help => 'The name of the file containing the primers. If not supplied, a default of "vtx_subtyping_primers.txt" is used.');
$ap->add_arg('--reportName', '-r', required => 0, default => 'vtx_subtyping_results', help => 'The name to use for the report. Defaults to "vtx_subtyping_results"'); 

# Create a variable to store the arguments
my $ns = $ap->parse_args();

# Populate the variables from the arguments
my $path = $ns->path;
my $primers = $ns->primerFile;
my $reportName = $ns->reportName;

# This start time will be used in calculating the total time of the run
my $start_time = time;

# Move to the path
chdir($path);

# Initialise variables
my ($subtype, $strain, $strand, $start, $finish, $c, $d, $length, $e, $f, $g, $i, $accession, $sub, @subject_id, @unique_names, @blast, @subtype, @subtypes);

# Print a welcome message
print "Welcome to the CFIA V-Typer.\n";

# Open the primer file
# Quick check to see if the file exists
open (INPUT, "<", $primers) or die "File does not exist! \n";

# The files must be in the $path/"sequences" subfolder. This subfolder must only have sequences that you wish to examine, or the program won't be able to find them
# Additionally, the file extensions must start with an "f"
chdir ("$path/sequences");
my @files = glob("*.f*");

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
	# my $exec = "$repcr -s $hash -n $mismatches -g 0 $forwardPrimer $reversePrimer 100-2000 -G > $outfile";
	system ("$path/ePCR/re-PCR -S $name.hash -r + -m 10000 -n 1 -g 0 -G -q -o $name.blast $path/$primers 2>/dev/null");
	
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
	# Create a sorted array of the uniques values in @subtype 	
	my @unique = sort(uniq(@subtype));
	# Because a file can have either zero, one, or multiple verotoxin subtypes, this loop counts the length of the unique entries in @subtype, 
	# and acts accordingly based on the length. For each type, the results are printed to the terminal, and pushed to a summary array
	# No subtypes
	if (scalar @unique == 0) {
		print "$name has no match.\n";
		push (@subtypes, "$name\t-\n");
	# Greater than one subtype - the @subtype array is printed as a comma-separated list using the join command
	} elsif (scalar @unique > 1) {
		print "$name has subtypes ", join(", ", @unique) . ".\n";
		push (@subtypes, "$name\t" . join("\t", @unique) . "\n");
	# One subtype
	} else {
		print "$name has subtype $sub.\n";
		push (@subtypes, "$name\t$sub\n");
	}
	# Clear @subtype for the next iteration	
	undef @subtype;
close INPUT;
}

# Remove all the files created during the analysis
unlink glob ("*.blast*");
unlink glob ("*.hash*");
unlink glob ("*.famap*");

# Create and open the output file
chdir ("$path");
# Get the data for naming the output file
my $date = Time::Piece->new->strftime('%Y%m%d%_H%M%S');
open (OUTPUT, ">", "reports/$reportName" ."_$date.csv");
print OUTPUT "Strain\t", "Subtype\n";

# Print the values from the array to the file
foreach (uniq(@subtypes)) {
	print OUTPUT "$_";
}
close OUTPUT;

# Return the run time
my $end_time = time;
my $total_time = $end_time - $start_time;
my $legible_time = sprintf("%.1d", $total_time);

print "-------------------------------------------\n";
print "Verocytoxin subtyping completed.\n";
print "The total run time was $legible_time seconds.\n\n";

exit;