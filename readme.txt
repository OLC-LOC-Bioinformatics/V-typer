The CFIA VT-Typer v0.1

This program analyses fasta-formatted files for the presence of primer pairs that are diagnostic of verotoxin subtypes. It uses the ePCR tool developed by NCBI, and a custom primer list developed by Flemming Scheutz at the WHO.

To use the program:

1) Put your fasta-formatted files in the "sequences" folder - ensure that the names of the files are informative, as these names will be used in the generation of the report.
2) Ensure that both the perl script (vtx2_blast.pl) and all binary files in the ePCR subfolder have the appropriate permissions (e.g. allow executing file as program)
3) run "perl vtx2_blast.pl" on the command line
