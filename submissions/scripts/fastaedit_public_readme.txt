FASTAEDIT_PUBLIC 

The fastaedit_public program allows users to trim ambiguous bases using the same parameters 
that are used in NCBI submission tools for assembled GenBank records.  The input file is a 
fasta file and you must specify the output file location in the command-line.  The default 
output file format is fasta but can be changed.  A sample command line for removing ambiguous 
characters is:

fastaedit_public -in YOUR_INPUT_FASTA -trim_ambig_bases -out_seq_file FILE_TO_WRITE_EDITED_FASTA_TO -out FILE_TO_WRITE_XML_VALIDATION_TO



COMMAND-LINE ARGUMENTS

1.   Input/Output Arguments
 -in <File_In>
   full path name to input fasta file (standard input by default)
   Default = `-'
 -out <File_Out>
   full path name to output XML log file (standard output by default)
   Default = `-'
 -out_seq_fmt <String, `asn', `fasta', `json', `xml'>
   Output format of the edited sequence file
   Default = `fasta'
 -out_seq_file <File_Out>
   full path name to edited output file

2.   Checking/Editing Arguments
 -trim_ambig_bases
   Trim ambiguous bases from 5’ and 3’ ends
 -remove_short_seqs <Integer>
   remove sequences that are lesser than given length
 -remove_long_seqs <Integer>
   remove sequences that are greater than given length
 -remove_ambig_seqs <Real>
   remove sequences where percent of ambiguous bases exceed given number


3.  Optional Arguments
 -h
   Print USAGE and DESCRIPTION;  ignore all other parameters
 -help
   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
 -xmlhelp
   Print USAGE, DESCRIPTION and ARGUMENTS in XML format; ignore all other
   parameters
 -logfile <File_Out>
   File to which the program log should be redirected
 -conffile <File_In>
   Program's configuration (registry) data file
 -version
   Print version number;  ignore other arguments
 -version-full
   Print extended version data;  ignore other arguments
 -version-full-xml
   Print extended version data in XML format;  ignore other arguments
 -version-full-json
   Print extended version data in JSON format;  ignore other arguments
 -dryrun
   Dry run the application: do nothing, only test all preconditions
