#!/usr/bin/perl

use warnings;
use Pod::Usage;
use Getopt::Long;
use strict;
use Cwd 'abs_path';
use File::Path;
	
######################################################################
#
#	Variables initialization with default values
#
######################################################################

our ($record,@words,@unlinkfiles);

our ($Assembly,$Window,$Target_File_Out);

our ($Program_Folder_Path,$Target_Name,$Source_Data);

our ($Path_2_Wig,$Path_2_fasta);

our ($verbose,$help,$man);

our $Target_Filt_Path;

######################################################################
#
#	Defining options
#
######################################################################

GetOptions('verbose|v'=>\$verbose,'help|h'=>\$help,'man|m'=>\$man,) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 4 or pod2usage ("Syntax error");
	
######################################################################
#
#	Defining system variables
#
######################################################################

my ($myscriptname,$myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 1);

$Program_Folder_Path="$workingfolder";

print "Program folder is: $Program_Folder_Path\n";

$Source_Data=$ARGV[0];
$Target_Name=$ARGV[1];
$Window=$ARGV[2];
$Assembly=$ARGV[3];

$Target_Filt_Path="$Program_Folder_Path/data/targets/$Assembly/$Target_Name";




######################################################################
#
#	Removing empty lines from source file
#	and creating support temporary files
#
######################################################################

my $range=99999;
my $ID=int(rand($range));
my $filename="source.$ID";

if(-e $filename){ 
	$ID=$ID+100000;
}

system qq(awk NF $Source_Data > source.$ID);

######################################################################
#
#	Reading source file
#
######################################################################

open(CHECKBOOK,"source.$ID") || die "Couldn't open the source file!";

while($record=<CHECKBOOK>){
	chomp($record);
	@words=split(' ',$record);
	
	$Path_2_Wig=$words[0];
	$Path_2_fasta=$words[1];
			

######################################################################
#
#	Target initialization
#
######################################################################

print "Creating Consecutive and non-overlapping windows of $Window bp Size \n";
system qq(R --slave --args $Program_Folder_Path,$Target_Name,$Assembly,$Window,$Path_2_fasta < $Program_Folder_Path/lib/R/FilterTarget.R);


$Target_File_Out="$Program_Folder_Path/data/targets/$Assembly/$Target_Name/Filtered.txt";
print "Calculating Mappability and GC content...\n";
system qq($Program_Folder_Path/lib/bash/./TargetCreate.sh $Path_2_Wig $Target_File_Out $Program_Folder_Path $Assembly $Target_Name $Path_2_fasta);
print "...done!\n";

}


close(CHECKBOOK);
@unlinkfiles=("source.$ID");
unlink @unlinkfiles;

######################################################################
#
#	Documentation
#
######################################################################

=head1 SYNOPSIS 

 perl ReferenceWindowInitialize.pl [arguments] [options]

 Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.

 Function:
 
ReferenceWindowInitialize.pl initialises windows informations for further data processing with the XCAVATOR package. It requires 4 arguments (one source files - with space-delimited paths to source data for mappability and GC-content calculations), a label name, window size and assembly to run properly. A sub-folder with the specified target name will be created under "XCAVATOR/data/targets/hgXX".

 Example: 
 
 XCAVATOR> perl ReferenceWindowInitialize.pl SourceTarget.txt LabelName 50000 hg19
 
=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.


=back

=head1 DESCRIPTION

ReferenceWindowInitialize.pl is a Perl script which is part of the XCAVATOR package. It includes all of the first step operations of the XCAVATOR package. It calculates GC-content and Mappability for consecutive and non-overlapping windows of size L.

It requires, as arguments, the path to a source file (the default source file is "SourceTarget.txt" which is placed in the main XCAVATOR folder) containing the paths to source data (for the calculations of mappability and GC-content), the path to the target input file, a "target name", the window size and the assembly. Setting the label name as "MyLabel", all data calculated will be saved in the "MyLabel" folder in (if you are using the hg19 assembly) XCAVATOR/data/targets/hg19/MyLabel.

The allowed assemblies are hg19 and hg38.

XCAVATOR is freely available to the community for non-commercial use. For questions or comments, please contact "albertomagi@gmail.com".

=cut

