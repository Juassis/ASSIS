#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;
use Getopt::Long;

my $script_dir = dirname($0);
$script_dir = "." if($script_dir eq "");

require( $script_dir . "/lib/Constant.pm" );
require( $script_dir . "/lib/Exception.pm" );
require( $script_dir . "/lib/Checker.pm" );
require( $script_dir . "/lib/Argument.pm" );
require( $script_dir . "/lib/Util.pm" );


my $input_file;
my $output_file;
my $help;
my $verbose;


###########################################################################
#receiving parameters
# to do: implement parameter h or help
###########################################################################
GetOptions(
  'f=s'    => \$input_file,
  'o=s'    => \$output_file,
  'v|verbose' => \$verbose,
  'h|help' => \$help
);$|=1;


main();

sub main{
	if(defined $help){
		print Util::usage($0, "-f fasta_input -o fasta_output [-v|verbose -h|help]");
		die "\n";
	}
	
	my $error = Checker::check_argument(
		Argument->new(
			name     => "-f",
			value     => $input_file,
			mandatory => Constant->YES,
			type      => Constant->FILE
		),
		Argument->new(
			name     => "-o",
			value     => $output_file,
			mandatory => Constant->YES,
			type      => Constant->FILE_TO_CREATE
		)
	);
	if($error){
		print Util::usage($0, "-f fasta_input -g gff_input -o output_file"), "\n";
		Exception::throw("Validating $0 arguments");
	}
	
	print(Util::show_running("Parsing $input_file and generating report...")) if($verbose);
	my $in = Bio::SeqIO->new(-file => "<$input_file");
	
	if(open(OUT, ">".$output_file)){
		print OUT sprintf("%s\n", join("\t", qw/ID	Length StartCheck StopCheck InternalStops/));
		while(my $s = $in->next_seq){
			my $id = $s->primary_id;
			my $alpha = $s->alphabet;
			my $aa = $s->seq;

			my $start = ($aa =~ /^M/) || 0;
			my $end = ($aa =~ /([^A-Z])$/) || 0;

			my @stops = ($aa =~ /(\.|\*)/g);
			my $inner_stop; #in-frame stop codons

			if($end){
				$inner_stop = scalar @stops - 1;
			}else{
				$inner_stop = scalar @stops;
			}

			print OUT sprintf("%s\t%d\t%s\t%s\t%s\n", $id, $s->length, $start, $end, $inner_stop);
		}
		close(OUT);
	}else{
		if($? != 0){
			Exception::throw("An error has occurred when the script $0 tried to create/open the file $output_file. Please check the messages returned by the operating system:", Exception->EXECUTION_ERROR, $?);
		}else{
			Exception::throw("An error has occurred when the script $0 tried to create/open the file $output_file", Exception->UNKNOWN_FILE_WRITING, $output_file);
		}
	}
	print(Util::show_finished("Parsing $input_file and generating report...\n")) if($verbose);
	print(Util::show_finished("The script execution has been finished.\n")) if($verbose);
}
