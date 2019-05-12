#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

my $script_dir = dirname($0);
$script_dir = "." if($script_dir eq "");

require( $script_dir . "/lib/Constant.pm" );
require( $script_dir . "/lib/Exception.pm" );
require( $script_dir . "/lib/Checker.pm" );
require( $script_dir . "/lib/Argument.pm" );
require( $script_dir . "/lib/Util.pm" );

my $PATTERN_TRANSCRIPT = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id


my $gff_file;
my $output_file;
my $accepted_features;
my $help;
my $verbose;

my %hash_features = ();

GetOptions(
	'i=s'    => \$gff_file,
	'a=s'    => \$accepted_features,
	'o=s'    => \$output_file,
	'v|verbose' => \$verbose,
	'h|help' => \$help
	
);$| = 1;    # Flush output

main();

sub main{
	if(defined $help){
		print Util::usage($0, "-i gff_input -a accepted_features -o report_output [-v|verbose -h|help]");
		die "\n";
	}
	
	my $error = Checker::check_argument(
		Argument->new(
			name     => "-i",
			value     => $gff_file,
			mandatory => Constant->YES,
			type      => Constant->FILE
		),
		Argument->new(
			name     => "-o",
			value     => $output_file,
			mandatory => Constant->YES,
			type      => Constant->FILE_TO_CREATE
		),
		Argument->new(
			name     => "-a",
			value     => $output_file,
			mandatory => Constant->YES,
			type      => Constant->STRING
		)
	);
	if($error){
		print Util::usage($0, "-i gff_input -a accepted_features -o output_file [-v|verbose -h|help]"), "\n";
		Exception::throw("Validating $0 arguments");
	}
	%hash_features = map{$_ => 1} split(/,/, $accepted_features);
	print(Util::show_running("Parsing $gff_file...")) if($verbose);
	read_gff_file();
	print(Util::show_finished("Parsing $gff_file...\n")) if($verbose);
	print(Util::show_finished("The script execution has been finished.\n")) if($verbose);
}

#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{
	if(open(IN, $gff_file)){
		open(OUT, ">$output_file");
		print OUT sprintf("%s\t%s\t%s\n", "CHROMOSOME", "GENE_ID", "BED_VALUE");
		while(<IN>){
			$_ =~ s/[\r\n]+$//g;
			unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#'
				my @cols = split /\t/;
				
				my $scaffold = $cols[0];
				my $attributes = $cols[-1];
				my @array_last_col = split(";", $attributes);
				my $bedvalues = $array_last_col[-1];
				$bedvalues =~ s/bedtools=//g;

				my $type = $cols[2];

				if (exists $hash_features{$type}) {
				  my $feature_id = get_feature_id($attributes, $PATTERN_TRANSCRIPT);
				  print OUT sprintf("%s\t%s\t%s\n", $scaffold, $feature_id,$bedvalues);
				}
			} # end unless(/^#/){
		  }# end while(<IN>){
		  close(OUT);
		  close(IN);
	}else{
		print STDERR "\n";
		if($? != 0){
			Exception::throw("An error has occurred when the script $0 tried to open the file $gff_file. Please check the messages returned by the operating system:", Exception->EXECUTION_ERROR, $?);
		}else{
			Exception::throw("An error has occurred when the script $0 tried to open the file $gff_file.", Exception->UNKNOWN_FILE_READING, $gff_file);
		}
	}
}

#########################################################################
# Extract the gene id based on the global variable pattern
#########################################################################
sub get_feature_id{
	my $attribute_text = shift;
	my $pattern = shift;
	my $regex = qr/$pattern/;
	if($attribute_text =~ $regex){
		$attribute_text =~ s/$regex/$2/g;
	}else{
		$attribute_text = "";
	}
	return $attribute_text;
}
