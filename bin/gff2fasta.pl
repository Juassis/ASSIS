#!/usr/bin/env perl

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Basename;
use Getopt::Long;
use Term::ANSIColor;

my $script_dir = $0;
my $script_dir = dirname($0);
$script_dir = "." if($script_dir eq "");

require( $script_dir . "/lib/Constant.pm" );
require( $script_dir . "/lib/Exception.pm" );
require( $script_dir . "/lib/Checker.pm" );
require( $script_dir . "/lib/Argument.pm" );
require( $script_dir . "/lib/Util.pm" );


my $PATTERN_CDS = "\^.*ID=(CDS:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $PATTERN_TRANSCRIPT = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $PATTERN_PARENT_TRANSCRIPT = "\^.*Parent=(transcript:)?([a-zA-Z0-9\-_\.]+);?.*\$"; # pattern to extract the gene id

my $help;
my $verbose;

my $file_fasta;
my $gff_file;
my $output_file;

my %hash_transcript = ();
my %hash_gff = ();
my $input_database;
my $output_obj;

GetOptions(
  'f=s'    => \$file_fasta,
  'g=s'    => \$gff_file,
  'o=s'    => \$output_file,
  'v|verbose' => \$verbose,
  'h|help' => \$help
); $| = 1;

main();

sub main{
	if($help){
		print Util::usage($0, "-f fasta_input -g gff_input -o output_file [-v|verbose -h|help]");
		die "\n";
	}

	# Before the pipeline logic we need to validate the arguments sent by the user
	my $error = Checker::check_argument(
		Argument->new(
			name     => "-f",
			value     => $file_fasta,
			mandatory => Constant->YES,
			type      => Constant->FILE
		),
		Argument->new(
			name     => "-g",
			value     => $gff_file,
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
	
#    my $output_file_name = $file_fasta;
#    $output_file_name =~ s/^(.+\/)+//g;
	$output_obj = Bio::SeqIO->new( -format => 'fasta', -file => ">$output_file", -verbose=>-1 );

	###### Output type description ######
	# cds - translated sequence (starting with ATG and ending with a stop codon included)
	# protein - cds translated (includes a * as the stop codon)
	print(Util::show_running("Parsing $file_fasta file...")) if($verbose);
	$input_database = Bio::DB::Fasta->new($file_fasta);
	print(Util::show_finished("Parsing $file_fasta file...\n")) if($verbose);
	
	print(Util::show_running("Parsing $gff_file file...")) if($verbose);
	read_gff_file();
	print(Util::show_finished("Parsing $gff_file file...\n")) if($verbose);
	
	print(Util::show_running("Exporting translated sequences to $output_file...")) if($verbose);    
	sort_gff();
	print(Util::show_finished("Exporting translated sequences to $output_file...\n")) if($verbose);
	print(Util::show_finished("The script execution has been finished.\n")) if($verbose);
}

sub sort_gff{
  my @scaffolds = keys %hash_gff;
  for(@scaffolds){
	my $scaffold = $_;
	my @array_genes = sort {$hash_gff{$scaffold}{$a}{"MIN"} <=> $hash_gff{$scaffold}{$b}{"MIN"}} keys %{$hash_gff{$scaffold}};
	for(@array_genes){
	  my $parent_id = $_;
		my @array_cds = sort {$a->{"START"} <=> $b->{"START"}} @{$hash_gff{$scaffold}{$parent_id}{"CDS"}};
		my $mergedCDS_seq = "";
		my $strand = "+";
		for(@array_cds){
		  $strand = $_->{"STRAND"};
		  my $cds_seq = $input_database->seq( $scaffold, $_->{"START"}, $_->{"END"} );
		  $mergedCDS_seq .= $cds_seq;
		}
		
		my $output_cds = Bio::Seq->new(
			-seq        => $mergedCDS_seq,
			-id         => $parent_id,
			-display_id => $parent_id,
			-alphabet   => 'dna',
		);
		if ($strand eq "-") {
			$output_cds = $output_cds->revcom();
		}
		#translate CDS to peptide for protein sequence
		my $output_pep = $output_cds->translate();
		#write to file
		if (length($mergedCDS_seq) != 0) {
			  $output_obj->write_seq($output_pep);
		}

	} # end 	for(@array_genes){
  }# end for(@scaffolds){
}

#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{

  if(open(IN, $gff_file)){
	my $transcript_id = "";
	  while(<IN>){
		$_ =~ s/[\r\n]+$//g;
		
		unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#'
			my $line = $_;
			my @cols = split /\t/;

			my $gene_name = $cols[-1]; # The gene name is at last column
			my $scaffold = $cols[0]; # The scaffold name is the first column
			my $type = $cols[2];
			my $gene_start = $cols[3];
			my $gene_end = $cols[4];
			my $strand = $cols [6];

			my $parent_id = get_feature_id($gene_name, $PATTERN_PARENT_TRANSCRIPT);

			if (lc($type) eq 'cds') {
				my $feature_id = get_feature_id($gene_name, $PATTERN_CDS);

				my %hash_cds = ();
				$hash_cds{"START"} = int($gene_start);
				$hash_cds{"END"} = int($gene_end);
				$hash_cds{"STRAND"} = $strand;
				$hash_cds{"ID"} = $feature_id;
				$hash_cds{"TYPE"} = $transcript_id;

				if(exists $hash_gff{$scaffold}{$parent_id}){
				my $array_cds = $hash_gff{$scaffold}{$parent_id}{"CDS"};
				if($gene_start < $hash_gff{$scaffold}{$parent_id}{"MIN"}){
				  $hash_gff{$scaffold}{$parent_id}{"MIN"} = $gene_start;
				}
				push(@$array_cds, \%hash_cds);

				}else{
					$hash_gff{$scaffold}{$parent_id}{"MIN"} = $gene_start;
					my @array_cds = (\%hash_cds);
					$hash_gff{$scaffold}{$parent_id}{"CDS"} = \@array_cds;
				}
				}else{
				my $feature_id = get_feature_id($gene_name, $PATTERN_TRANSCRIPT);
				if($feature_id ne ""){
					$transcript_id = $type;
				}
			}
		} # end unless(/^#/){
	  }# end while(<IN>){
	  close(IN);
  }else{
	if($? != 0){
			Exception::throw("An error has occurred when the script $0 tried to create/open the file $gff_file. Please check the messages returned by the operating system:", Exception->EXECUTION_ERROR, $?);
		}else{
			Exception::throw("An error has occurred when the script $0 tried to create/open the file $gff_file", Exception->UNKNOWN_FILE_WRITING, $gff_file);
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