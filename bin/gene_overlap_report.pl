#!/usr/bin/env perl

# GFF
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

my $script_dir = dirname($0);
$script_dir = "." if($script_dir eq "");

require( $script_dir . "/lib/Constant.pm" );
require( $script_dir . "/lib/Exception.pm" );
require( $script_dir . "/lib/Checker.pm" );
require( $script_dir . "/lib/Argument.pm" );
require( $script_dir . "/lib/Util.pm" );

# Some constants to easily change the columns index if necessary
use constant{
	CHROMOSOME => 0,
	SOURCE => 1,
	FEATURE => 2, # the column from the GFF that represents the feature (cdna, mrna, gene, etc)
	POS_START => 3, # the column that represents the start position of the feature
	POS_END => 4, # the column that represents the end position of the feature
	SCORE => 5,
	STRAND => 6, # the column from the GFF that represents the dna strand (plus or minus),
	PHASE => 7,
	ATTRIBUTES => 8,
	GENE_START => 's', # this variable is the start key of the hash with the genes
	GENE_END => 'e',# this variable is the end key of the hash with the genes
	GENE_STRAND => 't', # this variable is the strand key of the hash with the genes
	GENE_PARENT => 'p', # this variable is the strand key of the hash with the genes
};

#########################################################################
## Global variables to use in this script
#########################################################################

my $input_file;
my $output_file;
my $accepted_features;
my $help;
my $verbose;

my %hash_mrna_by_scaffolds;
my %hash_start_end_cds = ();
my %hash_pseudogene = ();


my $PATTERN_CDS = "\^.*ID=(CDS:)?([a-zA-Z0-9\-_]+(\.\\d+)?)(\.cds)?;?.*\$"; # pattern to extract the gene id
my $PATTERN_TRANSCRIPT = "\^.*ID=(transcript:)?([a-zA-Z0-9\-_\.]+);.*\$"; # pattern to extract the gene id
my $PATTERN_PARENT_TRANSCRIPT = "\^.*Parent=(transcript:)?([a-zA-Z0-9\-_\.]+);?.*\$"; # pattern to extract the gene id

my %hash_features = ();

###########################################################################
#receiving parameters
###########################################################################
GetOptions(
	'i=s'    => \$input_file,
	'o=s'    => \$output_file,
	'a=s'    => \$accepted_features,
	'v|verbose' => \$verbose,
	'h|help' => \$help
	
);$| = 1;

main();


###########################################################################
sub main {
	if(defined $help){
		print Util::usage($0, "-i gff_input -o report_output -a accepted_features [-v|verbose -h|help]");
		die "\n";
	}
	
	my $error = Checker::check_argument(
		Argument->new(
			name     => "-i",
			value     => $input_file,
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
			value     => $accepted_features,
			mandatory => Constant->YES,
			type      => Constant->STRING
		)
	);
	
	if($error){
		print Util::usage($0, "-i gff_input -o report_output -a accepted_features [-v|verbose -h|help]"), "\n";
		Exception::throw("Validating $0 arguments");
	}
	
	%hash_features = map{$_ => 1} split(/,/, $accepted_features);
	print(Util::show_running("Parsing $input_file...")) if($verbose);
	read_gff_file(); # Calling the subroutine to read the GFF file
	print(Util::show_finished("Parsing $input_file...\n")) if($verbose);
	print(Util::show_running("Sorting coding regions...")) if($verbose);
	sort_cds(); # Calling the subroutine to sort the CDS in each gene
	print(Util::show_finished("Sorting coding regions...\n")) if($verbose);
	print(Util::show_running("Reporting gene overlap...")) if($verbose);
	geting_overlap_in_gene();
	print(Util::show_finished("Reporting gene overlap...\n")) if($verbose);
	print(Util::show_finished("The script execution has been finished.\n")) if($verbose);

}

#########################################################################
# This subroutine  sorts the features from all genes to keep the features in
# crescent order.
#########################################################################
sub sort_cds{
	my @scaffolds = keys %hash_mrna_by_scaffolds;
	for(@scaffolds){
		my $scaffold = $_;
		my $array_mrnas = $hash_mrna_by_scaffolds{$scaffold};
		my @array_ordered = sort {$a->{GENE_START} <=> $b->{GENE_START} || $a->{GENE_END} <=> $b->{GENE_END} }  @$array_mrnas;
		for(@array_ordered){
			my $hash_mrna = $_;
			if(defined $hash_start_end_cds{$hash_mrna->{ID}}{START}){
				if($hash_mrna->{GENE_START} < $hash_start_end_cds{$hash_mrna->{ID}}{START}){
					$hash_mrna->{GENE_START} = $hash_start_end_cds{$hash_mrna->{ID}}{START};
				}
				if($hash_mrna->{GENE_END} > $hash_start_end_cds{$hash_mrna->{ID}}{END}){
					$hash_mrna->{GENE_END} = $hash_start_end_cds{$hash_mrna->{ID}}{"END"};
				}    
			}else{
				print($hash_mrna->{ID}, "\n");
				Exception::throw("There was a problem with the file $input_file.", Exception->LOGIC_ERROR, "There are feature(s) 'gene' without the corresponding 'CDS' feature(s). Please fix the GFF file before running the script again.");
			}
		}
		$hash_mrna_by_scaffolds{$scaffold} = \@array_ordered;
	} # end for(@scaffolds)

}# end sort_cds
#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
#########################################################################
sub read_gff_file{
# primeiro a gente abre o arquivo

  if(open(IN, $input_file)){
	while(<IN>){
		$_ =~ s/[\r\n]+$//g;
		unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#'
			my $line = $_;
			my @cols = split /\t/;

			my $gene_name = $cols[-1]; # The gene name is at last column
			my $scaffold = $cols[0]; # The scaffold name is the first column
			my $strand = $cols [STRAND];

			my $gene_id = get_feature_id($gene_name, $PATTERN_TRANSCRIPT);
			my $parent_id = get_feature_id($gene_name, $PATTERN_PARENT_TRANSCRIPT);

			my $feature = $cols[FEATURE];
					my $attributes = $cols[ATTRIBUTES];

			# just in case: when the start of the gene is greater than the end position we invert
			if($cols [POS_START] > $cols [POS_END]){
			  my $aux = $cols [POS_START];
			  $cols [POS_START] = $cols [POS_END];
			  $cols [POS_END] = $aux;
			}
			if(lc($feature) eq "pseudogene"){
				$hash_pseudogene{$gene_id} = 1;
			}elsif(exists $hash_features{$feature} && not exists $hash_pseudogene{$parent_id}){
			  my %hash_mrna = ();
			  $hash_mrna{ID} = $gene_id;
			  $hash_mrna{STRAND} = $strand;
			  $hash_mrna{GENE_START} = $cols [POS_START];
			  $hash_mrna{GENE_END} = $cols [POS_END];
			  $hash_mrna{GENE_PARENT} = $parent_id;

			  if (not exists $hash_mrna_by_scaffolds {$scaffold}){
				my @array_mrnas = ();
				push @array_mrnas, \%hash_mrna;
				$hash_mrna_by_scaffolds{$scaffold} = \@array_mrnas;
			  }else{
				my $array_mrnas = $hash_mrna_by_scaffolds {$scaffold};
				push @$array_mrnas, \%hash_mrna;
			  }
			}elsif(lc($feature) eq "cds"){
				# my $cds_id = get_feature_id($gene_name, $PATTERN_CDS);
				if(exists $hash_start_end_cds{$parent_id}){
					if($cols [POS_START] < $hash_start_end_cds{$parent_id}{"START"}){
						$hash_start_end_cds{$parent_id}{"START"} = $cols [POS_START] ;
					}
					if($cols [POS_END] > $hash_start_end_cds{$parent_id}{"END"}){
						$hash_start_end_cds{$parent_id}{"END"} = $cols [POS_END];
					}
				}else{
					$hash_start_end_cds{$parent_id}{"START"} = $cols [POS_START];
					$hash_start_end_cds{$parent_id}{"END"} = $cols [POS_END];
				}
			}
		} # end unless(/^#/){
	  }# end while(<IN>){
	  close(IN);
  }else{
	print STDERR "\n";
	if($? != 0){
		Exception::throw("An error has occurred when the script $0 tried to open the file $input_file. Please check the messages returned by the operating system:", Exception->EXECUTION_ERROR, $?);
	}else{
		Exception::throw("An error has occurred when the script $0 tried to open the file $input_file.", Exception->UNKNOWN_FILE_READING, $input_file);
	}
  }
}

#########################################################################
# Sub overlap in genes
#########################################################################
sub geting_overlap_in_gene{
	if(open(OUT, ">".$output_file)){
		print OUT "CHROMOSOME\tGENE_ID\tOVERLAP_WITH\n";
		while(my ($scaffold, $array_mrna) = each(%hash_mrna_by_scaffolds)){
			for(my $i = 0; $i<scalar @{$array_mrna}-1; $i++){
				my $gene_i = $array_mrna->[$i];
#                print(Dumper(%{$gene_i}));
				for(my $j = $i+1; $j<scalar @{$array_mrna}; $j++){
					my $gene_j = $array_mrna->[$j];
					if($gene_i->{"GENE_END"} >= $gene_j->{"GENE_START"}){
						if($gene_i->{"GENE_PARENT"} ne $gene_j->{"GENE_PARENT"}){
							print OUT sprintf("%s\t%s\t%s\n",$scaffold,$gene_i->{"ID"},$gene_j->{"ID"});
						}
					}
				}
			}#end for(my $i = 0; $i< scalar @{$array_mrna}; $i++){
		}#end while(my ($scaffold, $array_mrna) = each(%hash_mrna_by_scaffolds)){
		close(OUT);
	}else{
		print STDERR "\n"; 
		if($? != 0){
			Exception::throw("An error has occurred when the script $0 tried to create/open the file $output_file. Please check the messages returned by the operating system:", Exception->EXECUTION_ERROR, $?);
		}else{
			Exception::throw("An error has occurred when the script $0 tried to create/open the file $output_file", Exception->UNKNOWN_FILE_WRITING, $output_file);
		}
	
	}
	
}

########################################################################
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
