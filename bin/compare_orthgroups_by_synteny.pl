#!/usr/bin/env perl

# GFF e Orthofinder
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
	# development constants, no need to change
	PRINT_DEV => 0,

	# code constants
	CODE=>'c',
	EQUAL => 0,
	EQUAL_BASE_ISOFORM => 1,
	EQUAL_ORTH_ISOFORM => 2,

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
	GENE_END => 'e', # this variable is the end key of the hash with the genes
	GENE_STRAND => 't', # this variable is the strand key of the hash with the genes
	MAX_ERRORS => 2, # the number maximum of errors that we accept without breaking the block
	MAX_WINDOW => 5, # the number of genes in the sliding window
	MIN_SCAFFOLD_LENGTH => 10, # a ju vai estudar esse tamanho mínimo aceitável
	GENE_BASE => 0,
	ISOFORM => 1,
	INDEX => 'i',
	TYPE => 't',

	# score values:
	SCORE_NO_UNIQ => 1,
	SCORE_NO_UNIQ_TEXT => "NO_UNIQ",
	SCORE_NO_SYNTENY => 2,
	SCORE_NO_SYNTENY_TEXT => "NO_SINTENY",
	SCORE_SCAFFOLD_MIN => 2,
	SCORE_SCAFFOLD_MIN_TEXT => "SCAFFOLD_MIN",
	SCORE_NO_OG => 2,
	SCORE_NO_OG_TEXT => "NO_OG",
	SCORE_WRONG_ISOFORM => 2,
	SCORE_WRONG_ISOFORM_TEXT => "WRONG_ISOFORM",
	SCORE_NO_BLOCK_SYNTENY => 2,
	SCORE_NO_BLOCK_SYNTENY_TEXT => "NO_SYNTENY_BLOCK",
	SCORE_LIGHT_SINTENY_ERROR => 3,
	SCORE_LIGHT_SINTENY_ERROR_TEXT => "LIGHT_SINTENY_ERROR",
	SCORE_CONSERVED_BLOCK => 4,
	SCORE_CONSERVED_BLOCK_TEXT => "CONSERVER_BLOCK",
	SCORE_EQUAL => 5,
	SCORE_EQUAL_TEXT => "EQUAL"

};

#########################################################################
## Global variables to use in this script
#########################################################################

my $file_orthofinder;    # input file from the orthofinder
my @array_input_gff; # each item in this array is one GFF file
my $output_file; # the output file with the results
my $verbose;
my $accepted_features;
my $help;

my %hash_features;
my %hash_orthogroup_by_gene = (); # chave: gene_id; valor: grupo de ortologo
my %hash_genes_by_orthogroup = (); # chave: orthologous group; valor: array de genes

my %hash_orthogroup_by_gene_uniq = (); # key 1: specie index; key 2: reference gene_id; value: ortholog group

my @array_file_in_memory = ();

my @array_scaffolds_by_species = (); # array with hashes with key: scaffold; value: array of genes
my @array_scaffolds_order_by_species = ();
my @array_isoforms_by_species = ();

my %hash_score_genes = ();

my $pattern = "\^.*ID=(transcript:)?([^;]+);?.*\$"; # pattern to extract the gene id
my $pattern_parent_name = "\^.*Parent=([^;]+);?.*\$"; # pattern to extract the gene id

###########################################################################
#receiving parameters
###########################################################################
GetOptions(
	'i=s' => \$file_orthofinder,
	'g=s'    => \@array_input_gff,
	'o=s'    => \$output_file,
	'a=s'    => \$accepted_features,
	'v|verbose' => \$verbose,
	'h|help' => \$help
); $| = 1;


main();

##################################################################################
# The script is divided into subroutines. The main() subroutine call the
# another subroutines. Tha main subroutine is just to organize our script
##################################################################################
sub main {
	if(defined $help){
		print Util::usage($0, "-i orthofinder_file -g (gff_file1 gff_file2...gff_fileN) -a accepted_features -o output_file [-v|verbose -h|help]");
		die "\n";
	}

	my $error = Checker::check_argument(
		Argument->new(
			name     => "-i",
			value     => $file_orthofinder,
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
		),
		Argument->new(
			name     => "-g",
			value     => \@array_input_gff,
			mandatory => Constant->YES,
			type      => Constant->FILE
		)
	);

	if($error){
		print Util::usage($0, "-i orthofinder_file (-g gff_file1 -g gff_file2...-g gff_fileN) -a accepted_features -o output_file [-v|verbose -h|help]"), "\n";
		Exception::throw("Validating $0 arguments");
	}

	%hash_features = map{$_ => 1} split(/,/, $accepted_features);

	print(Util::show_running("Parsing orthofinder file...")) if($verbose);
	read_orthofinder(); # reading orthofinder input file and mRNA rating output
	print(Util::show_finished("Parsing orthofinder file...\n")) if($verbose);
	
	print(Util::show_running("Parsing GFF inputs...")) if($verbose);
	read_gff_file(); # reading GFF file and populating @array_scaffolds_by_species
	print(Util::show_finished("Parsing GFF inputs...\n")) if($verbose);

	print(Util::show_running("Sorting GFF inputs...")) if($verbose);
	sort_cds();
	print(Util::show_finished("Sorting GFF inputs...\n")) if($verbose);

	print(Util::show_running("Analysing synteny...")) if($verbose);
	compare_orthologous_by_synteny_with_isoforms();
	print(Util::show_finished("Analysing synteny...\n")) if($verbose);
	
	print(Util::show_running("Exporting report...")) if($verbose);
	export_report();
	print(Util::show_finished("Exporting report...\n")) if($verbose);
	print(Util::show_finished("The script execution has been finished.\n")) if($verbose);
}

sub export_report{
	my $hash_scaffolds_reference = $array_scaffolds_by_species[0];

	open(OUT, sprintf(">%s", $output_file));

	my $is_fasta = 0;
	for(my $idx = 0; $idx < scalar @array_file_in_memory; $idx++){
		my $line = $array_file_in_memory[$idx];
		if($line =~ /^#/){
			# print 
			print OUT $line, "\n";
		}elsif($line =~ />/ || $is_fasta == 1){
			$is_fasta = 1;
			# print 
			print OUT $line, "\n";
		}else{
			my @cols = split(/\t/, $line);
			my $feature = $cols[FEATURE];
			
			if (exists $hash_features{$feature}){
				my $attributes = $cols[-1]; # The gene name is at last column
				my $scaffold = $cols[0]; # The scaffold name is the first column

				my $gene_id = get_gene_id($attributes);
				my $parent_id = get_parent_gene($attributes);
				if(exists $hash_score_genes{1}{$gene_id}){
					my $score = $hash_score_genes{1}{$gene_id}{"SCORE"};
					$line = sprintf("%s;score=%d", $line, $score);
					for(my $sp_idx = 2; $sp_idx < scalar @array_input_gff; $sp_idx++){
						my $score = $hash_score_genes{$sp_idx}{$gene_id}{"SCORE"};
						$line = sprintf("%s;%d", $line, $score);
					}
				}
			}
			print OUT $line, "\n";
		}

	}

	close(OUT);
}

sub change_gene_base_by_isoform{
	my $species_number = shift;
	my $scaffold_name = shift;
	my $gene_base_index = shift;
	my $gene_isoform = shift;
	
	my $hash_gff_genes = $array_scaffolds_by_species[$species_number];
	my $array_genes_by_sccafold = $hash_gff_genes->{$scaffold_name};

	my $hash_isoforms_by_scaffolds = $array_isoforms_by_species[$species_number];

	my $gene_base = $array_genes_by_sccafold->[$gene_base_index]->{"ID"};

	$hash_isoforms_by_scaffolds->{$gene_isoform} = delete($hash_isoforms_by_scaffolds->{$gene_base});
 
	$array_genes_by_sccafold->[$gene_base_index]->{"ID"}=$gene_isoform;
}

sub compare_orthologous_by_synteny_with_isoforms {
	for(my $species_number = 1; $species_number < scalar @array_scaffolds_by_species; $species_number++){
		my $hash_gff_reference = $array_scaffolds_by_species[0];

		my @array_scaffolds_names_orth = keys %{$array_scaffolds_by_species[$species_number]};

		my @array_scaffolds_gff_ref = keys %{$hash_gff_reference};

		my $hash_isoforms_by_scaffolds_ref = $array_isoforms_by_species[0];
		
		my $array_scaffolds_order = $array_scaffolds_order_by_species[0];

		for(@$array_scaffolds_order){
			my $scaffold_reference = $_;
			print("-----------------------------------------------------------------------------\n") if(PRINT_DEV);
			print("Scaffold $scaffold_reference\n") if(PRINT_DEV);
			my $array_genes_by_sccafold_ref = $hash_gff_reference->{$scaffold_reference};

			if(scalar @$array_genes_by_sccafold_ref > MIN_SCAFFOLD_LENGTH){
				# my $hash_uniq_gene = get_first_uniq_gene($species_number, $scaffold_referenc, $array_genes_by_sccafold_ref);
				# my $idx_uniq_gene = $hash_uniq_gene->{INDEX};
				my $idx_uniq_gene = get_first_uniq_gene($species_number, $scaffold_reference, $array_genes_by_sccafold_ref);
				
				if($idx_uniq_gene == -1){
					set_array_score_genes_with_isoform($species_number, $array_genes_by_sccafold_ref, SCORE_NO_UNIQ, SCORE_NO_UNIQ_TEXT);
					next;
				}else{
					my @initial_genes = ();

					if($idx_uniq_gene <= 1){
						@initial_genes = ($array_genes_by_sccafold_ref->[0]);
					}else{
						@initial_genes = @$array_genes_by_sccafold_ref[0..($idx_uniq_gene-1)];
					}
					
					set_array_score_genes_with_isoform($species_number, \@initial_genes, SCORE_NO_UNIQ, SCORE_NO_UNIQ_TEXT);
				}
				
				my $hash_gene_ref_uniq = $array_genes_by_sccafold_ref->[$idx_uniq_gene];
				
				my $orth_group = $hash_orthogroup_by_gene_uniq{$species_number}{$hash_gene_ref_uniq->{"ID"}};

				my $array_genes = $hash_genes_by_orthogroup{$orth_group}{$species_number};
				
				my $gene_specie_ref = $array_genes->[0]; # it is the same value as $uniq_ref_gene

				my $gene_specie_orth = $array_genes->[1]; # name of the ortholog

				# ortholog_tuple = ( scaffold_id, gene_id )
				my $ortholog_tuple = get_ortholog_tuple_with_isoforms($species_number, $gene_specie_orth);

				my $scaffold_name = $array_scaffolds_names_orth[$ortholog_tuple->[0]];
				
				my $hash_gff_n = $array_scaffolds_by_species[$species_number];

				my $array_genes_by_sccafold_orth = $hash_gff_n->{$scaffold_name};

				my $number_genes_reference = scalar @$array_genes_by_sccafold_ref;

				my @array_queue_reference = @{$array_genes_by_sccafold_ref}[($idx_uniq_gene)..($number_genes_reference-1)];
				
				# numero de genes do scaffold no homologo
				my $number_genes_orth = scalar @$array_genes_by_sccafold_orth;

				# montamos uma pilha com os genes no scaffold do homologo
				my @array_queue_orthologs = @{$array_genes_by_sccafold_orth}[($ortholog_tuple->[1])..($number_genes_orth-1)];

				my $minor_queue = ($number_genes_reference < $number_genes_orth) ? $number_genes_reference : $number_genes_orth;

				my $flag_pop_reference = 1;
				my $flag_pop_ortholog = 1;
				my $gene_reference;
				my $gene_ortholog;

				# ideia inicial colocar num hash todos os genes com erros e IGUAIS
				# aqueles que não estiverem no hash, provavelmente sao genes que nao
				# quebram o bloco
				my %hash_report = ();
				
				while($minor_queue > 0){
					
					if($flag_pop_reference){
						$gene_reference = shift(@array_queue_reference);
					}
					if($flag_pop_ortholog){
						$gene_ortholog = shift(@array_queue_orthologs);
					}

					print(sprintf("%s\t%s\n", $gene_reference->{"ID"}, $gene_ortholog->{"ID"})) if(PRINT_DEV);
					
					my $hash_status = compare_genes_with_isoforms($species_number, $gene_reference, $gene_ortholog, 1);
					if($hash_status->{POP_BASE} == 1 && $hash_status->{POP_ORTHOLOG} == 1){
						$flag_pop_reference = 1;
						$flag_pop_ortholog = 1;
					}elsif($hash_status->{POP_BASE} == 1){
						$flag_pop_reference = 1;
						$flag_pop_ortholog = 0;
					}elsif($hash_status->{POP_ORTHOLOG} == 1){
						$flag_pop_ortholog = 1;
						$flag_pop_reference = 0;
					}else{
						if(not $hash_status->{EQUAL} ){
							# aqui comparamos o gene base atual com a janela dos cinco genes subsequentes do ortólogo
							my @array_orth_window = ($gene_ortholog);
							for(my $idx = 0; ($idx < scalar @array_queue_orthologs) && ($idx < scalar @array_queue_reference && scalar @array_orth_window < MAX_WINDOW);$idx++){
								push(@array_orth_window, $array_queue_orthologs[$idx]);
							}
							my $base_to_orth = compare_from_one_to_sec_with_isoforms($gene_reference, \@array_orth_window, $species_number, 1);
							
							# aqui comparamos o gene atual com a janela dos cinco genes subsequentes do ortólogo
							my @array_base_window = ($gene_reference);
							for(my $idx = 0; ($idx < scalar @array_queue_reference) && ($idx < scalar @array_queue_orthologs && scalar @array_base_window < MAX_WINDOW);$idx++){
								push(@array_base_window, $array_queue_reference[$idx]);
							}
							my $orth_to_base = compare_from_one_to_sec_with_isoforms($gene_ortholog, \@array_base_window, $species_number, 0);

							if(!$base_to_orth && !$orth_to_base){
								set_score_genes($species_number, $gene_reference->{"ID"}, SCORE_NO_SYNTENY, SCORE_NO_SYNTENY_TEXT);
								$flag_pop_reference=1;
								$flag_pop_ortholog=1;
							}elsif(!$base_to_orth){
								set_score_genes($species_number, $gene_reference->{"ID"}, SCORE_NO_SYNTENY, SCORE_NO_SYNTENY_TEXT);
								$flag_pop_reference=1;
								$flag_pop_ortholog=0;
							}elsif(!$orth_to_base){
								$flag_pop_reference=0;
								$flag_pop_ortholog=1;
							}else{
								my $quant_equal = compare_two_lists_with_isoforms(\@array_base_window, \@array_orth_window, $species_number);
								if($quant_equal == scalar @array_base_window){
									set_array_score_genes_with_isoform($species_number, \@array_base_window, SCORE_CONSERVED_BLOCK, SCORE_CONSERVED_BLOCK_TEXT);
								}elsif($quant_equal == scalar @array_base_window - 1 ){
									set_array_score_genes_with_isoform($species_number, \@array_base_window, SCORE_LIGHT_SINTENY_ERROR, SCORE_LIGHT_SINTENY_ERROR_TEXT);
								}else{
									# pontuar todo mundo com 2?
									set_array_score_genes_with_isoform($species_number, \@array_base_window, SCORE_NO_BLOCK_SYNTENY, SCORE_NO_BLOCK_SYNTENY_TEXT);
								}
								for(my $idx = 0; $idx < scalar @array_base_window; $idx++){
									shift(@array_queue_reference);
									shift(@array_queue_orthologs);
									$flag_pop_reference=1;
									$flag_pop_ortholog=1;
								}
							}   
						}else{
							print("It should never be here\n");
						}

					}

					$number_genes_reference = scalar @array_queue_reference;
					$number_genes_orth = scalar @array_queue_orthologs;
					$minor_queue = ($number_genes_reference < $number_genes_orth) ? $number_genes_reference : $number_genes_orth;
				} # while($minor_queue > 0){
				if(scalar @array_queue_reference > 0){
					set_array_score_genes_with_isoform($species_number, \@array_queue_reference, SCORE_NO_SYNTENY, SCORE_NO_SYNTENY_TEXT);
				}

			} # if(scalar @$array_genes_by_sccafold_ref > MIN_SCAFFOLD_LENGTH){
			else{
				set_array_score_genes_with_isoform($species_number, $array_genes_by_sccafold_ref, SCORE_SCAFFOLD_MIN, SCORE_SCAFFOLD_MIN_TEXT);
				print("$scaffold_reference has too few genes\n") if(PRINT_DEV);
			}
		} # for(@array_scaffolds_gff_ref){

	} #for(my $species_number = 1; $species_number < scalar @array_scaffolds_by_species; $species_number++){

}


sub set_array_score_genes_id_with_isoform {
	my $species_number = shift;
	my $array_genes_id = shift;
	my $score = shift;
	my $category = shift;
	for(@$array_genes_id){
		set_score_genes($species_number, $_, $score, $category);
	}
}

sub set_array_score_genes_with_isoform {
	my $species_number = shift;
	my $array_genes_id = shift;
	my $score = shift;
	my $category = shift;
	for(@$array_genes_id){
		set_score_genes($species_number, $_->{"ID"}, $score, $category);
	}
}

sub set_score_genes{
	my $species_number = shift;
	my $gene_base_id = shift;
	my $score = shift;
	my $category = shift;
	my $ortholog_id = shift;
	$hash_score_genes{$species_number}{$gene_base_id}{"SCORE"} = $score;
	$hash_score_genes{$species_number}{$gene_base_id}{"CATEGORY"} = $category;

	if(defined $ortholog_id){
		$hash_score_genes{$species_number}{$gene_base_id}{"GENE_HOMO"} = $ortholog_id;
	}

}

# if there is ortholog in at least one gene or isoforms
# for both reference and homologous
sub there_is_ortholog{
	my $gene_ref =  shift;
	my $gene_orth = shift;
	my $orth_index = shift;
	my $flag_orth_ref_exists = exists $hash_orthogroup_by_gene{0}{$gene_ref};
	my $flag_orth_ort_exists = exists $hash_orthogroup_by_gene{$orth_index}{$gene_orth};

	if($flag_orth_ref_exists && $flag_orth_ort_exists){
		return 1;
	} elsif(!$flag_orth_ref_exists && !$flag_orth_ort_exists){
		my $array_isoforms_ref = $array_isoforms_by_species[0]->{$gene_ref};
		for my $isoform_name ($array_isoforms_ref){
			if(exists $hash_orthogroup_by_gene{0}{$isoform_name}){
				$flag_orth_ref_exists = 1;
				last;
			}
		}
		my $array_isoforms_orth = $array_isoforms_by_species[$orth_index]->{$gene_orth};
		for my $isoform_name ($array_isoforms_orth){
			if(exists $hash_orthogroup_by_gene{$orth_index}{$isoform_name}){
				$flag_orth_ort_exists = 1;
				last;
			}
		}
	} elsif(!$flag_orth_ref_exists){
		my $array_isoforms_ref = $array_isoforms_by_species[0]->{$gene_ref};
		for my $isoform_name ($array_isoforms_ref){
			if(exists $hash_orthogroup_by_gene{0}{$isoform_name}){
				$flag_orth_ref_exists = 1;
				last;
			}
		}
	}else{
		my $array_isoforms_orth = $array_isoforms_by_species[$orth_index]->{$gene_orth};
		for my $isoform_name ($array_isoforms_orth){
			if(exists $hash_orthogroup_by_gene{$orth_index}{$isoform_name}){
				$flag_orth_ort_exists = 1;
				last;
			}
		}
	}
	return ($flag_orth_ref_exists && $flag_orth_ort_exists);
}


sub compare_two_lists_with_isoforms{
	my $array_queue_base = shift;
	my $array_queue_ortholog = shift;
	my $species_number = shift;

	my %hash_gh1 = ();
	my %hash_gh2 = ();
	my $is_equal = 0;

	my $hash_result;
	my %hash_idx_equal = ();
	my $quant_equal = 0;
	for(my $i = 0; $i<scalar @$array_queue_base; $i++){
		my $first_gene = $array_queue_base->[$i];
		for(my $j = 0; $j<scalar @$array_queue_ortholog; $j++){
			# TODO Tem algo de errado nesta lógica
			# Preciso de casos onde entre para testar melhor
			if(not exists $hash_idx_equal{$j}){
				my $second_gene = $array_queue_ortholog->[$j];
				$hash_result = compare_genes_with_isoforms($species_number, $first_gene, $second_gene, 0);
				if($hash_result->{EQUAL} == 1){
					$hash_idx_equal{$j} = 1;
					$quant_equal++;
				}
			}
		}
	}
	
	return $quant_equal;

	
}

sub compare_two_lists{
	my $array_queue_one = shift;
	my $array_queue_second = shift;
	my $species1 = shift;
	my $species2 = shift;

	my %hash_gh1 = ();
	my %hash_gh2 = ();

	for(my $i = 0; $i<scalar @$array_queue_second; $i++){
		my $gene = $array_queue_second->[$i];

		my $orthogroup_second = $hash_orthogroup_by_gene{$species2}{$gene->{"ID"}};
		$hash_gh2{$orthogroup_second}{$i}= $i;
	}

	my %hash_result = ("status"=>0);
	my @array_equals = ();
	for(my $i = 0; $i<scalar @$array_queue_one; $i++){
		my $gene = $array_queue_one->[$i];
		my $orthogroup_first = $hash_orthogroup_by_gene{$species1}{$gene->{"ID"}};
		if(not exists $hash_gh2{$orthogroup_first}){
			$hash_result{"status"}=-1;
			if(exists $hash_result{"diff"}){
				$hash_result{"diff"}++;
			}else{
				$hash_result{"diff"} = 1;
			}
		}elsif(exists $hash_gh2{$orthogroup_first}{$i}){
			$hash_result{"equalTo"}{$i}=1;
		}
	}

	return \%hash_result;
}


# the idea of this function is compare two arrays of genes
# and verify if the first gene of the first queue has the ortolog group
# in the second queue
sub compare_from_one_to_sec_with_isoforms{
	my $first_gene = shift;
	my $array_queue_second = shift;
	my $species_number = shift;
	my $ref_to_orth = shift;

	my $is_equal = 0;
	for(my $i = 0; $i < scalar @$array_queue_second; $i++){
		my $second_gene = $array_queue_second->[$i];
		my $hash_result;
		if(!$ref_to_orth){
			$hash_result = compare_genes_with_isoforms($species_number, $second_gene, $first_gene, 0);
		}
		else{
			$hash_result = compare_genes_with_isoforms($species_number, $first_gene, $second_gene, 0);
		}
		
		if($hash_result->{EQUAL} == 1){
			# print("Inside compare_from_one_to_sec_with_isoforms:\n");
			# print("\tIs equal: ", $first_gene->{"ID"}, "\t", $second_gene->{"ID"}, "\n");
			$is_equal = 1;
			last;
		}
	}
	return $is_equal;
}

sub compare_genes_with_isoforms{
	my $species_number = shift;
	my $gene_base = shift;
	my $gene_ortholog = shift;
	my $set_score = shift;
	# my @array_isoforms_base;
	# my @array_isoforms_orth;


	my $hash_isoforms_base = $array_isoforms_by_species[0];
	my $hash_isoforms_orth = $array_isoforms_by_species[$species_number];
	my $array_isoforms_base = $hash_isoforms_base->{$gene_base->{"ID"}};
	my $array_isoforms_orth = $hash_isoforms_orth->{$gene_ortholog->{"ID"}};

	my $orthogroup_base;
	my $orthogroup_orth;
	my $base_with_og = 0;
	my $orth_with_og = 0;
	my $number_orth_equal = 0;
	my %hash_base =  ();
	my %hash_orth = ();
	my %hash_equals = ();

	my %hash_no_og = ();
	my %hash_gene_og = ();
	my %hash_result = ();
	
	# print("\t\tInside compare_genes_with_isoforms:\n");
	for my $gene_base_isoform (@$array_isoforms_base){
		# print("\t\t\tGene base: $gene_base_isoform\n");
		if(exists $hash_orthogroup_by_gene{0}{$gene_base_isoform}){
			$orthogroup_base = $hash_orthogroup_by_gene{0}{$gene_base_isoform};
			# print("\t\t\tOG base: $orthogroup_base\n");
			$hash_base{$orthogroup_base} = 1;
			$hash_gene_og{$gene_base_isoform}=$orthogroup_base;
		}else{
			$hash_no_og{$gene_base_isoform} = $orthogroup_base;
		}
	}
	for my $gene_orth_isoform (@$array_isoforms_orth){
		# print("\t\t\tGene orth: $gene_orth_isoform\n");
		if(exists $hash_orthogroup_by_gene{$species_number}{$gene_orth_isoform}){
			$orthogroup_orth = $hash_orthogroup_by_gene{$species_number}{$gene_orth_isoform};
			# print("\t\t\tOG orth: $orthogroup_orth\n");
			if(exists $hash_base{$orthogroup_orth}){
				if(exists $hash_equals{$orthogroup_orth}){
					$hash_equals{$orthogroup_orth}++;
				}else{
					$hash_equals{$orthogroup_orth} = 1;
				}
			}
			$hash_orth{$orthogroup_orth} = 1;
		}
	}

	my $number_og_base = scalar keys %hash_base;
	my $number_og_orth = scalar keys %hash_orth;
	my $number_og_equals = scalar keys %hash_equals;
	
	# se ambos sao maiores que zero pode ter algo em comum
	if($number_og_base == 0 && $number_og_orth == 0){
		# todos nao possuem grupo de ortologo
		if($set_score){
			set_array_score_genes_id_with_isoform($species_number, $array_isoforms_base, 2, "NO_GH");
		}
		# remover gene da pilha base
		# remover gene da pilha de ortologos
		$hash_result{POP_BASE} = 1;
		$hash_result{POP_ORTHOLOG} = 1;
		$hash_result{EQUAL} = 0;
	} elsif($number_og_base == 0){
		if($set_score){
			set_array_score_genes_id_with_isoform($species_number, $array_isoforms_base, 2, "NO_GH");
		}

		$hash_result{POP_BASE} = 1;
		$hash_result{POP_ORTHOLOG} = 0;
		$hash_result{EQUAL} = 0;
	} elsif($number_og_orth == 0){
		$hash_result{POP_BASE} = 0;
		$hash_result{POP_ORTHOLOG} = 1;
		$hash_result{EQUAL} = 0;
	} else{
		if( $number_og_equals >= 1 ){

			$hash_result{POP_BASE} = 1;
			$hash_result{POP_ORTHOLOG} = 1;
			$hash_result{EQUAL} = 1;
			# tem pelo menos um og igual
			if($set_score){
				for my $gene_base_isoform (@$array_isoforms_base){
					if(exists $hash_no_og{$gene_base_isoform}){
						# pontuar como sem orthogroup
						set_score_genes($species_number, $gene_base_isoform, 2, "NO_GH");
					}else{
						my $og = $hash_gene_og{$gene_base_isoform};
						if(exists $hash_equals{$og}){
							set_score_genes($species_number, $gene_base_isoform, 5, "EQUAL");
						}else{
							set_score_genes($species_number, $gene_base_isoform, 2, "WRONG_ISOFORM");
						}
					}
				}
			}
		}else{
			$hash_result{POP_BASE} = 0;
			$hash_result{POP_ORTHOLOG} = 0;
			$hash_result{EQUAL} = 0;
		}
	}
	return \%hash_result;
}


sub get_ortholog_tuple_with_isoforms {
	my $species_number = shift;
	my $gene_name = shift;

	my $hash_gff_n = $array_scaffolds_by_species[$species_number];

	my @array_scaffolds_gff_n = keys %$hash_gff_n;
	my $idx_ortholog = -1;
	my $scaffold_idx = -1;

	my $scaffold = "";
	for(my $sc_idx = 0; $sc_idx < scalar @array_scaffolds_gff_n; $sc_idx++){
		$scaffold = $array_scaffolds_gff_n[$sc_idx];

		my $array_genes_by_sccafold_n = $hash_gff_n->{$scaffold};
		my $gene_key = "";
		for(my $idx = 0; $idx<scalar @$array_genes_by_sccafold_n; $idx++){
			my $gene = $array_genes_by_sccafold_n->[$idx];
			if($gene->{"ID"} eq $gene_name){
				$idx_ortholog = $idx;
				$scaffold_idx = $sc_idx;
			}else{
				my $hash_isoforms = $array_isoforms_by_species[$species_number];
				my $array_isoforms = $hash_isoforms->{$gene->{"ID"}};
				my $isoform_gene;
				for(@$array_isoforms){
					if($_ eq $gene_name){
						$idx_ortholog = $idx;
						$scaffold_idx = $sc_idx;
						$gene_key = $gene->{"ID"};
						last;
					}
				}
			}
			last if($idx_ortholog != -1);
		}
		if($idx_ortholog != -1){
			$array_genes_by_sccafold_n->[$idx_ortholog]->{"ID"}=$gene_name;
			
			if($gene_key ne ""){
				$array_isoforms_by_species[$species_number]->{$gene_name} =
				delete($array_isoforms_by_species[$species_number]->{$gene_key});
			}
			last;
		}
		# last if($idx_ortholog != -1);
		# g has the index of the homologous gen
	} # end @array_scaffolds_gff_n
	my @array_return = ($scaffold_idx, $idx_ortholog);
	return \@array_return;
}


sub get_ortholog_tuple {
	my $species_number = shift;
	my $gene_name = shift;
	my $hash_gff_n = $array_scaffolds_by_species[$species_number];
	# print(Dumper(%$hash_gff_n));
	my @array_scaffolds_gff_n = keys %$hash_gff_n;
	my $idx_ortholog = -1;
	my $scaffold_idx = -1;

	my $scaffold = "";
	for(my $sc_idx = 0; $sc_idx < scalar @array_scaffolds_gff_n; $sc_idx++){
		$scaffold = $array_scaffolds_gff_n[$sc_idx];

		my $array_genes_by_sccafold_n = $hash_gff_n->{$scaffold};

		for(my $idx = 0; $idx<scalar @$array_genes_by_sccafold_n; $idx++){
			my $gene = $array_genes_by_sccafold_n->[$idx];
			# print($gene->{"ID"}, "\t");
			# print($gene_name, "\n");
			# print Data::Dumper::qquote($second_gene);
			# print Data::Dumper::qquote($gene->{"ID"});
			if($gene->{"ID"} eq $gene_name){
				$idx_ortholog = $idx;
				$scaffold_idx = $sc_idx;
				last;
			}
		}
		# g has the index of the homologous gen
	} # end @array_scaffolds_gff_n
	my @array_return = ($scaffold_idx, $idx_ortholog);
	return \@array_return;
}

sub get_first_uniq_gene {
	my $species_number = shift;
	my $scaffold_name = shift;
	my $array_genes = shift;

	my $idx_uniq_reference = 0;
	my $flag_no_og = 0;


	# here we are going to find the first gene who belongs to a
	# orthogroup with only uniq genes
	#not exists $hash_orthogroup_by_gene_uniq{$species_number}{$nth_gene->{"ID"}}
	
	my $gene_is_uniq = 0;
	my $isoform_is_uniq = 0;
	my $limit_no_uniq = scalar @$array_genes-2;
	my $nth_gene;
	
	do{
		$nth_gene = $array_genes->[$idx_uniq_reference++];

		$gene_is_uniq = gene_is_uniq($species_number, $nth_gene->{"ID"});
		if(!$gene_is_uniq){
			$isoform_is_uniq = isoform_is_uniq($species_number, $nth_gene->{"ID"});	
		}
	} while(!$gene_is_uniq && !$isoform_is_uniq && $idx_uniq_reference <= $limit_no_uniq);

	if(!$gene_is_uniq && !$isoform_is_uniq){
		$idx_uniq_reference = -1;
	}else{
		$idx_uniq_reference = $idx_uniq_reference-1;
		if($isoform_is_uniq){
			change_gene_base_by_isoform(0, $scaffold_name, $idx_uniq_reference,
						$isoform_is_uniq);
		}
		# return $idx_uniq_reference-1;
	}
	return $idx_uniq_reference;
}

sub gene_is_uniq{
	my $species_number = shift;
	my $gene_id = shift;
	my $is_uniq = 0;
	if(exists $hash_orthogroup_by_gene_uniq{$species_number}{$gene_id}){
		$is_uniq = 1;
	}
	return $is_uniq;
}

sub isoform_is_uniq {
	my $species_number = shift;
	my $gene_id = shift;

	my $is_uniq = 0;
	my $hash_isoforms = $array_isoforms_by_species[0];
	my $array_isoforms_ref = $hash_isoforms->{$gene_id};
	for my $isoform_name (@$array_isoforms_ref){
		if(exists $hash_orthogroup_by_gene_uniq{$species_number}{$isoform_name}){
			$is_uniq = $isoform_name;
			last;
		}
	}
	return $is_uniq;
}


sub verify_if_error_queue_is_full{
	my $array_queue = shift;
	if(scalar @$array_queue > MAX_ERRORS){
		return 1;
	}
	return 0;

}

############################################################################
# Fazer a leitura do arquivo Orthofinder
############################################################################
sub read_orthofinder {

	# opening the input file
	if(open( IN, $file_orthofinder )){
		my $count     = 0;
		my $wrap_line = 0;
		my $line = <IN>;

		$line =~ s/[\r\n]+$//g;
		my @columns = split /\t/, $line;
		my $nspecies = scalar  @columns -1;

		while (<IN>) {
			# chomp;
			$_ =~ s/[\r\n]+$//g;
			#   0   1	2	3	4	5	6	7	...
			#Orthogroup  mRNAs
			# OG0004424       PCOAH_0006314001        PKNOH_S140255000-t35_1  PVP01_1436100.1
			my @columns = split /\t/;
			my $orthogroup = $columns[0];
			my $flag = 1;


			my @genes_per_species = @columns[1..$nspecies];

			if(defined $genes_per_species[0] && $genes_per_species[0] ne ""){
				# genes for the reference species
				my @array_genes_reference = split(",\\s*", $genes_per_species[0]);
				my $is_uniq = (scalar @array_genes_reference == 1 ) ? 1 : 0;

				for(my $i = 1; $i < scalar @genes_per_species; $i++){
					if(defined($genes_per_species[$i]) && $genes_per_species[$i] ne ""){
						my @array_genes_by_orthogroup = ($genes_per_species[0], $genes_per_species[$i]);
						$hash_genes_by_orthogroup{$orthogroup}{$i} = \@array_genes_by_orthogroup;

						my @array_genes = split(",\\s*", $genes_per_species[$i]);

						$is_uniq = ($is_uniq && scalar @array_genes == 1 ) ? 1 : 0;

						for my $gene (@array_genes_reference){
							$hash_orthogroup_by_gene{0}{$gene} = $orthogroup;
						}
						for my $gene (@array_genes){
							$hash_orthogroup_by_gene{$i}{$gene} = $orthogroup;
						}
						if($is_uniq){
							# this variable tell us if a gene from the reference
							# belongs to a ortholog group with only uniq genes
							# $i -1 refers to the specie who we are comparing to
							# we prefer to use $i-1 to keep the GFF parameter order
							$hash_orthogroup_by_gene_uniq{$i}{$genes_per_species[0]} = $orthogroup;
							# $hash_orthogroup_by_gene_uniq{$i-1}{$gene_specie_n} = $orthogroup;
						}

						}
				} # end for(my $i = 2; $i < @genes_per_species; $i++){
			} # end if(defined $genes_per_species[1] && $genes_per_species[1] ne ""){

		}
		close(IN);
	}else{
		if($? != 0){
			Exception::throw("An error has occurred when the script $0 tried to open the file $file_orthofinder. Please check the messages returned by the operating system:", Exception->EXECUTION_ERROR, $?);
		}else{
			Exception::throw("An error has occurred when the script $0 tried to open the file $file_orthofinder", Exception->UNKNOWN_FILE_READING, $file_orthofinder);
		}
	}

}

#########################################################################
# Extract the gene id based on the global variable pattern
#########################################################################
sub get_gene_id{
	my $gene_text = shift;
	my $rex = qr/$pattern/;
	if($gene_text =~ $rex){
		$gene_text =~ s/$rex/$2/g;
	}
	return $gene_text;
}

sub get_parent_gene{
	my $attribute_text = shift;
	my $regex = qr/$pattern_parent_name/;
	if($attribute_text =~ $regex){
		$attribute_text =~ s/$regex/$1/g;
	}
	return $attribute_text;
}

#########################################################################
# This subroutine  sorts the features from all genes to keep the features in
# crescent order.
#########################################################################
sub sort_cds{
	for(@array_scaffolds_by_species){
		my $hash_mrna_by_scaffolds = $_;
		my @scaffolds = keys %$hash_mrna_by_scaffolds;
		for(@scaffolds){
			my $scaffold = $_;
			my $array_mrnas = $hash_mrna_by_scaffolds->{$scaffold};
			my @array_ordered = sort {$a->{GENE_START} <=> $b->{GENE_START} || $a->{GENE_END} <=> $b->{GENE_END} }  @$array_mrnas;
			$hash_mrna_by_scaffolds->{$scaffold} = \@array_ordered;
		} # end for(@scaffolds)
	}
}# end sort_cds

#########################################################################
# Reads a GFF file and organize it in some hashes and arrays
# The first GFF in the hash is the reference
#########################################################################
sub read_gff_file{

	for(my $idx = 0; $idx < scalar @array_input_gff; $idx++){
		my $file_input_gff = $array_input_gff[$idx];
		open(IN, $file_input_gff);
		my %hash_mrna_by_scaffolds = ();
		my %hash_isoforms_by_scaffolds = ();
		my %hash_parents_by_scaffolds = ();
		my @scaffold_names = ();
		while(<IN>){
			$_ =~ s/[\r\n]+$//g;
			if($idx == 0){
				push(@array_file_in_memory, $_);
			}
			if(/^>/){ # when the GFF has fasta sequences inside we stop reading the file
				last;
			}
			unless(/^#/ || /^\s*$/){  # To avoid lines which start with '#' and empty lines
					my $line = $_;
					my @cols = split /\t/;

					my $attributes = $cols[-1]; # The gene name is at last column
					my $scaffold = $cols[0]; # The scaffold name is the first column
					my $strand = $cols [STRAND];

					my $gene_id = get_gene_id($attributes);
					my $parent_id = get_parent_gene($attributes);
					my $feature = $cols[FEATURE];

					# just in case: when the start of the gene is greater than the end position we invert
					if($cols [POS_START] > $cols [POS_END]){
						my $aux = $cols [POS_START];
						$cols [POS_START] = $cols [POS_END];
						$cols [POS_END] = $aux;
					}

					if (exists $hash_features{$feature}){
							# flag to control the isoforms
							my %hash_mrna = ();
							$hash_mrna{ID} = $gene_id;
							$hash_mrna{STRAND} = $strand;
							if($cols [POS_START] > $cols [POS_END]){
								$hash_mrna{GENE_START} = $cols [POS_END];
								$hash_mrna{GENE_END} = $cols [POS_START];
							}else{
								$hash_mrna{GENE_START} = $cols [POS_START];
								$hash_mrna{GENE_END} = $cols [POS_END];
							}

							if(!exists $hash_parents_by_scaffolds{$parent_id}){
								$hash_parents_by_scaffolds{$parent_id} = $gene_id;
								my @array_isoforms = ($gene_id);
								$hash_isoforms_by_scaffolds{$gene_id} = \@array_isoforms;

								if (not exists $hash_mrna_by_scaffolds {$scaffold}){
									my @array_mrnas = ();
									push @array_mrnas, \%hash_mrna;
									$hash_mrna_by_scaffolds {$scaffold} = \@array_mrnas;
									# just to keep the order of the scaffolds at the gff
									push(@scaffold_names, $scaffold);
								} else{
									my $array_mrnas = $hash_mrna_by_scaffolds {$scaffold};
									push @$array_mrnas, \%hash_mrna;
								}
							} else{
								my $gene_model_id = $hash_parents_by_scaffolds{$parent_id};
								my $array_isoforms = $hash_isoforms_by_scaffolds{$gene_model_id};
								push(@$array_isoforms, $gene_id);
								# if (not exists $hash_isoforms_by_scaffolds {$gene_model_id}){
								#     my @array_isoforms = {$gene_model_id, $gene_id};
								#     $hash_isoforms_by_scaffolds{$gene_model_id} = \@array_isoforms;
								# }else{
								# }
							}
					}
				}
		}
		
		push(@array_isoforms_by_species, \%hash_isoforms_by_scaffolds);
		push(@array_scaffolds_by_species, \%hash_mrna_by_scaffolds);
		push(@array_scaffolds_order_by_species, \@scaffold_names);
	}
}