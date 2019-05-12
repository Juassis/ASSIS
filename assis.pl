#!/usr/bin/env perl

# use strict;
use warnings;
use Getopt::Long;
use File::Basename;

use constant(
	SHOW_LOG => 1
);

my $script_dir = dirname($0);
$script_dir = "." if($script_dir eq "");

require( $script_dir . "/bin/lib/Constant.pm" );
require( $script_dir . "/bin/lib/Exception.pm" );
require( $script_dir . "/bin/lib/Checker.pm" );
require( $script_dir . "/bin/lib/Argument.pm" );
require( $script_dir . "/bin/lib/Util.pm" );

my %hash_params;

my $gff_reference;
my $fasta_reference;
my $gff_evaluated;
my $fasta_evaluated;

# my @fasta_files;
# my @gff_files;
# my %hash_fasta = ();

my $project_dir;
my $pipeline_dir;
my $params_file;
my $help;
my $GFF_DIR;
my $FASTA_DIR;
my $OUTPUT_DIR;

my $orthofinder_folder; # name of the Orthofinder folder result

GetOptions(
	'd=s' => \$project_dir,
	'p=s' => \$params_file,
	'h|help' => \$help,
);

$| = 1;

main();

# The main steps of the pipeline are executed insite the main function below
sub main {
	if($help){
		print Util::usage($0, "[-d project_dir -p parameter_file]");
		die "\n";
	}

	if(defined $project_dir){
		$project_dir = fix_directory($project_dir);
		$params_file = (defined $params_file) ? "$params_file" : "$project_dir/params";
	}else{
		$project_dir = ".";
		$params_file = (defined $params_file) ? "$params_file" : "params";
	}

	# Before the pipeline logic we need to validate the arguments sent by the user
	my $error = Checker::check_argument(
		Argument->new(
	  name     => "-d",
			value     => $project_dir,
			mandatory => Constant->YES,
			type      => Constant->DIRECTORY_NOT_EMPTY
		),
		Argument->new(
	  name     => "-p",
			value     => $params_file,
			mandatory => Constant->NO,
			type      => Constant->FILE
		)
	);
	if($error){
			print Util::usage($0, "[-d project_dir -p parameter_file]"), "\n";
			Exception::throw("Validating $0 arguments");
	}

	$pipeline_dir = $script_dir . "/bin";

	print(Util::show_running("Reading pipeline parameters file...")) if(SHOW_LOG);
	read_params();
	print(Util::show_finished("Reading pipeline parameters file...\n"))  if(SHOW_LOG);

	print(Util::show_running("Verifying user defined pipeline directories...")) if(SHOW_LOG);
	validate_defined_directories();
	print(Util::show_finished("Verifying user defined pipeline directories...\n")) if(SHOW_LOG);

	print(Util::show_running("Reading input files...")) if(SHOW_LOG);
	read_gff_directory();
	print(Util::show_finished("Reading input files...\n")) if(SHOW_LOG);


	print(Util::show_running("Extracting and translating mRNA sequences...")) if(SHOW_LOG);
	extract_protein();
	print(Util::show_finished("Extracting and translating mRNA sequences...\n")) if(SHOW_LOG);

	print(Util::show_running("Generating protein coding sequences report...")) if(SHOW_LOG);
	generate_codon_report();
	print(Util::show_finished("Generating protein coding sequences report...\n")) if(SHOW_LOG);

	print(Util::show_running("Generating gene overlap report...")) if(SHOW_LOG);
	generate_overlap_report();
	print(Util::show_finished("Generating gene overlap report...\n")) if(SHOW_LOG);

	if(lc($hash_params{"RNASEQ"}) eq "yes"){
		print(Util::show_running("Generating RNASeq report...")) if(SHOW_LOG);
		generate_rnaseq_report();
		print(Util::show_finished("Generating RNASeq report...\n")) if(SHOW_LOG);
	}
	print(Util::show_running("Running Orthofinder. This could take a while...")) if(SHOW_LOG);
	run_orthofinder();
	print(Util::show_finished("Running Orthofinder. This could take a while...\n")) if(SHOW_LOG);

	print(Util::show_running("Calculating synteny score...")) if(SHOW_LOG);
	generate_syntery_report();
	print(Util::show_finished("Calculating synteny score...\n")) if(SHOW_LOG);

	print(Util::show_finished("The pipeline execution has been finished. Check out the output files inside the directory $OUTPUT_DIR.\n")) if(SHOW_LOG);
}

sub generate_syntery_report{
	unless ( -e "$OUTPUT_DIR/REPORT" ) {
		`mkdir -p $OUTPUT_DIR/REPORT`;
	}

	my $gff_argument = sprintf("-g %s -g %s", $gff_evaluated, $gff_reference);
#	$orthofinder_folder = "Results_Apr05";
	my $command = sprintf(
	"$pipeline_dir/compare_orthgroups_by_synteny.pl -i $OUTPUT_DIR/PROTEINS/%s/Orthogroups.csv $gff_argument -a %s -o $OUTPUT_DIR/REPORT/%s.gff",
	$orthofinder_folder, $hash_params{"ACCEPTED_TRANSCRIPTS"}, $hash_params{"EVALUATED_ANNOTATION"});
	print($command, "\n");
	my $output = `$command`;

	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}
}

sub run_orthofinder {
	unless ( -e "$OUTPUT_DIR/REPORT" ) {
		`mkdir -p $OUTPUT_DIR/REPORT`;
	}

	# create folder name
	my ( $day, $month ) = (localtime)[ 3, 4 ];
	$month = (qw/ Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec /)[$month];

	my $t_param = "";
	if(exists $hash_params{"ORTHOFINDER_PARAMETER_T"}){
		$t_param = "-t ".$hash_params{"ORTHOFINDER_PARAMETER_T"};
	}
	my $a_param = "";
	if(exists $hash_params{"ORTHOFINDER_PARAMETER_A"}){
		$a_param = "-a ".$hash_params{"ORTHOFINDER_PARAMETER_A"};
	}

	my $command .= sprintf(
		"mv %s/$OUTPUT_DIR/PROTEINS/%s.fasta %s/$OUTPUT_DIR/PROTEINS/000_%s.fasta", $project_dir, $hash_params{"EVALUATED_ANNOTATION"},
		$project_dir, $hash_params{"EVALUATED_ANNOTATION"}
	);
	print("$command\n");
	my $output = `$command`;
	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}

	$orthofinder_folder = sprintf( "Results_%s%02d", $month, $day );
	$command = "";
	if ( exists $hash_params{"CONDA_ENVIRONMENT"}
		&& defined $hash_params{"CONDA_ENVIRONMENT"} )
	{
		$command = sprintf(
			"bash -c \"source activate %s && %s -f %s/$OUTPUT_DIR/PROTEINS %s %s \"",
			$hash_params{"CONDA_ENVIRONMENT"},
			$hash_params{"ORTHOFINDER_COMMAND"},
			$project_dir, $t_param,
			$a_param
		);
	}else{
		$command = sprintf(
			"%s -f %s/$OUTPUT_DIR/PROTEINS %s %s",
			$hash_params{"ORTHOFINDER_COMMAND"},
			$project_dir, $t_param,
			$a_param
		);
	}
	print("$command\n");

	$output = `$command`;

	$command = sprintf(
	"mv %s/$OUTPUT_DIR/PROTEINS/000_%s.fasta %s/$OUTPUT_DIR/PROTEINS/%s.fasta", $project_dir, $hash_params{"EVALUATED_ANNOTATION"},
	$project_dir, $hash_params{"EVALUATED_ANNOTATION"}
	);
	$output = `$command`;
	if ($? != 0) {
			Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}
}



sub generate_rnaseq_report {
	unless ( -e "$OUTPUT_DIR/REPORT" ) {
		my $output = `mkdir -p $OUTPUT_DIR/REPORT 2>&1`;
	}
	my $command = sprintf(
		"$pipeline_dir/rnaseq_report.pl -i %s -a %s -o $OUTPUT_DIR/REPORT/report_rnaseq",
		$gff_evaluated, $hash_params{"ACCEPTED_TRANSCRIPTS"});
	my $output = `$command`;
	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}
}

sub generate_overlap_report {
	unless ( -e "$OUTPUT_DIR/REPORT" ) {
		my $output = `mkdir -p $OUTPUT_DIR/REPORT 2>&1`;
		if ($? != 0) {
			Exception::throw("It was not possible to create the directory $OUTPUT_DIR/REPORT. The pipeline will be terminated.", Exception->EXECUTION_ERROR, $output);
		}
	}
	my $command = sprintf(
		"$pipeline_dir/gene_overlap_report.pl -i %s -a %s -o $OUTPUT_DIR/REPORT/report_overlap",
		$gff_evaluated, $hash_params{"ACCEPTED_TRANSCRIPTS"} );
	my $output = `$command`;
	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}

}

sub generate_codon_report {
	unless ( -e "$OUTPUT_DIR/REPORT" ) {
		my $output = `mkdir -p $OUTPUT_DIR/REPORT 2>&1`;
		if ($? != 0) {
			Exception::throw("It was not possible to create the directory $OUTPUT_DIR/REPORT. The pipeline will be terminated.", Exception->EXECUTION_ERROR, $output);
		}
	}
	my $command = sprintf(
		"$pipeline_dir/protein_coding_report.pl -f $OUTPUT_DIR/PROTEINS/%s.fasta -o $OUTPUT_DIR/REPORT/report_protein_coding",
		$hash_params{"EVALUATED_ANNOTATION"});
	my $output = `$command`;
	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $pipeline_dir/protein_coding_report.pl script", Exception->EXECUTION_ERROR, $output);
	}

}

sub extract_protein {
	unless ( -e "$OUTPUT_DIR/PROTEINS" ) {
		my $output = `mkdir -p $OUTPUT_DIR/PROTEINS 2>&1`;
		if ($? != 0) {
			Exception::throw("It was not possible to create the directory $OUTPUT_DIR/PROTEINS. The pipeline will be terminated.", Exception->EXECUTION_ERROR, $output)
		}
	}

	my $command = sprintf(
		"$pipeline_dir/gff2fasta.pl -f %s -g %s -o $OUTPUT_DIR/PROTEINS/%s.fasta",
		$fasta_evaluated, $gff_evaluated, $hash_params{"EVALUATED_ANNOTATION"});

	my $output = `$command 2>&1`;

	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}

	$command = sprintf(
		"$pipeline_dir/gff2fasta.pl -f %s -g %s -o $OUTPUT_DIR/PROTEINS/%s.fasta",
		$fasta_reference, $gff_reference, $hash_params{"REFERENCE_ANNOTATION"});

	$output = `$command 2>&1`;

	if ($? != 0) {
		Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output);
	}

}

sub read_gff_directory {
	my $expression = $GFF_DIR."/".$hash_params{"EVALUATED_ANNOTATION"}."*";
	my @files = glob $expression;
	if(scalar @files == 1){
		$gff_evaluated   = $files[0];
	}else{
		Exception::throw("Insige $GFF_DIR directory there is zero or more than one file with the name ".$hash_params{"EVALUATED_ANNOTATION"},
			Exception->EXECUTION_ERROR, "Please check the file $params_file to check if it's configured correctly.");
	}
	$expression = $GFF_DIR."/".$hash_params{"REFERENCE_ANNOTATION"}."*";

	@files = glob $expression;
	if(scalar @files == 1){
		$gff_reference   = $files[0];
	}else{
		Exception::throw("Insige $GFF_DIR directory there is zero or more than one file with the name ".$hash_params{"REFERENCE_ANNOTATION"},
			Exception->EXECUTION_ERROR, "Please check the file $params_file to check if it's configured correctly.");
	}

	$expression = $FASTA_DIR."/".$hash_params{"EVALUATED_ANNOTATION"}."*";
	@files = glob $expression;
	@files = grep !/\.index$/, @files;

	if(scalar @files == 1){
		$fasta_evaluated   = $files[0];
	}else{
		Exception::throw("Insige $FASTA_DIR directory there is zero or more than one file with the name ".$hash_params{"EVALUATED_ANNOTATION"},
			Exception->EXECUTION_ERROR, "Please check the file $params_file to check if it's configured correctly.");
	}

	$expression = $FASTA_DIR."/".$hash_params{"REFERENCE_ANNOTATION"}."*";
	@files = glob $expression;
	@files = grep !/\.index$/, @files;
	if(scalar @files == 1 ){
		$fasta_reference   = $files[0];
	}else{
		Exception::throw("Insige $FASTA_DIR directory there is zero or more than one file with the name ".$hash_params{"REFERENCE_ANNOTATION"},
			Exception->EXECUTION_ERROR, "Please check the file $params_file to check if it's configured correctly.");
	}
}

sub validate_defined_directories{
#    $fixed_dir="";
	my $fixed_dir =  $project_dir ne "." ? $project_dir."/" : "";
	$GFF_DIR    = $fixed_dir.fix_directory( $hash_params{"GFF_DIR"} );
	$FASTA_DIR  = $fixed_dir.fix_directory( $hash_params{"FASTA_DIR"} );
	$OUTPUT_DIR = $fixed_dir.fix_directory( $hash_params{"OUTPUT_DIR"} );
	if(-d $OUTPUT_DIR){
		Exception::throw("Something went wrong with the user defined GFF directory.") if(Checker::check_directory($GFF_DIR, Constant->NO));
	}else{
		my $output = `mkdir $OUTPUT_DIR 2>&1`;
		if ($? != 0) {
			Exception::throw("Finishing pipeline execution due to an error at $0 script", Exception->EXECUTION_ERROR, $output)
		}
	}
	Exception::throw("Something went wrong with the user defined FASTA directory.") if(Checker::check_directory($FASTA_DIR, Constant->NO));
	Exception::throw("Something went wrong with the user defined OUTPUT directory.") if(Checker::check_directory($OUTPUT_DIR, Constant->YES));
}

sub read_params {
	Exception::throw("Reading the pipeline parameter file.") if(Checker::check_file($params_file));
	if ( open( IN, $params_file ) ) {
		while (<IN>) {
			$_ =~ s/[\r\n]+$//g;
			unless ( /^#/ || /^\s*$/ ) {
				my @split_name = split( /\s*=\s*/, $_ );
				$hash_params{ $split_name[0] } = $split_name[1];
			}
		}
		close(IN);
	}
	else {
		Exception::throw("The script could not open the file.", Exception->UNKNOWN_FILE_READING, $params_file);
	}
}

sub fix_directory {
	my $dir = shift;
	$dir =~ s/\/$//g;
	return $dir;
}
