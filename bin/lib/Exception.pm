package Exception;

# my $script_dir = $0; $script_dir =~ s/\w+.pl//g;
# require($script_dir."bin/lib/Constants.pm");

use strict;
use warnings;
use Data::Dumper;
use base 'Exporter';
use Term::ANSIColor;

use constant{
	FILE_NOT_FOUND => 1,
    FILE_NAME_OUTPUT_INVALID => 2,
    FILE_NAME_INPUT_INVALID => 3,
	FILE_INVALID => 4,
	FILE_NOT_READABLE => 5,
	FILE_EMPTY => 6,
	DIRECTORY_NOT_FOUND=>7,
	DIRECTORY_INVALID=>8,
	DIRECTORY_NOT_READABLE=>9,
	DIRECTORY_EMPTY=>10,
    UNKNOWN_FILE_READING=>11,
    UNKNOWN_FILE_WRITING=>12,
    UNKNOWN_DIRECTORY_READING=>13,
    MISSING_PARAMETER=>14,
    UNKNOWN_FILE_CREATION=>15,
    UNKNOWN_DIRECTORY_CREATION=>16,
    EXECUTION_ERROR=>17,
    LOGIC_ERROR =>18
};

our @EXPORT_OK = ('FILE_NOT_FOUND', 'FILE_NAME_OUTPUT_INVALID','FILE_NAME_INPUT_INVALID',              'FILE_INVALID', 'FILE_NOT_READABLE',
	'FILE_EMPTY', 'DIRECTORY_NOT_FOUND', 'DIRECTORY_INVALID', 'DIRECTORY_NOT_READABLE',
    'DIRECTORY_EMPTY', 'UNKNOWN_FILE_READING','UNKNOWN_FILE_WRITING', 'UNKNOWN_DIRECTORY_READING','MISSING_PARAMETER',
    'UNKNOWN_FILE_CREATION','UNKNOWN_DIRECTORY_CREATION','EXECUTION_ERROR', 'LOGIC_ERROR');

sub throw{
	my $aditionalMessage = shift;
    my $exceptionType = shift;
	my $exceptionValue = shift;

	if(exists $ENV{"ASSIS_EXCEPTION"} ||
							(defined $exceptionType && defined $exceptionValue) ){
        $exceptionType = $ENV{"ASSIS_EXCEPTION"} if(not defined $exceptionType);
        $exceptionValue = $ENV{"ASSIS_EXCEPTION_VALUE"} if(not defined $exceptionValue);
		my $message;
		my $header_message = Constant->PREFIX_ERROR." ".$aditionalMessage;
		if($exceptionType == FILE_NOT_FOUND){
			$message = sprintf("%s File %s was not found. Please check if the file exists.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == FILE_NAME_OUTPUT_INVALID){
			$message = sprintf("%s %s is not a valid output file name. Please give a valid output file name.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}
        elsif($exceptionType == FILE_NAME_INPUT_INVALID){
			$message = sprintf("%s %s is not a valid input file name. Please give a valid input file name.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}
        elsif($exceptionType == FILE_INVALID){
			$message = sprintf("%s File %s is not a valid file. Please check if it's a plain text file.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == FILE_NOT_READABLE){
			$message = sprintf("%s File %s has no reading permission. Please give reading permission to the file.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == FILE_EMPTY){
			$message = sprintf("%s File %s is empty. Please check if the file is the correct one.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == DIRECTORY_NOT_FOUND){
			$message = sprintf("%s Directory %s/ was not found. Please check if the directory exists.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == DIRECTORY_INVALID){
			$message = sprintf("%s Directory %s/ is not a valid directory. Please check if it's a directory.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == DIRECTORY_NOT_READABLE){
			$message = sprintf("%s Directory %s/ has no reading permission. Please give reading permission to the directory.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == DIRECTORY_EMPTY){
			$message = sprintf("%s Directory %s/ is empty. Please check if you configured the parameter file correctly.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == UNKNOWN_FILE_READING){
			$message = sprintf("%s Sorry, but we could not read the file %s.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == UNKNOWN_FILE_WRITING){
			$message = sprintf("%s Sorry, but we could not open the file %s for writing.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}
        elsif($exceptionType == UNKNOWN_DIRECTORY_READING){
			$message = sprintf("%s Sorry, but we could not read the directory %s.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == MISSING_PARAMETER){
			$message = sprintf("%s The parameter %s is missing, although it is mandatory.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == UNKNOWN_FILE_CREATION){
			$message = sprintf("%s Sorry, but the script was not able to create the file %s.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == UNKNOWN_DIRECTORY_CREATION){
			$message = sprintf("%s Sorry, but the script was not able to create the directory %s.",
				Constant->PREFIX_ERROR, $exceptionValue);
		}elsif($exceptionType == EXECUTION_ERROR){
            my @system_message = split("\n",$exceptionValue);
            $message = "";
            for(@system_message){
                $message .= sprintf("%s %s\n",Constant->PREFIX_EXECUTION_ERROR, $_);
            }
		}elsif($exceptionType == LOGIC_ERROR){
            $message = sprintf("%s %s",Constant->PREFIX_ERROR, $exceptionValue);
        }

		if($message){
            print STDERR color('bold red');
			print STDERR $header_message;
			print STDERR  color('reset');
            print STDERR "\n";
			print STDERR $message;
			delete($ENV{"ASSIS_EXCEPTION"});
			delete($ENV{"ASSIS_EXCEPTION_VALUE"});
			delete($ENV{"ASSIS_EXCEPTION_SOURCE"});
			die "\n";
		}
	}
}

sub register_exception {
	my $exceptionType = shift;
	my $value = shift;
	$ENV{"ASSIS_EXCEPTION"}=$exceptionType;
	$ENV{"ASSIS_EXCEPTION_VALUE"}=$value;
}

1;
