package Checker;

use strict;
use warnings;
use Data::Dumper;

# import FileChecker;
# import Argument;


sub check_argument {
	# a list of Argument objects
	my @params = @_;
	my $checker;
    
	for my $argument (@params){
        $checker = check_mandatory($argument->mandatory, $argument->name, $argument->value);
        last if($checker);
        if(defined $argument->value || $argument->mandatory == Constant->YES){
            if($argument->type == Constant->FILE){
                if(ref($argument->value) eq 'ARRAY') {
                    for(@{$argument->value}){
                        $checker = check_file($_);
                        last if($checker);
                    }
                }else{
                    $checker = check_file($argument->value);
                }
            }elsif($argument->type == Constant->DIRECTORY){
                $checker = check_directory($argument->value, Constant->YES);
            }elsif($argument->type == Constant->DIRECTORY_NOT_EMPTY){
                $checker = check_directory($argument->value, Constant->NO);
            }elsif($argument->type == Constant->FILE_TO_CREATE){
                $checker = check_output_file($argument->value);
            }
            last if($checker);
        }
	}
	return $checker;
}

sub check_mandatory{
    my $mandatory = shift;
    my $name = shift;
    my $value = shift;
    my $exceptionType;
    if(ref($value) eq 'ARRAY') {
        if($mandatory == Constant->YES){
            for(@{$value}){
                if(not defined $_){
                    $exceptionType = Exception->MISSING_PARAMETER;
                    Exception::register_exception($exceptionType, $name);
                    last;
                }
            }
        }
    }else{
        if($mandatory == Constant->YES && not defined $value){
            $exceptionType = Exception->MISSING_PARAMETER;
            Exception::register_exception($exceptionType, $name);
        }
    }
    
    return($exceptionType);
}

sub file_name_is_valid{
    my $file = shift;
	if($file =~ /\/$/){
		return 0;
	}
    return 1;
}

sub check_output_file{
    my $file = shift;
	my $exceptionType;
	if(!file_name_is_valid($file)){
		$exceptionType = Exception->FILE_NAME_OUTPUT_INVALID;
        Exception::register_exception($exceptionType, $file);
	}else{
        if($file =~ /^(.+[\/\\])+.+$/){
            my $dir = $1; chop($dir);
            $exceptionType = check_directory($dir, Constant->YES);
            Exception::register_exception($exceptionType, $dir);
        }
    }
	return($exceptionType);
}

sub check_file {
	my $file = shift;
	my $exceptionType;
	if(!file_name_is_valid($file)){
        $exceptionType = Exception->FILE_NAME_INPUT_INVALID;
    }elsif(!-e $file){
		$exceptionType = Exception->FILE_NOT_FOUND;
	}elsif(!-f $file){
		$exceptionType = Exception->FILE_INVALID;
	}elsif(!-r $file){
		$exceptionType = Exception->FILE_NOT_READABLE;
	}elsif(-z $file){
		$exceptionType = Exception->FILE_EMPTY;
	}
    
	if($exceptionType){
		my ($package, $filename, $line) = caller;
		Exception::register_exception($exceptionType, $file);
	}
	return($exceptionType);
}

sub check_directory {
	my $directory = shift;
    my $allow_empty = shift;
	my $exceptionType;
	if(!-e $directory){
		$exceptionType = Exception->DIRECTORY_NOT_FOUND;
	}elsif(!-d $directory){
		$exceptionType = Exception->DIRECTORY_INVALID;
	}elsif(!-r $directory){
		$exceptionType = Exception->DIRECTORY_NOT_READABLE;
	}elsif($allow_empty == Constant->NO && is_dir_empty($directory)){
        $exceptionType = Exception->DIRECTORY_EMPTY;
    }
    
	if($exceptionType){
		Exception::register_exception($exceptionType, $directory);
	}
	return($exceptionType);
}

sub is_dir_empty {
    my ($dir) = @_;

    opendir my $h, $dir
        or die "Cannot open directory: '$dir': $!";

    while ( defined (my $entry = readdir $h) ) {
        return unless $entry =~ /^[.][.]?\z/;
    }

    return 1;
}

1;
