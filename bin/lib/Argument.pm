package Argument;

# my $script_dir = $0; $script_dir =~ s/\w+.pl//g;
# require($script_dir."bin/lib/Constants.pm");

use strict;
use warnings;
use Data::Dumper;


sub new {
  my ($class, %args) = @_;
  return bless \%args, $class;
}

sub name{
  my $self = shift;
  my $param = shift;
  if(defined $param){
    $self->{"name"}=$param;
  }else{
    return $self->{"name"};
  }
}

sub value{
  my $self = shift;
  my $param = shift;
  if(defined $param){
    $self->{"value"}=$param;
  }else{
    return $self->{"value"};
  }
}

sub mandatory{
  my $self = shift;
  my $value = shift;
  if(defined $value){
    $self->{"mandatory"}=$value;
  }else{
    return $self->{"mandatory"};
  }
}

sub type{
  my $self = shift;
  my $value = shift;
  if(defined $value){
    $self->{"type"}=$value;
  }else{
    return $self->{"type"};
  }
}



1;
