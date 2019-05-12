package Util;

use strict;
use warnings;
use Term::ANSIColor;

sub usage {
	my $script_name = shift;
    my $params_text = shift;
    return sprintf("%s %s%s %s%s", Constant->PREFIX_USAGE, color('bold'), $script_name, $params_text, color('reset'))
}


sub show_command{
    my $text = shift;
    return sprintf("%s %s", Constant->PREFIX_COMMAND, $text)
}

sub show_running{
    my $text = shift;
    return sprintf("%s %s", Constant->PREFIX_RUNNING, $text)
}
sub show_finished{
    my $text = shift;
    return sprintf("\r%s %s", Constant->PREFIX_FINISHED, $text)
}

sub bold {
	my $text = shift;
    return sprintf("%s%s%s", color('bold'), $text, color('reset'));
}

1;
