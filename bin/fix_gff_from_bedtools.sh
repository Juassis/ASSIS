#!/usr/bin/env sh

FILE=$1

if [ -f $FILE ]
then
	cat $FILE | perl -n -e 'unless(/^#/){ chomp; @cols = split /\t/; print(join("\t", @cols[0..8]), ";bedtools=",$cols[12],"\n");}else{print $_;} '
else
	echo "File does not exist\n"
fi


#bedtools coverage -a file.bam -b file.gff > file.bedtools
