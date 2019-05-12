
my $file_orthogroups = shift or die "perl $0 orthogroup_file GFF1 GFF2";
my $file_gff1 = shift;
my $file_gff2 = shift;


my %hash_genes = ();


open(IN, $file_gff1);

while(<IN>){
	next if /^#/;
	chomp;
	my @cols = split /\t/;
	my $attributes_colum = $cols[-1];
	if($cols[2] eq "mRNA"){
		my $attributes_colum = $cols[-1];
		if($attributes_colum =~ /^ID=(transcript\:)?(.+);Parent=.+$/){
			$hash_genes{$2} = 1;
		}
	}
}

close(IN);


open(IN, $file_gff2);

while(<IN>){
	next if /^#/;
	chomp;
	my @cols = split /\t/;
	my $attributes_colum = $cols[-1];
	if($cols[2] eq "mRNA"){
		my $attributes_colum = $cols[-1];
		if($attributes_colum =~ /^ID=(transcript\:)?(.+);Parent=.+$/){
			$hash_genes{$2} = 1;
		}
	}
}

close(IN);

open(IN, $file_orthogroups);
<IN>;
while(<IN>){
	next if /^#/;
	chomp;
	my $line = $_;
	my @cols = split /\t/;
	my $og = $cols[0];
	my @genes_reference = split(/, /, $cols[1]) if($cols[1] ne "");
	my @genes_ortholog = split(/, /, $cols[2]) if($cols[2] ne "");
	my $is_printed = 0;
	for(@genes_reference){
		if(exists $hash_genes{$_}){
			print $line, "\n";
			$is_printed = 1;
			last;
		}
	}
	if(not $is_printed){
		for(@genes_ortholog){
			if(exists $hash_genes{$_}){
				print $line, "\n";
				last;
			}
		}
	}

}

close(IN);






