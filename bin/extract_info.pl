
use Data::Dumper;

my $file_gff = shift or die "perl $0 file_gff file_orthofinder species_index";
my $file_orthofinder = shift;
my $species_index = shift;
my $has_header = shift;

my %hash_orthofinder = ();
my %hash_uniq = ();


my %hash_scaffolds = ();
my @array_scaffolds_order = ();


read_orthofinder();
read_gff();
sort_gff();
print_output();


sub print_output{
	for my $scaffold (@array_scaffolds_order){
		my $array_genes = $hash_scaffolds{$scaffold};
		for (@$array_genes){
			my $gene = $_->{ID};
			my $parent = $_->{PARENT};
			if(exists $hash_uniq{$gene}){
				print("$scaffold\t$gene\t",$hash_uniq{$gene},"\tUniq\t$parent\n");
			}else{
				if(exists $hash_orthofinder{$species_index}{$gene}){
					print("$scaffold\t$gene\t",$hash_orthofinder{$species_index}{$gene},"\tNotUniq\t$parent\n");
				}else{
					print("$scaffold\t$gene\tNo OG\tNotUniq\t$parent\n");
				}
			}
		}
	}
}

sub sort_gff{
	for my $scaffold (@array_scaffolds_order){
		my $array_genes = $hash_scaffolds{$scaffold};
		my @array_ordered = sort {$a->{GENE_START} <=> $b->{GENE_START} || $a->{GENE_END} <=> $b->{GENE_END} }  @$array_genes;
		$hash_scaffolds{$scaffold} = \@array_ordered;
	}
}


sub read_gff{
	open(IN, $file_gff);
	while(<IN>){
		chomp;
		unless(/^#/){
			my @cols = split /\t/;
			if($cols[2] eq "mRNA"){
				my $scaffold = $cols[0];
				my $start = $cols [3];
				my $end = $cols [4];
				my $strand = $cols[6];
				my $attributes = $cols[8];
				my $gene;
				my $parent;
				
				if($attributes =~ m/^ID=(transcript:)?(.+);Parent=(gene:)?([^;]+).*$/){
					$gene = $2;
					$parent = $4;

					my %hash_gene = ();
                    $hash_gene{ID} = $gene;
                    $hash_gene{STRAND} = $strand;
                    $hash_gene{PARENT} = $parent;
                    if($start > $end){
                        $hash_gene{GENE_START} = $end;
                        $hash_gene{GENE_END} = $start;
                    }else{
                        $hash_gene{GENE_START} = $start;
                        $hash_gene{GENE_END} = $end;
                    }

					if(exists $hash_scaffolds{$scaffold}){
						my $array_genes = $hash_scaffolds{$scaffold};
						push(@$array_genes, \%hash_gene);
					}else{
						my @array_genes = (\%hash_gene);
						$hash_scaffolds{$scaffold} = \@array_genes;
						push(@array_scaffolds_order, $scaffold);

					}
				}
				
			}
		}
	}
	close(IN);
}

sub read_orthofinder {
	open(IN, $file_orthofinder);

	if($has_header){
		<IN>;
	}
	while(<IN>){
		chomp;
		my @cols = split /\t/;
		my $og = $cols[0];
		my $is_uniq = 1;
		if($cols[1] ne ""){
			my @array_genes = split(/, /, $cols[1]);
			if(scalar @array_genes > 1){
				$is_uniq = 0;
			}
			for(@array_genes){
				$hash_orthofinder{1}{$_} = $og;
			}
		}
		if($cols[2] ne ""){
			my @array_genes = split(/, /, $cols[2]);
			if(scalar @array_genes > 1){
				$is_uniq = 0;
			}
			for(@array_genes){
				$hash_orthofinder{2}{$_} = $og;
			}
		}
		if($is_uniq){
			$hash_uniq{$cols[$species_index]}=$og;
		}
	}
	close(IN);
}