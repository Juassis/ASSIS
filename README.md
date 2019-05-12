<!-- /TOC -->

-----
#Hardware/Software Requirements
	- 64 bit Linux or Mac OS X
	- x86-64 compatible processors

#Tutorial
 	https://github.com/Juassis/iedb PDF

	RELEASEnotes contains detailed information about the latest major release

#Directory Contents
<b> main dir <b>

		assis.pl, params_model

<b>bin:  scripts<b>

	 compare_orthogroups_by_synteny.pl, export_rnaseq.pl, fix_gff_from_bedtools.sh, gene_overlap_report.pl, gff2fasta.pl, overlap_gff.pl, protein_check.pl, protein_coding_report.pl, rnaseq_report.pl

<b>tests: DataSet Example<b>

	D_simulans_Assessment

		FASTA: Drosophila_simulans.fasta, Drosophila_melanogaster.fasta
		GFF:   Drosophila_simulans.gff3, Drosophila_melanogaster.gff3

		OUTPUT:
			PROTEINS: Drosophila_simulans_proteins.fasta
			REPORT: report_overlap, report_protein_coding

<b>doc: documentation<b>

	Tutorial:PDF
	Paper: Link


#Installation

	 Interpreters: Perl, Python

 - [Perl](https://www.perl.org/get.html)
    * Libraries
         * [BioPerl](https://metacpan.org/pod/BioPerl) (current version)


* [Python 2.7](https://www.python.org/downloads/)
	 * Libraries

		* [Numpy](http://www.numpy.org) (python 2.7)


<span style="padding-left:20px"></span>

	Install dependencies: OrthoFinder, BedTools


####[OrthoFinder:](https://github.com/davidemms/OrthoFinder)
-  <u>[Tutorial](https://github.com/davidemms/OrthoFinder)<u>


<b>Installation <b>
		-  [Latest Development Version](https://github.com/davidemms/OrthoFinder/releases/download/v2.2.7/OrthoFinder-2.2.7.tar.gz)

* MAC:
		 -  [Install the source package](https://github.com/davidemms/OrthoFinder/releases/download/v2.2.7/OrthoFinder-2.2.7_source.tar.gz)  


> All SO
>- ######  Install dependencies: [MCL](http://micans.org/mcl/), [FASTME](http://www.atgc-montpellier.fr/fastme/binaries.php) and [DIAMOND](https://github.com/bbuchfink/diamond/releases) or [MMseqs2](https://github.com/soedinglab/MMseqs2/releases/tag/7-4e23d)
>- ###### Install the [source package](https://github.com/davidemms/OrthoFinder/archive/v2.2.7.tar.gz) if some of the dependencies do not work  

<span style="padding-left:20px"></span>

####[BedTools](http://bedtools.readthedocs.io/en/latest/)

<b>Installation<b>
		- [Latest Version](https://bedtools.readthedocs.io/en/latest/)

* MAC:
[brew](https://brew.sh) install bedtools

- Linux:
apt-get install bedtools

<span style="padding-left:20px"></span>


	ASSIS Software Installation


##[ASSIS](https://github.com/Juassis/iedb)

	# Get latest ASSIS releases from github

	wget https://github.com/Juassis.tar.gz
	tar -xzf XXX.tar.gz
	cd ASSIS

	# Alternatively, get ASSIS using git

	git clone https://github.com/Juassis/iedb.git



## Running the program:
To Run ASSIS use the Example DataSet:

	Project/
		bin
		params_model
		assis.pl
		D_simulans_Assessment
		README.md
		Tutorial.md

		1- cd D_simulans_Assessment		

			FASTA/ Drosophila_simulans.fasta, Drosophila_melanogaster.fasta
			GFF/  Drosophila_simulans.gff3, Drosophila_melanogaster.gff3
			OUTPUT/

		2- perl assis.pl -d D_simulans_Assessment -p params_model



<b>Running with the own DataSet:<b>

Command Line and Options:

	1- mkdir Project

	2- cd Project

		3- mkdir FASTA
			3.1 ln -s genome.fasta (Evaluated and Reference (Closely related species)) FASTA;

		4- mkdir GFF
			4.1 ln -s genome.gff (Evaluated and Reference) GFF;

			* FASTA and GFF should have the same suffix (See Example DataSet for fruit fly genome); *

		 	4.2 This Only for RNASeq Mapped Data (See the tutorial to generate a bam/bed file)


	5- Edit the params file (see the tutorial)

	6- mkdir OUTPUT

	7- perl bin/Assis.pl -d Project -p params_model


---
###Limitations

Splicing leader/ Trans splicing

Fragmented genomes (Finished genomes annotation to help)

This release was tested with the default parameters for drosophila, mouse and parasites genomes. Mammal genomes require at least 4GB of RAM, ideally 16GB.


<!-- Este é um comentário -->
<!--<p>Este <br> é um pará<br>grafo com quebras de linha</p>-->
