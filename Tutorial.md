
- [ASSIS](#titulo)
	- [Evaluation the annotation](#subtitulo)

	<!-- /TOC -->



<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Genome Annotation](#titulo)
- [Tasks](#subtitulo)
- [Motivation](#subtitulo)
- [Purpose](#subtitulo)
- [What does ASSIS do?](#titulo)
- [Input Files](#titulo)
- [Output Files](#titulo)
  - [Recommendation](#subtitulo)
- [The three modules of ASSIS](#titulo)
	- [Syntenic](#subtitulo)
	- [Structural](#subtitulo)
	- [Transcriptional](#subtitulo)



<!-- /TOC -->


# ASSIS
## Evaluation the genome annotation


##Tasks

<p align="justify">Implement and improve an automatic procedure based on orthologs, synteny and transcript evidences to evaluate the whole <b>genome annotation</b><sup>1<p align="justify">


##Motivation

<p align="justify">The main goal of this pipeline is to perform the evaluation of the  <b>whole genome annotation</b> without too much humans efforts. For two reasons: humans evaluation can introduce some bias and the time is very important in a big data era that we are living.


##Purpose

<p align="justify">ASSIS is a genome annotation evaluation pipeline. Its purpose is to allow some eukaryotic genome annotation projects to evaluate their genes quality. ASSIS identifies problematics genes encoding proteins and automatically synthesizes these misannotation genes in a GFF3 format.

<p align="justify">The system allow to perform misannotation detection using orthologs syntenic information from closely related species producing higher quality and confidential score by gene. They can also be perform two more analysis divided in different categories: 1- structural issues reports and 2- completeness RNAseq coverage statements. The pipeline permits the RNAseq module to be worked as additional analysis, assigning an extra score per module. In total they have three modules: <u>Syntenic, Structural and
Transcriptional</u>.



<p align="justify">ASSIS should prove especially useful to evaluate new annotation genomes projects, however it can be applied to any other available genome annotation. No great experience in bioinformatics is required, nevertheless we strongly encourage to study the genome characteristic before start the analysis to increase the quality scores in a simple report by gene.

<p align="justify">This Pipeline focuses on structural annotation from the whole genome annotation. <sup>1

<p align="justify">ASSIS inputs are a genome file and the outputs can be directly loaded into the Apollo or Artemis to perform a manual curation. The misannotation genes receive a low scores, while the trusted genes are assigned with a higher score.


-----
##What does ASSIS do?

	-> Identifies misannotation genes encoding proteins;

	-> Produces evidence-based quality values for each gene in the genome annotation;

	-> Synthesizes these information into final report;

	-> Produces new annotation file with problematic genes to be correct;


<p align="justify">No gene is removed from the annotation, only a report is generated containing all the scores assigned to each gene. The output files generated are in Tabular format of type text and GFF3 format.

-------

## The three modules of the ASSIS:

<p align="justify">The general idea to classify the genes encoding protein by quality score. To achieve this goal, three different methodologies where developed. The main methodology is the syntenic information from orthologs. A high level of synteny is expected in closely related species although we know that it is not true for all those orthologs genes.
<p align="justify"> Theoretically two orthologs genes should share
homologous neighboring genes, however, there is a chance of non homology matches
occurring in collinearity. In order to account the score  for these expected events, two
other modules have been developed in ASSIS pipeline: Structural and
Transcriptional to increase the accuracy of the scores attributed.

The score  goes to 0-5

### <u>1 Syntenic:


><p align="justify">Evaluated genome     -----> Provided by the user to evaluate the quality of the genes.
><p align="justify">Reference genome  -----> Closely related species provided by the user to evaluate the reference genome.


<p align="justify">Closely related species are added for evaluation using the correspondence of orthologs. Prediction of orthologs is performed through the OrthoFinder program,

<!-- Este é um comentário  nevertheless if the user prefers, the orthofinder's output can be provided as a input file for the ASSIS pipeline (see the tutorial).-->

<p align="justify"> The output from <b>OrthoFinder</b>: Orthogroups.csv, is the work input file to ASSIS pipeline. Orthogroups.csv is a tab separated text file. Each row contains the genes belonging to a single orthogroup. The genes from each orthogroup are organized into columns, one per species. Since the orthologos can be one-to-one, one-to-many or many-to-many in ASSIS pipeline these categories are divided in: <b>Unique (one-to-one), Multiple (one-to-many and many-to-many) and Unassigned (genes that were not assigned to any orthogroup)</b>.

<p align="justify"> We basically start from the premise that two orthologs genes should share homologous neighboring genes. The idea is keep the continuous genes inside the block. The comparison is performed gene by gene using the annotation file (GFF3) for each specie. The GFF3 files are sorted by Chromosome/Scaffold, based on the genomic coordinates of the mRNA feature. The first gene ordered into the genome reference is crossed to the set of orthologous in the OrthoFinder's csv file.

<p align="justify">The syntenic block starts when the first match between the gene's reference ordered and the set of Unique orthologs occurring. The next step is to continue the comparison of the orthologous neighborhood, including all orthologs information: Unique, Multiple and Unassigned groups.  The Syntenic block is broken when two continuous mismatch occur.

####The pipeline rules are as follows:


	1. Search for the first single gene in the reference chr/scaffold. Search in the ortholog file, the gene that belongs to the same ortholog group.

	2. Move to the next gene in the corresponding chr/scaffold.

	2.1. If both reference and ortholog genes have the same ortholog group, we score 5 and move on to the next gene.
	2.2. If the reference gene has no orthologous group, we give a score 2 for this gene and move the cursor to the next reference gene.
	2.3. If the orthologous gene has no orthologous group, we give a score 2 for this gene and move the cursor to the next orthologous gene.
	2.4  If both genes have orthologous groups and this ortholog group is not equal:
		2.4.1. We check if the reference ortholog group exists in the window of the next 4 genes of the ortholog, if it does not exist, we score this gene with 2. The same logic is for the orthologous gene.
		2.4.2. If the orthologous reference gene group exists in the ortholog window and the orthologous gene exists in the reference window, we build a block of 5 genes (the current comparison gene, plus the next four) in both orthologs and we score this group as follows:
			2.4.2.1 If all orthogroups in both windows are identical, we score 4 for all genes in this window
			2.4.2.2 If there is only one group of orthologs in the different reference, we scored 3; if they have more than one different, we scored 2 for the whole window


...... (definir mismatch) If they not find the correspondence the new block start to be built. It is expected that the larger the block have a better the annotation quality.

Micro inverssões,desenhar etc

This approach is able to detect underpredicted and overpredicted genes/regions generated by automatic annotation.  

The score is attributed by gene isoform? . The weight assignment of the score is a... the weight is equally distributed among all categories.????

####The pipeline Score are assigning as follows:

	Score 0: Scaffolds with less than 10 genes
	Score 1: Genes from the beginning of the scaffold belonging to the multigenic families
	Score 2: Reference genes that do not belong to the next ortholog window
	Score 2: Genes in the reference comparison window with more than 2 genes outside the ortholog window
	Score 3: Genes in the reference comparison window with only one gene outside the ortholog window
	Score 4: Genes in the reference comparison window where all genes are in the orthologous comparison window, has one or more genes out of order
	Score 5: Gene identical to orthologous compared gene

-------

### <u>2 Structural
<p align="justify">The structural evaluation consists in the analysis of the genes structure.
Remember: Structural genome annotation is the process of identifying genes and their intron–exon structures.

The pipeline is able to detected:


	->  Start and Stop codons missing;

		-> Frameshift;

		-> Gene Length (length of the gene annotated - RBH in the UniProt database)

	-> Nested genes (genes with overlap);


<p align="justify">All likely problematic genes are reported to the researchers for inspect them.


-------

### <u>3 Transcriptional

<p align="justify">Transcript evidence are used to improve the score quality by gene.

<p align="justify">The RNASeq data must be previously mapped by the user against the genome reference. The ASSIS input file needs to be in a BED format with at least "mRNA" and "CDS" features, generated by the Bedtools software.

<p align="justify">Remember: The sequence IDs in BED file must correspond (name and coordinate) to the GFF file. <b>See how to generated a BED file</b>

<p align="justify">Gene coverage is accounted for through the BedCoverage of the Bedtools package.

-----



## Input Files:

Genome Sequence. (a genome is a .fasta file and a .gff file)



>**- GFF3**
**- FASTA**





 ### <span style="color:red;">Expected sequence format:</span>

Sequences must be in FASTA format, with one-word headers without special characters such as ?<";|\{}[]+=- etc. The width of the individual lines does not matter (i.e. whole sequences on one line are fine, as are sequences separated over many lines).

######GFF3 preferably needs to be in an [ENSEMBL format](https://www.ensembl.org/info/website/upload/gff3.html)



Expected GFF3 structure:

The features in the input GFF3 file should have the following minimal type structure (e.g. should be connected like this in via Parent attributes):

		gene				  pseudogene              
		└───mRNA                            └─── pseudogenic_transcript                
		    └───CDS                             └───pseudogenic_exon




The GFF3 file should contain features identified as "gene", "mRNA, "exon", "CDS" "intron", "start_codon", "stop_codon". However it's mandatory only the "gene", "mRNA" and "CDS" feature information.

The file must be strictly valid GFF3! Make sure your file is validated, or when in doubt clean it up with GenomeTools:

	                gt gff3 -sort -retainids -tidy input.gff3 > input.clean.gff3

Then, continue working with the cleaned file input.clean.gff3. The command above will also provide error messages if the file is broken in such a way that it cannot be automatically fixed.


 ### <span style="color:red;">Recommendation:</span>

<p align="justify">We wisely recommend GFF file validation before parsing starts. We suggest using GFFvalidator from genometools. The genometools software can also be used to correct simple errors in gff3 files prior to utilizing ASSIS.

######GFFvalidator
http://genometools.org/cgi-bin/gff3validator.cgi




##Output Files:

>**- GFF3**
>**- CSV File**


---


##Generate BAM/Bedtools file coverage


bedtools coverage -s -a file.gff3 -b file.bam  > newfile.gff3


---

##Glossary

**1 Genome Annotation:**

<p align="justify"> A term used to describe two distinct processes. Structural: genome annotation is the process of identifying genes and their intron–exon structures. Functional: genome annotation is the process of attaching meta-data such as gene ontology terms to structural annotations (Yandel, 2012).

-----
##License

<!--https://github.com/sanger-pathogens/companion/wiki/Preparing-reference-data-sets-->



<!-- Este é um comentário -->


<!-- mudando a cor -->

<!-- Roses are <span style="color:red; font-family:Georgia; font-size:2em;">red.</span> -->





-----------
<!-- **Notas Avulsas**

 Este é um comentário
its permits the manual curation to be faster

Advantage above GeneValidator (whole genome)

bedtools coverage -s -a file.bam -b file.gff > file.bedtools
-g???? -d -hist
inverter a com b



overprediction of genes
orthologous clusters
 each with a single copy in each of the organisms

First Header | Second Header
------------ | -------------
Content from cell 1 | Content from cell 2
Content in the first column | Content in the second column
-->
