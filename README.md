===============================
MetaShot (Metagemomics Shotgun)
===============================

MetaShot (Metagenomics Shotgun) is a pipeline designed for the complete taxonomic assessment of the human microbiota.
In MetaShot, third party tools and new developed Python and Bash scripts are integrated to analyze paired-end (PE) Illumina reads, offering an automated procedure covering all the analysis steps from raw data management to taxonomic profiling. MetaShot is designed to analyze both DNA-Seq and RNA-Seq data. 

Pipeline
========
The MetaShot analysis procedure can be divided in three main processes:
1.	Pre-processing procedures: input sequences containing low-quality/complexity regions and reads shorter than 50 nucleotides are removed. This is performed by applying FaQCs [1]. Moreover, is also possible to remove the phage PhiX sequences. This procedure is optional because it is time and computational extensive (FaQC with PhiX removal is 7-time slower). Only reads passing all the quality check filters are directed to the following steps.
2.	Comparison with the human genome: cleaned data are mapped against the human genome (hg19, GRCh37, 2009) [2] and transcriptome (UCSC RefSeq) by applying STAR [3].
3.	Comparison with reference databases and taxonomic annotation: in MetaShot are implemented four reference collections for Prokaryotes, Viruses, Fungi and Protista, obtained by processing data from GenBank and RefSeq databases (for more information see the “Division data creation” section). The taxonomic annotation is a three step procedure: 
a.	Identification of candidate microbial sequences: in order to reduce the processing time and the computational requirements, and to improve the assignment accuracy, the PE reads passing the denoising procedure are compared to the four reference collections, by applying Bowtie2 [4] with the default options (a compromise between speed and mapping accuracy). The mapping data are filtered and PE reads obtaining at least one match with a hamming distance [5] below 10% of the sequence length are retained. In this way MetaShot identifies a subset of reads to be addressed to the microbial taxonomic classification.
b.	Comparison with reference databases and taxonomic annotation: microbial candidate sequences are compared to the four reference collections by using Bowtie2. In this step, the mapping accuracy is increased by modifying the length of seed substrings (-L option set to 20, the default value is 22) and the number of alignments per read (-k option set to 100, default 1). The resulting alignments are filtered according to identity percentage (threshold ≥ 97%) and query coverage (threshold ≥ 70%). The obtained results for the four different taxonomic divisions and for human genome and transcriptome are then intersected to remove the sequences that map ambiguously across them. Finally, the sequences are taxonomically annotated by using the TANGO (Taxonomic assignment in metagenomics) [6] tool.
c.	Human Endogenous Retrovirus (HERV) identification: in order to identify HERV sequences, the PE reads labeled as ambiguous are parsed to identify those mapping only to human genome/transcriptome and to the viral division. These PE reads are subsequently taxonomically annotated by applying TANGO on the mapping information obtained against the viral division. If the obtained classification is under the HERV group it is accepted, otherwise the PE will be considered ambiguous.
4.	Report generation: a CSV file, an HTML interactive table summarizing the taxonomic assignment and a Krona graph [7] of the obtained taxonomy are provided for each division.

Division data creation
======================
The reference collections for Prokaryotes, Viruses, Fungi and Protista have been built by following a common procedure, with some specific add-ons, implemented in a Bash and Python pipeline:
-	For each collection the GenBank and RefSeq flat-file and FASTA files were downloaded from the NCBI ftp site (ftp://ftp.ncbi.nlm.nih.gov). For Prokaryotes and Viruses two specific GenBank divisions (BCT and VRL) are available. For Fungi, the “Plantae” GenBank division (PLN) was downloaded and the fungal entries were identified by parsing the taxonomic information of each entry. For Protista, the “Invertebrate” GenBank division (INV) was downloaded and the Protista entries were collected by suitably parsing the taxonomic information of each entry. Regarding the RefSeq data, specific divisions are available for Prokaryotes, Viruses and Fungi. As for GenBank data, the Prostista collection was created by downloading and parsing the “invertebrate” RefSeq division.
-	The entire NCBI taxonomy dump data were downloaded from the NCBI ftp site  (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz). 
-	The Genbank flat-file was parsed in order to associate the accession number of each entry to a taxonomy identifier (taxid) in the NCBI taxonomy, by extracting the information annotated in the feature table “source” field. An entry was discarded if resulting associated to more than one microbial division or the associated taxid was not present in the reference taxonomy due to inconsistency between Genbank and Taxonomy databases. For the entries passing the taxonomy-based filter, the sequence data were annotated and the accession number – taxid association was stored in a dump.
-	FASTA sequences were indexed by applying the bowtie2-build program.
-	By using the dump containing the accession number – taxid association the TANGO reference taxonomy was built.
The current MetaShot reference collections were built starting from the releases 209 and 72 of the GenBank and RefSeq databases, respectively.

References
==========
1.	Lo, C.C. and P.S. Chain, Rapid evaluation and quality control of next generation sequencing data with FaQCs. BMC Bioinformatics, 2014. 15: p. 366.
2.	Rosenbloom, K.R., et al. The UCSC Genome Browser database: 2015 update. Nucleic Acids Res, 2015. 43(Database issue): p. D670-D681.
3.	Dobin, A., et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 2013. 29(1): p. 15-21.
4.	Sayers, E.W., et al., Database resources of the National Center for Biotechnology Information. Nucleic Acids Res, 2009. 37(Database issue): p. D5-15.
5.	He, M.X., S.V. Petoukhov, and P.E. Ricci, Genetic code, hamming distance and stochastic matrices. Bull Math Biol, 2004. 66(5): p. 1405-21.
6.	Alonso-Alemany, D., et al., Further Steps in TANGO: improved taxonomic assignment in metagenomics. Bioinformatics, 2014. 30(1): p. 17-23.
7.	Ondov, B.D., N.H. Bergman, and A.M. Phillippy, Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 2011. 12: p. 385.



