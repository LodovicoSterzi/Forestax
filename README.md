# Forestax: define bacterial taxonomy using the RecA protein marker

## What does Forestax do?
Forestax is a fast, easy to use tool to assign bacterial taxonomy starting from raw reads, FASTA nucleotide file (such as genome assemblies) or FASTA protein sequences files. Forestax takes the input, extracts RecA protein sequences and assigns the bacterial taxonomy to each sequence using a Random forest algorithm. If raw reads are given, the tool also provides bacterial abundances using the reads coverage. More information about the tool is available (for the moment) at:

  RecA is a reliable marker for bacterial taxonomy, even in the Candidate Phyla Radiation
Lodovico Sterzi, Simona Panelli, Clara Bonaiti, Stella Papaleo, Giorgia Bettoni, Enza D’Auria, Gianvincenzo Zuccotti, Francesco Comandatore
bioRxiv 2024.06.21.600076; doi: https://doi.org/10.1101/2024.06.21.600076

![image](https://github.com/user-attachments/assets/8a84e3e4-b3b7-4cc7-a15d-c2b578a72a65)

## How to install Forestax
Forestax is a Linux-based software. The tool is available for installation by conda or mamba (faster!).

**Conda**:
```
conda create -n forestax  #create virtual environment
conda activate forestax  #activate the virtual environment
conda install --override-channels -c lodovicosterzi -c r -c conda-forge -c bioconda forestax #install the software from the repository
```
**Mamba**:
```
mamba create -n forestax  #create virtual environment
mamba activate forestax  #activate the virtual environment
mamba install -c r -c lodovicosterzi -c conda-forge -c bioconda forestax  #install the software from the repository
```

## How to run Forestax

To see the usage parameters, run:
```
forestax -h
```

### Parameters:
```
Mandatory parameters:
  -o, --output                         Output folder
  -itype, --input_type                 Input type (reads, nucleotide or proteins)

Basic options:
  -r1, --reads_forward                 If the input type are reads, path of FASTQ forward reads file
  -r2, --reads_reverse                 If the input type are reads, path of FASTQ reverse reads file
  -f, --fasta                          If the input type are nucleotide or proteins, path of the FASTA file
  -p, --cpu                            Number of threads (default: 1)
  -tax, --tax_to_assign                Maximum taxonomic level to assign (default: genus)
  -u, --unclassified_threshold         The minimum percentag of confidence necessary to consider an assignment reliable (default: 50)
  -q, --quantify                       Use reads coverage to infer relative abundances of the taxa
  -D, --include_comparison_to_reca_db Disable DIAMOND comparison against the curated database (default: True)

Advanced options:
  -E, --tool_to_extract                Extract RecA reads using kmc (faster) or bowtie2 (slower) (default: kmc)
  -N, --kmc_num_kmers                  Minimum number of known RecA k-mers required to retain a read (default: 30)
  -b, --bowtie_mode                    Bowtie2 mode used to align reads to RecA bowtie database (default: very-sensitive)
  -A, --spades_size_kmers              The size of kmers for SPAdes assembly (default: 55)  
  -r FILE, --reads_paths FILE          Tab-delimited file without a header, containing three columns: input filename, forward reads path, and 
                                       reverse reads path. Required for quantification from nucleotide input. (default: not_given)
  -r, --reads_paths                    If the input is a nucleotide FASTA file (eg. a genome assembly) but wish to quantify the bacterial taxa, 
                                       the tool must be given a text file with THREE TAB-DELIMITED columns and NO HEADER: 
                                       name-of-input-file | path-of-forward-reads-file | path-of-reverse-reads-file
  -k, --keep_track                     Keep track of RAM usage/time for each step. Beware: it will slow down a bit the analysis!
```
### Basic usage
The tool must be given the type of input (reads,nucleotide or protein), the input file/s and the name of the output folder. 

#### From paired FASTQ reads (default)
```
forestax -r1 path/to/reads/*_1.fastq.gz -r2 path/to/reads/*_2.fastq.gz -o example_reads
```
#### From FASTA nucleotide
```
forestax -itype nucleotide -f path/to/nucleotide/*.fasta -o output_nucleotide
```
#### From FASTA proteins
```
forestax -itype proteins -f path/to/proteins/*.faa -o example_proteins
```

### More complex usage
It is also possible to customise some important options. 
For example, the user can ask the tool to:

1) Obtain the coverage of each annotated RecA sequence
```
forestax -r1 path/to/reads/*_1.fastq.gz -r2 path/to/reads/*_2.fastq.gz -o example_reads -q
```
3) Customise the taxonomic depth to reach in the assignment (default is genus).
```
forestax -r1 path/to/reads/*_1.fastq.gz -r2 path/to/reads/*_2.fastq.gz -o example_reads -tax family
```
3) Avoid the reannotation step of unclassified sequences using DIAMOND against a curated reference database
```
forestax -r1 path/to/reads/*_1.fastq.gz -r2 path/to/reads/*_2.fastq.gz -o example_reads -D
```
4) If the user wants to obtain the coverage of each annotated RecA sequence BUT the input is a nucleotide FASTA file (e.g. metagenome assembly), it is possible to do so by adding a tab-delimited plain text file with the name of the input and the paths to the corresponding reads files. 
```
forestax -itype nucleotide -f path/to/nucleotide/*.fasta -o example_nucleotide -q -r reads_path_file.tab
```
With reads_path_file.tab written like this:
```
example1.fasta  path/to/reads/example1_1.fastq.gz  path/to/reads/example1_2.fastq.gz
example2.fasta  path/to/reads/example2_1.fastq.gz  path/to/reads/example2_2.fastq.gz
```
## The output
Forestax will start by creating a folder for the run. In the main folder, there will be: i) an **Output** folder containing the three output files; ii) an **Analysis** folder containing files for each step of the pipeline;  iii) a **log file** containing information about each step of the pipeline.

The **Output** folder contains several files:

The main outputs are:
- *\*.RA.tab*: these files show the relative abundance of taxa at each taxonomic level in the samples (if quantification was asked by the user)     
- *RecA_Taxonomy_OUTPUT_with_unclassified_and_DIAMOND_reannotation.tab*: this file contains the taxonomic assignment and probability at each taxonomic level for each RecA sequence in each input sample. If (-q) has been flagged, it also gives the reads mean coverage for each sequence. Sequences with low taxonomic assignment confidence are marked as UNCLASSIFIED and sequences have been reannotated (when possible) via a DIAMOND sequence comparison step.
- *RecA_Taxonomy_OUTPUT_with_unclassified.tab*: compared to the  previous file, in this file the sequences with low taxonomic assignment confidence ( < than -t parameter) are marked as UNCLASSIFIED, but have NOT been reannotated via DIAMOND sequence comparison. If (-q) has been flagged, it also gives the reads mean coverage for each sequence.
- *RecA_Taxonomy_OUTPUT_raw.tab*: this is the basic output file after the Random Forest taxonomic assignment, without marking low confidence assignments as UNCLASSIFIED and without DIAMOND sequence comparison. If (-q) has been flagged, it also gives the reads mean coverage for each sequence.
- *RecA_Taxonomy_OUTPUT.extended.tab*: this is the basic output file with extra information such as the amino acid sequence, the dayhoff-compressed sequence and the DIAMOND hit coverage values. If (-q) has been flagged, it also gives the reads mean coverage for each sequence.
- *RecA_sequences.fasta*: this file with all the RecA protein sequences (nucleotide sequences are also available in the Analysis/5_RecA_sequences_and_coverage/ folder). This file can be useful for any downstream analysis.

