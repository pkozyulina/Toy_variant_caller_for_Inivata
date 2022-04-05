# Toy SNV Caller as a coding test
A coding task from 20220228.
There were some options on how to approach this task. I could use some more of biopython (classes etc) and make it shorter, however the task appeared to be interesting to me so I decided to do most of the coding from scratch. 

## How to run:
If you are using default annotation, such as:
   * gene1: (55019278, 55205617)
   * gene2: (140734770, 140924566)

Then you can proceed right away to running 
> ./ToySNVcaller.py -i [INPUT1 INPUT2 INPUT3 ...] -r [REFERENCE] -b [BED] -bp [OPTIONAL to generate bar plot]

If you want to create a different annotation, please run 
>./parse_annotation.py -i [INPUT gene1.bed genen2.bed gene3.bed ...] 

where gene1.bed - a file with bed coordinates for gene1. This script should be run prior to ToySNVcaller.py. And after that you can run ./ToySNVcaller.py with -a flag and a new annotation will be provided.


## What does it do:

1. Prepares reference
    * reads *.fa.gz and *.bed files
    * substracts amplicones from *.fasta reference by coordinates in *.bed

2. Evaluates amplicone sequencing data
    * reads *.fa.gz
    * matchs reads from amplicon data to reference
    * reports SNV

3. Deals with SNV
    * calculates VAF (variant allele frequency)
    * vusualizes SNV frequency as bar plot

4. Generates output table containing:
   * gene name
   * nature of the mutation (e.g. C->T etc)
   * frequency of the mutation
   * number of amplicons supporting the mutation
   * reference sequence for the amplicon
   * mutated sequence for the amplicon
