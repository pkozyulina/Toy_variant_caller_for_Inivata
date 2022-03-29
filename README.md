# Toy Variant Caller for Inivata
A coding task for Inivata from 20220228

What does it do:

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
        ◦ gene name
        ◦ nature of the mutation (e.g. C->T etc)
        ◦ frequency of the mutation
        ◦ number of amplicons supporting the mutation
        ◦ reference sequence for the amplicon
        ◦ mutated sequence for the amplicon
