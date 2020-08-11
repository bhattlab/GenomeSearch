# *GenomeSearch*
A command line tool to quickly identify closely related genomes using a marker-gene based approach.

# About
~500,000 RefSeq/GenBank (RefBank) and ~280,000 Unified Human Gastrointestinal Genome (UHGG) were downloaded and analyzed to identify [phylophlan](https://huttenhower.sph.harvard.edu/phylophlan/) marker genes. These marker genes are stored in a database and all mapped back to the original genome. Given a query genome, *GenomeSearch* extracts marker genes from the query genome and searches against the selected database to quickly identify closely related genomes. This can also be run in `meta` mode where each contig in the FASTA file is treated independently.

# Installation

Just type in 
    
    pip install genomesearch
    
And then

    genomesearch install
    
It will then prompt you to enter a number between 1 and 400, representing the number of marker gene databases you want to have available. If you plan on using *GenomeSearch* in *meta* mode, I recommend 400. If you plan on using it only to analyze isolates, the default of 150 should suffice.

That's it, you are ready to go!

# Usage

If you want to search a single isolate genome against the RefBank database, type:

    genomesearch isolate refbank <myIsolate.fna>
    
If you want to search against the UHGG database, which is a detailed database of Human Microbiome species, use:

    genomesearch isolate uhgg <myIsolate.fna>
    
If you want to perform a metagenomic search against refbank, use:

    genomesearch meta refbank <myMetagenome.fna>
    
And if you want to perform a metagenomic search against uhgg, use:

    genomesearch meta uhgg <myIsolate.fna>

# Citation
