# bioinfScripts
some useful Python scripts for various NGS analyses

These are scripts I regularly use when processing Stacks (or iPyrad) RADseq output data. Alongside the scripts there are some data files that can be processed to see how things work. Instructions for each script follows:

[1] vcf2csv.py – converts Stacks and iPyrad VCF files to CSV files (e.g. to conduct PCA analyses). 
Usage: `python vcf2csv.py batch_0.vcf ouFile.csv`

[2] popWriter.py will create a population map (format: taxon <tab> location) file from a VCF and a master file of taxa and locations (useful for many downstream programs).
Usage: `python popWriter.py batch_0.vcf popinfo.csv outFile.csv`

[3] stacksFastaWriter.py takes Stacks formatted .fa files and writes each locus as an alignment to a separate file in a new folder.
Usage: `python stacksFastaWriter.py batch_0.vcf batch_0.fa fastaFolder`

[4] Concatenator.py will concatenate fasta files that contain the same ordered taxa (i.e. outputs from [3]).
Usage: `python Concatenator.py fastaFolder alignment.fas`

[5] vcf2fineStr.py converts Stacks VCF outputs to fineRADstructure formatted files (NB does not consider phasing of data – see fineRADstructure manual).
Usage: `python vcf2fineStr.py batch_0.vcf popinfo.csv outFile.txt`

[6] ref2fasta.py takes a reference genome (from NCBI but any format should work), an input of chromosome number, and the sequence range (bases) you want chopping out `python ref2fasta.py`
