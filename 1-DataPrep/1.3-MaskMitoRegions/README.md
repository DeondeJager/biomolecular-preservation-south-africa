# Masking mitochondrial-like regions in host/target genomes

## Mapping simulated ancient mitochondrial reads
 - Reads were simulated using gargammel (see folder 1.1-Gargammel).
 - Scripts `1.3.1-*` map the simulated reads using paleomix and the `*.yaml` files in this folder.
 - BED files were then generated from the locations where the mitochondrial reads mapped to the nuclear genomes using samtools (the simulated reads were pre-fixed with "MT", making them easy to identify in the BAM file. Exact commands to do this were not saved.
 - Adjacent regions in the BED files were then merged to reduce the number of lines, using `bedtools merge` to generate the `.mergedbed` files used for masking.

## Masking identified mitochondrial-like regions
 - Scripts `1.3.2.-*` in combination with the `*.mergedbed` BED files were used to mask mito-like regions in the fasta reference genomes (i.e. convert these regions to `Ns`) in the concatenated human-target fasta files.