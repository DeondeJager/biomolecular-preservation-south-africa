# Generate stats of genome ito contig lengths, size of genome, etc.
/projects/lorenzen/apps/conda/bbmap-38.84/bin/stats.sh in=GCA_006416875.1_CME_genomic.fasta > GCA_006416875.1_CME_genomic.bbstats.txt

# Sort contigs from longest to shortest
/projects/lorenzen/apps/conda/bbmap-38.84/bin/sortbyname.sh in=GCA_006416875.1_CME_genomic.fasta out=./GCA_006416875.1_CME_genomic.sorted.fasta length descending

# From the sorted fasta, extract the scaffold names of the longest 1,526,878 scaffolds that represent 90% of the length of the genome (N90 as calculated using stats.sh from BBmap)
grep ">SJYK" GCA_006416875.1_CME_genomic.sorted.fasta --max-count=1526878 > GCA_006416875.1_CME_genomic.N90.scaffoldnames.txt

# Count the lines in the GCA_006416875.1_CME_genomic.N90.scaffoldnames.txt

wc -l GCA_006416875.1_CME_genomic.N90.scaffoldnames.txt
# Output: 1526878 GCA_006416875.1_CME_genomic.N90.scaffoldnames.fasta

# Remove everything after the first whitespace
cut -d" " -f1 GCA_006416875.1_CME_genomic.N90.scaffoldnames.txt > GCA_006416875.1_CME_genomic.N90.scaffoldnames.clean.txt

# Then remove the ">" symbol that is in front of each scaffold name:
sed -i 's/>//' GCA_006416875.1_CME_genomic.N90.scaffoldnames.clean.txt

# Extract top 1526878 scaffolds into new fasta file using samtools faidx in combination with xargs:
xargs samtools faidx GCA_006416875.1_CME_genomic.sorted.fasta < GCA_006416875.1_CME_genomic.N90.scaffoldnames.clean.txt > GCA_006416875.1_CME_genomic.N90.fasta

# Count the number of scaffolds in the output file
grep -c ">" GCA_006416875.1_CME_genomic.N90.fasta
# Output: 1526878, which is the expected number

# Then count the number of bases in the file:
grep -v ">" GCA_006416875.1_CME_genomic.N90.fasta | wc | awk '{print $3-$1}'
# Output: 2147716187

# Count the number of bases in the original genome
grep -v ">" GCA_006416875.1_CME_genomic.fasta | wc | awk '{print $3-$1}'
# Output: 2837493372
# And calculate whether the maths adds up to 90% of the original genome
# 2147716187/2837493372 = 0.7569. So it actually represents 75.7% of the original genome, but that's fine for now

# Concatenate these scaffolds into Super_Scaffolds using ScaffoldStitcher and then use that genome for mapping and endogenous content calculations
python /projects/lorenzen/people/pzx702/software/ScaffoldStitcher-1.0/ScaffoldStitcher.py -fasta GCA_006416875.1_CME_genomic.N90.fasta -identifier SJYK -short 36 -nlength 1000 -maxlength 150000000 > GCA_006416875.1_CME_genomic.N90.SS.fasta
# This gives a stitched genome with 25 Super_Scaffolds, which is close to the number of chromosomes eland have (31/32 for male/female).

