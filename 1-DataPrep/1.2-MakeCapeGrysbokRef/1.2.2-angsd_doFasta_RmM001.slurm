#!/bin/bash
#SBATCH --job-name=doFasta_RmM001     # Job name
#SBATCH --mail-type=BEGIN,END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user= # Where to send mail	
#SBATCH --nodes=1                      # Run all processes on a single node
#SBATCH --ntasks=2                     # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G               # Memory required per CPU
#SBATCH --time=02:00:00                # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/userid/palaeobovids/raphicerus/fasta/angsd_doFasta_RmM001_%j.out    # Standard output log
#SBATCH --error=/projects/lorenzen/people/userid/palaeobovids/raphicerus/fasta/angsd_doFasta_RmM001_%j.err     # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used 
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules
module load angsd/0.921
module load samtools/1.15 # For bgzip
module load bbmap/39.01 # For stats

# Navigate to working directory and set up variables
cd /projects/lorenzen/people/userid/palaeobovids/raphicerus/bams
OUTDIR=/projects/lorenzen/people/userid/palaeobovids/raphicerus/fasta

# Insert commands here:
## Write RmM001 fasta file based on the alignment to the sheep genome
angsd -i RmM001_sheep.MD.MQ30.rmMasked.bam \
 -doFasta 2 \
 -doCounts 1 -setMinDepth 5 -setMaxDepth 67 -minQ 25 \
 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
 -out $OUTDIR/RmM001_sheep_rmMasked_20231221

## Unzip, rename and index
bgzip -d -@ 2 $OUTDIR/RmM001_sheep_rmMasked_20231221.fa.gz
mv $OUTDIR/RmM001_sheep_rmMasked_20231221.fa $OUTDIR/RmM001_sheep_rmMasked_20231221.fasta
samtools faidx $OUTDIR/RmM001_sheep_rmMasked_20231221.fasta
bbstats.sh in=$OUTDIR/RmM001_sheep_rmMasked_20231221.fasta out=$OUTDIR/RmM001_sheep_rmMasked_20231221.bbstats.txt
