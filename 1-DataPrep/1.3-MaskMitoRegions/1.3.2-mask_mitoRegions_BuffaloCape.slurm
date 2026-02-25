#!/bin/bash
#SBATCH --job-name=mask_mitoRegions_BuffaloCape         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2               # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G                # Memory required per CPU
#SBATCH --time=03:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/userid/palaeobovids/refgenomes/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/userid/palaeobovids/refgenomes/%x_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load module and define environment variables
module load bedtools/2.31.1
module load samtools/1.15
module load bwa/0.7.17
BED="/projects/lorenzen/people/userid/palaeobovids/gargammel/paleomix/BuffaloCape/nc_host/BuffaloCape_nc_mitoRegions.mergedbed"

# Navigate to working directory
## For genome concatenated with human genome
cd /projects/lorenzen/people/userid/palaeobovids/refgenomes/competitive_mapping

# Mask mito regions in nuclear reference already concatenated with human reference
echo "$(date) Masking mitochondrial regions in nuclear reference..."
bedtools maskfasta -fi BuffaloCape_Human_nc.fasta -bed $BED -fo BuffaloCape_Human_nc_mitoMask.fasta
# Index new reference
echo "$(date) Making samtools index file..."
samtools faidx BuffaloCape_Human_nc_mitoMask.fasta
# Index new reference with bwa
echo "$(date) Making bwa index files..."
bwa index BuffaloCape_Human_nc_mitoMask.fasta
echo "$(date) Done."

# Remove old reference genome files
rm BuffaloCape_Human_nc.fasta.*
# Compress old reference genome
pigz -p 10 BuffaloCape_Human_nc.fasta
