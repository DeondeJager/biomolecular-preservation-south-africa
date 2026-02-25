#!/bin/bash
#SBATCH --job-name=MCWorkflow_ElandCommon         # Job name
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=24               # Number of CPU cores per task
#SBATCH --mem-per-cpu=20G                # Memory required per CPU
#SBATCH --time=72:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/userid/palaeobovids/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/userid/palaeobovids/%x_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Activate conda environment (if module is not available)
source /home/userid/.bashrc
conda activate MCWorkflow

# Load modules and set up enviroment variables
module load samtools/1.15
module load gcc/11.2.0 # For R
module load R/4.3.3
module load bowtie2/2.5.4

# Navigate to working directory
cd /projects/lorenzen/people/userid/palaeobovids/MCWorkflow

# Insert commands here:
./micr_cont_detect.sh GCA_006416875.1_CME_genomic.N90.SS_mitoMask.fna.gz /projects/lorenzen/people/userid/palaeobovids/MCWorkflow/data GTDB 24 GTDB_sliced_seqs_sliding_window.fna.gz 10

# Deactivate conda environment (if conda is used instead of modules)
conda deactivate

