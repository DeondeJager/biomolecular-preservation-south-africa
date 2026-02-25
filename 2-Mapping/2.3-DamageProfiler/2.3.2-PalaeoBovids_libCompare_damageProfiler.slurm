#!/bin/bash
#SBATCH --job-name=damageProfiler_libCompare         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/userid/palaeobovids/paleomix/libCompare/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/userid/palaeobovids/paleomix/libCompare/%x_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules and set up enviroment variables
module purge
module load samtools/1.15

# Navigate to working directory and set up variables
cd /projects/lorenzen/people/userid/palaeobovids/paleomix/libCompare
REF=/projects/lorenzen/people/userid/palaeobovids/refgenomes/competitive_mapping

## Cape buffalo samples
for file in ScA*.noMicro.bam; do
    base=${file%.bam}
    out=${base/.BuffaloCape_Human_nc.Target.rmdup.noMicro/}
    java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/BuffaloCape_Human_nc_mitoMask.fasta
done

## Reedbuck samples
for file in RxA*.noMicro.bam; do
    base=${file%.bam}
    out=${base/.ReedbuckBh_Human_nc.Target.rmdup.noMicro/}
    java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/ReedbuckBh_Human_nc_mitoMask.fasta
done
