#!/bin/bash
#SBATCH --job-name=damageProfiler_compMap_nc         # Job name
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=08:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/pzx702/palaeobovids/paleomix/compMap_nc/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/pzx702/palaeobovids/paleomix/compMap_nc/%x_%j.err         # Standard error log

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
cd /projects/lorenzen/people/pzx702/palaeobovids/paleomix/compMap_nc
REF=/projects/lorenzen/people/pzx702/palaeobovids/refgenomes/competitive_mapping

## Split bam files for samples with multiple libraries, so that DamageProfiler can be done on lobs separately
### These are (list saved in file multiLib_samples.txt and list of bams saved in multiLib_samples_bams.txt):
### PaA004 - no true endo reads
### PaA011 - no true endo reads
### PaA012 - no true endo reads
### PaA013 - no true endo reads
### PcA005 - no true endo reads
### RxA007 - true endo reads
### ScA005 - true endo reads
### ScA009 - no true endo reads
### ToA008 - no true endo reads
while read -r line; do samtools split -v -f '%*_%!.%.' ${line}; done < multiLib_samples_bams.txt

# In slurm job
## Long-horned buffalo
### All samples (incl. the combined bams for those with multiple libraries):
for file in PaA*.noMicro.bam; do java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${file:0:6}.damageProfiler -r $REF/BuffaloCape_Human_nc_mitoMask.fasta; done
### Samples with multiple libraries:
for file in PaA*.noMicro_*.bam; do 
	base=${file%.bam}
	out=${base/.BuffaloCape_Human_nc.Target.rmdup.noMicro/}
	java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/BuffaloCape_Human_nc_mitoMask.fasta
done

## Cape buffalo samples
### All samples (incl. the combined bams for those with multiple libraries):
for file in ScA*.noMicro.bam; do java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${file:0:6}.damageProfiler -r $REF/BuffaloCape_Human_nc_mitoMask.fasta; done
### Samples with multiple libraries:
for file in ScA*.noMicro_*.bam; do
    base=${file%.bam}
    out=${base/.BuffaloCape_Human_nc.Target.rmdup.noMicro/}
    java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/BuffaloCape_Human_nc_mitoMask.fasta
done

## Grey rhebok samples
### All samples (incl. the combined bams for those with multiple libraries):
for file in PcA*.noMicro.bam; do java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${file:0:6}.damageProfiler -r $REF/ReedbuckBh_Human_nc_mitoMask.fasta; done
### Samples with multiple libraries:
for file in PcA*.noMicro_*.bam; do
    base=${file%.bam}
    out=${base/.ReedbuckBh_Human_nc.Target.rmdup.noMicro/}
    java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/ReedbuckBh_Human_nc_mitoMask.fasta
done

## Reedbuck samples
### All samples (incl. the combined bams for those with multiple libraries):
for file in RxA*.noMicro.bam; do java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${file:0:6}.damageProfiler -r $REF/ReedbuckBh_Human_nc_mitoMask.fasta; done
### Samples with multiple libraries:
for file in RxA*.noMicro_*.bam; do
    base=${file%.bam}
    out=${base/.ReedbuckBh_Human_nc.Target.rmdup.noMicro/}
    java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/ReedbuckBh_Human_nc_mitoMask.fasta
done

## Eland samples
### All samples (incl. the combined bams for those with multiple libraries):
for file in ToA*.noMicro.bam; do java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${file:0:6}.damageProfiler -r $REF/ElandCommon_Human_nc_mitoMask.fasta; done
### Samples with multiple libraries:
for file in ToA*.noMicro_*.bam; do
    base=${file%.bam}
    out=${base/.ElandCommon_Human_nc.Target.rmdup.noMicro/}
    java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${out}.damageProfiler -r $REF/ElandCommon_Human_nc_mitoMask.fasta
done

## Cape grysbok samples
### All samples (there are none with multiple libraries):
for file in RmA*.noMicro.bam; do java -Xmx8g -jar ~/programs/damage_profiler/DamageProfiler-1.1-java11.jar -i $file -o ${file:0:6}.damageProfiler -r $REF/GrysbokCape_Human_nc_mitoMask.fasta; done

