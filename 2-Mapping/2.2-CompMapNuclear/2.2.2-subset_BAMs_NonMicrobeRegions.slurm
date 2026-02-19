#!/bin/bash
#SBATCH --job-name=subset_BAMs_nonMicrobe         # Job name
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=2               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/pzx702/palaeobovids/paleomix/libCompare/%x_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/pzx702/palaeobovids/paleomix/libCompare/%x_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules and set up enviroment variables
module load bedtools/2.31.0
module load samtools/1.15

# Subset bam files to remove reads mapped to microbial-like regions
## As identified by Oskolkov et al. (2025) https://doi.org/10.1101/2025.03.19.644176.

## 1. Subset bams for only the target species portion (i.e. exclude human chrs) - also to be used for damage evaluation after removal of microbial-like reads 
cd /projects/lorenzen/people/pzx702/palaeobovids/paleomix/compMap_nc/
HUMANBED="/projects/lorenzen/people/pzx702/palaeobovids/refgenomes/competitive_mapping/Human_Nc.bed"
# TEST: bedtools intersect -a ScA004.BuffaloCape_Human_nc.bam -b $HUMANBED -v > ScA004.BuffaloCape_Human_nc.Buffalo.bam
for file in *.bam; do bedtools intersect -a $file -b $HUMANBED -v > ${file/%.bam/.Target.bam}; done

### Remove duplicate reads and unmapped reads (only marked duplicates with paleomix, didn't remove, so doing it now)
# TEST: samtools view -F DUP,UNMAP -b -o ScA004.BuffaloCape_Human_nc.Buffalo.rmdup.bam ScA004.BuffaloCape_Human_nc.Buffalo.bam
# TEST: samtools view ScA004.BuffaloCape_Human_nc.Buffalo.rmdup.bam | wc -l # Check read count (1723360) matches that of unique reads from paleomix (1723360)...it does :)
for file in *.Target.bam; do samtools view -F DUP,UNMAP -b -o ${file/%.bam/.rmdup.bam} $file; done

## 2. Subset to non-microbial-like regions using bedtools intersect: -v for inverse of microbial regions
# TEST: bedtools intersect -a ScA004.BuffaloCape_Human_nc.Buffalo.rmdup.bam -b $REGIONS/microbeLike_Buffalo.bed -v > ScA004.BuffaloCape_Human_nc.Buffalo.rmdup.noMicro.bam
# TEST: samtools view ScA004.BuffaloCape_Human_nc.Buffalo.rmdup.noMicro.bam | wc -l > ScA004.BuffaloCape_Human_nc.Buffalo.rmdup.noMicro.readCount.txt # Reads = 1716997, thus 1723360-1716997 = 6363 reads removed that overlapped with microbial-like regions.
REGIONS=/projects/lorenzen/people/pzx702/palaeobovids/refgenomes
### Buffalo samples
for file in ScA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Buffalo.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in ScA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads
for file in PaA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Buffalo.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in PaA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads
### Reedbuck & rhebok samples
for file in RxA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Reedbuck.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in RxA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads
for file in PcA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Reedbuck.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in PcA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads
### Eland samples
for file in ToA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Eland_N90.SS.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in ToA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads
### Grysbok samples
for file in RmA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Grysbok.sheep.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in RmA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads

### Combine all readCount files into one, with count, sample name, file name in tab-separated text file
for f in *readCount.txt; do printf "%s\t%s\t%s\n" "$(cat "$f")" "${f:0:6}" "$f"; done > all.noMicro.readcount.txt
### Delete intermediate bam files
rm *.Target.bam
rm *.rmdup.bam

## 3. For libCompare samples (SCR vs BEST libraries analysis)
cd /projects/lorenzen/people/pzx702/palaeobovids/paleomix/libCompare

### 3.1. Subset bams for target only
HUMANBED="/projects/lorenzen/people/pzx702/palaeobovids/refgenomes/competitive_mapping/Human_Nc.bed"
for file in *.bam; do bedtools intersect -a $file -b $HUMANBED -v > ${file/%.bam/.Target.bam}; done

### 3.2. Remove duplicate and unmapped reads
for file in *.Target.bam; do samtools view -F DUP,UNMAP -b -o ${file/%.bam/.rmdup.bam} $file; done

### 3.3. Subset to non-microbial-like regions
REGIONS=/projects/lorenzen/people/pzx702/palaeobovids/refgenomes
#### Reedbuck libs
for file in RxA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Reedbuck.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in RxA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads
#### Buffalo libs
for file in ScA*.Target.rmdup.bam; do bedtools intersect -a $file -b $REGIONS/microbeLike_Buffalo.bed -v > ${file/%.bam/.noMicro.bam}; done
for file in ScA*.noMicro.bam; do samtools view $file | wc -l > ${file/%.bam/.readCount.txt}; done # count reads

### 3.4. Combine all readCount files into one, with count, sample name, file name in tab-separated text file
for f in *readCount.txt; do printf "%s\t%s\t%s\n" "$(cat "$f")" "${f:0:15}" "$f"; done > all.noMicro.readcount.txt
### Delete intermediate bam files
rm *.Target.bam
rm *.rmdup.bam

