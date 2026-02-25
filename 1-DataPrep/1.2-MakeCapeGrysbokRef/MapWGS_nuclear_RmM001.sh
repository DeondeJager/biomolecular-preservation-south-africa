#!/bin/bash
#source /home/userid/.bashrc
# Map modern WGS reads

## Modules are loaded in the slurm script used to execute this script
### Modules required: bwa,samtools,gatk4,java-jdk/8.0.112

## Reads should be cleaned before mapping with fastp

## $1 - Threads
## $2 - Clean data path
## $3 - Clean reads parent name (without _1.clean.fq.gz and _2.clean.fq.gz)
## $4 - Code name for sample
## $5 - Results folder
## $6 - Path to reference (including file name)
## $7 - Reference code name

## Sanity check for arguments
echo "Sanity check for provided arguments:"
echo "Threads: $1"
echo "Clean data path: $2"
echo "Clean reads parent name: $3"
echo "Code name for sample: $4"
echo "Results folder: $5"
echo "Reference genome: $6"
echo "Reference genome code name: $7"

cd $5

#----------STEP 1----------#
## Get read group information from fastq files for GATK
### Need: ID,PU,SM,PL,LB (see: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671)
header=$(zcat $2/$4/$3_1.clean.fq.gz | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g')
barcode=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
echo "Read Group @RG\tID:$id\tPU:$id"."$barcode\tSM:$4\tPL:ILLUMINA\tLB:$barcode"
#
### Index reference genome and map clean reads - produce sorted bam with mate score tags
bwa index $6
samtools faidx $6
date "+%F %T"; echo "Mapping reads to reference genome..."
bwa mem -M -a -t $1 \
 -R $(echo "@RG\tID:$id\tPU:$id"."$barcode\tSM:$4\tPL:ILLUMINA\tLB:$barcode") \
 $6 $2/$4/$3_1.clean.fq.gz $2/$4/$3_2.clean.fq.gz | \
 samtools fixmate -@ $1 -m -u - - | \
 samtools sort -@ $1 -m 7G -O bam -o $5/${3}_${7}.sorted.bam -
date "+%F %T"; echo "Done mapping reads and merging temp bam files."
### Index
date "+%F %T"; echo "Indexing bam file..."
samtools index -@ $1 $5/${3}_${7}.sorted.bam
date "+%F %T"; echo "Done indexing bam file."

#----------STEP 2----------#
# For RmM001, the above steps were done separately for each lane, because the library was sequenced over two lanes.
# Then the bam files were merged before duplicate marking, because it is still the same library.
# So everything below this line was commented-out at first, and then everything above this line was commented out to continue.
date "+%F %T"; echo "Merging bam files from the same library but different lanes..."
samtools merge -@ $1 -o $5/${4}_${7}.sorted.merge.bam $5/RmM001_FKDN230042835-1A_HNCK2DSX3_L4_sheep.sorted.bam $5/RmM001_FKDN230042835-1A_HYNL7DSX3_L3_sheep.sorted.bam
date "+%F %T"; echo "Done merging bam files."
date "+%F %T"; echo "Indexing merged bam..."
samtools index -@ $1 $5/${4}_${7}.sorted.merge.bam
date "+%F %T"; echo "Done indexing merged bam."

# Mark duplicates
#date "+%F %T"; echo "Marking duplicate reads..."
gatk --java-options "-Xmx18G" MarkDuplicates \
 -I $5/${4}_${7}.sorted.merge.bam \
 -O $5/${4}_${7}.MD.bam \
 -M $5/${4}_${7}.MD.metrics.txt \
 --TAGGING_POLICY All \
 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 # 2500 for NovaSeq reads (uses a patterned flowcell).
MDBAM="$5/${4}_${7}.MD.bam"
date "+%F %T"; echo "Done marking duplicates."

# Filter bam file for high quality mapped reads
date "+%F %T"; echo "Filtering bam file; discarding reads with MQ<30..."
samtools view -@ $1 -q 30 -b --write-index -o $5/${4}_${7}.MD.MQ30.bam $MDBAM
MQ30BAM="$5/${4}_${7}.MD.MQ30.bam"
date "+%F %T"; echo "Done filtering bam file."

# Index and calculate stats
date "+%F %T"; echo "Calculating mapping statistics for MQ30 bam with samtools view and depth..."
samtools view -@ $1 -c $MQ30BAM | awk '{print "# combined reads " $1}' > ${MQ30BAM/%.bam/.mapping.results.txt}
samtools depth -a $MQ30BAM | awk '{sum+=$3;cnt++}END{print " Coverage " sum/cnt " Total mapped bp " sum}' >> ${MQ30BAM/%.bam/.mapping.results.txt}
date "+%F %T"; echo "Calculating mapping statistics for MQ30 bam with gatk CollectWgsMetrics..."
gatk --java-options "-Xmx90G" CollectWgsMetrics \
 -I $MQ30BAM \
 -O ${MQ30BAM/%.bam/.WgsMetrics.txt} \
 -R $6
date "+%F %T"; echo "Done calculating mapping statistics. Remember to remove intermediate bam files."

## Get stats for "raw" bam file
date "+%F %T"; echo "Calculating mapping statistics for MD bam with samtools view and depth..."
samtools view -@ $1 -c $MDBAM | awk '{print "# combined reads " $1}' > ${MDBAM/%.bam/.mapping.results.txt}
samtools depth -a $MDBAM | awk '{sum+=$3;cnt++}END{print " Coverage " sum/cnt " Total mapped bp " sum}' >> ${MDBAM/%.bam/.mapping.results.txt}
date "+%F %T"; echo "Calculating mapping statistics for MD bam with gatk CollectWgsMetrics..."
gatk --java-options "-Xmx90G" CollectWgsMetrics \
 -I $MDBAM \
 -O ${MDBAM/%.bam/.WgsMetrics.txt} \
 -R $6
date "+%F %T"; echo "Done calculating mapping statistics for MD bam. Remember to remove intermediate bam files."

## Remove intermediate files
rm RmM001_FKDN230042835-1A_HYNL7DSX3_L3_sheep.sorted.bam RmM001_FKDN230042835-1A_HNCK2DSX3_L4_sheep.sorted.bam $5/${4}_${7}.sorted.merge.bam 
