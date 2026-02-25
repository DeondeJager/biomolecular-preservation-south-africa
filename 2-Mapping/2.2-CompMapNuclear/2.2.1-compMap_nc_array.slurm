#!/bin/bash
#SBATCH --job-name=paleomix_compMap_nc         # Job name
#SBATCH --exclude=mjolnircomp05fl,mjolnircomp06fl # Exclude these two nodes, as they're a bit older and slower than the others.	
#SBATCH --array=1-55,69-161%20   # The range of indices for batch array (how many jobs)%No. of jobs to run in parallel.
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=10               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=12:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/userid/palaeobovids/paleomix/%x_%A.%a.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/userid/palaeobovids/paleomix/%x_%A.%a.err         # Standard error log

# Path to job files
echo "Processing input file compMap_nc_${SLURM_ARRAY_TASK_ID}.sh"
bash /home/userid/scripts/paleomix/palaeobovids/compMap_nc/job_files/compMap_nc_${SLURM_ARRAY_TASK_ID}.sh
wait
# Remove old summary file (otherwise it will be appended and not overwritten with the new data)

cd /projects/lorenzen/people/userid/palaeobovids/paleomix/compMap_nc
rm paleomix_compMap_nc.summary.txt

# Combine results across all samples into summary table
awk 'BEGIN{second=0;one="";two="";three="";header="";line=""}{if($0!~/^#/ && $0!=""){if($0!~/^Target/){if($1!=one || $2!=two || $3!=three){if(header!="" && second==0){print "Target\tSample\tLibrary\t"header;second=1};if(line!=""){print one"\t"two"\t"three"\t"line;line=""}};if(second==0){split($0,h,"#");header=header h[2] "\t"}line=line $5 "\t";one=$1;two=$2;three=$3}}}END{print one"\t"two"\t"three"\t"line}' *summary | grep "\*\s*\*\s*" >> paleomix_compMap_nc.summary.txt

