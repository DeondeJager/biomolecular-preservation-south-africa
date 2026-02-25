#!/bin/bash
#SBATCH --job-name=paleomix_ReedbuckSn_nc_host         # Job name
#SBATCH --mail-type=BEGIN,END,FAIL      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=  # Where to send mail	
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=10               # Number of CPU cores per task
#SBATCH --mem-per-cpu=10G                # Memory required per CPU
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/userid/palaeobovids/gargammel/paleomix/ReedbuckSn/paleomix_ReedbuckSn_nc_host_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/userid/palaeobovids/gargammel/paleomix/ReedbuckSn/paleomix_ReedbuckSn_nc_host_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script
hostname; date

echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules/activate conda environments
## Initiate conda environment on the compute node first
source /home/userid/.bashrc
### Activate package
conda deactivate
conda activate paleomix # This is version 1.3.9

# Navigate to working directory where yaml file is
cd /projects/lorenzen/people/userid/palaeobovids/gargammel/paleomix/ReedbuckSn/nc_host

# Execute paleomix
paleomix bam run paleomix_ReedbuckSn_nc_host.yaml --max-threads 10 --adapterremoval-max-threads 10 --bwa-max-threads 10

# Remove old summary file (otherwise it will be appended and not overwritten with the new data)
rm paleomix_ReedbuckSn_nc_host_summary.txt

# Combine results into summary table
awk 'BEGIN{second=0;one="";two="";three="";header="";line=""}{if($0!~/^#/ && $0!=""){if($0!~/^Target/){if($1!=one || $2!=two || $3!=three){if(header!="" && second==0){print "Target\tSample\tLibrary\t"header;second=1};if(line!=""){print one"\t"two"\t"three"\t"line;line=""}};if(second==0){split($0,h,"#");header=header h[2] "\t"}line=line $5 "\t";one=$1;two=$2;three=$3}}}END{print one"\t"two"\t"three"\t"line}' *summary | grep "\*\s*\*\s*" >> paleomix_ReedbuckSn_nc_host_summary.txt

