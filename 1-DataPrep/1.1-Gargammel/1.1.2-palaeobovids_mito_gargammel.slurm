#!/bin/bash
#SBATCH --job-name=gargammel_palaeobovids_mito         # Job name
#SBATCH --mail-type=BEGIN,END,FAIL      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dejager4@gmail.com  # Where to send mail	
#SBATCH --nodes=1                       # Run all processes on a single node
#SBATCH --ntasks=1                      # Run on a single task
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem-per-cpu=12G                # Memory required per CPU
#SBATCH --time=00:10:00                 # Time limit hrs:min:sec
#SBATCH --output=/projects/lorenzen/people/pzx702/palaeobovids/gargammel/gargammel_mito_%j.out        # Standard output log
#SBATCH --error=/projects/lorenzen/people/pzx702/palaeobovids/gargammel/gargammel_mito_%j.err         # Standard error log

# The %j in the --output & --error lines tells SLURM to substitute the job ID in the name of the output file
# See: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts#Sample_SLURM_Scripts 
# and https://help.rc.ufl.edu/doc/Annotated_SLURM_Script

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Activate conda environment (if module is not available)
#source /home/pzx702/.bashrc
#conda activate /projects/lorenzen/apps/conda/...

# Load modules and set up enviroment variables
module purge
module load gargammel/1.1.2
module load bbmap/39.17 

# Navigate to working directory
cd /projects/lorenzen/people/pzx702/palaeobovids/gargammel/data
DAMAGE=/projects/lorenzen/people/pzx702/palaeobovids/gargammel/paleomix/RmA009_SCR086/nuclear/SCR_086.GrysbokCape.mapDamage/SCR_086

# Insert commands here:
# Generate 150PE reads
## Cape buffalo mito (20X coverage, relative to nuclear data which is at 1X coverage). Can subsample data later for lower coverage if needed.
echo "Processing Cape buffalo..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o BuffaloCape/mito/BuffaloCape_mt.damaged BuffaloCape/mito/
### Prepend "MT" to read names:
rename.sh in=BuffaloCape/mito/BuffaloCape_mt.damaged_s1.fq.gz in2=BuffaloCape/mito/BuffaloCape_mt.damaged_s2.fq.gz out=BuffaloCape/mito/BuffaloCape_mt.damaged.rename_s1.fq.gz out2=BuffaloCape/mito/BuffaloCape_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
### Remove previous read files
rm BuffaloCape/mito/BuffaloCape_mt.damaged_s1.fq.gz BuffaloCape/mito/BuffaloCape_mt.damaged_s2.fq.gz
date

## Cape buffalo divergent mitogenome
echo "Processing Cape buffalo divergent..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged BuffaloCapeDiv/mito/
rename.sh in=BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged_s1.fq.gz in2=BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged_s2.fq.gz out=BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged.rename_s1.fq.gz out2=BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
rm BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged_s1.fq.gz BuffaloCapeDiv/mito/BuffaloCapeDiv_mt.damaged_s2.fq.gz
date

## Common eland
echo "Processing eland..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o ElandCommon/mito/ElandCommon_mt.damaged ElandCommon/mito/
rename.sh in=ElandCommon/mito/ElandCommon_mt.damaged_s1.fq.gz in2=ElandCommon/mito/ElandCommon_mt.damaged_s2.fq.gz out=ElandCommon/mito/ElandCommon_mt.damaged.rename_s1.fq.gz out2=ElandCommon/mito/ElandCommon_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
rm ElandCommon/mito/ElandCommon_mt.damaged_s1.fq.gz ElandCommon/mito/ElandCommon_mt.damaged_s2.fq.gz
date

## Cape grysbok
echo "Processing Cape grysbok..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o GrysbokCape/mito/GrysbokCape_mt.damaged GrysbokCape/mito/
rename.sh in=GrysbokCape/mito/GrysbokCape_mt.damaged_s1.fq.gz in2=GrysbokCape/mito/GrysbokCape_mt.damaged_s2.fq.gz out=GrysbokCape/mito/GrysbokCape_mt.damaged.rename_s1.fq.gz out2=GrysbokCape/mito/GrysbokCape_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
rm GrysbokCape/mito/GrysbokCape_mt.damaged_s1.fq.gz GrysbokCape/mito/GrysbokCape_mt.damaged_s2.fq.gz
date

## Mountain reedbuck
echo "Processing mountain reedbuck..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o ReedbuckMt/mito/ReedbuckMt_mt.damaged ReedbuckMt/mito/
rename.sh in=ReedbuckMt/mito/ReedbuckMt_mt.damaged_s1.fq.gz in2=ReedbuckMt/mito/ReedbuckMt_mt.damaged_s2.fq.gz out=ReedbuckMt/mito/ReedbuckMt_mt.damaged.rename_s1.fq.gz out2=ReedbuckMt/mito/ReedbuckMt_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
rm ReedbuckMt/mito/ReedbuckMt_mt.damaged_s1.fq.gz ReedbuckMt/mito/ReedbuckMt_mt.damaged_s2.fq.gz
date

## Southern reedbuck
echo "Processing southern reedbuck..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o ReedbuckSn/mito/ReedbuckSn_mt.damaged ReedbuckSn/mito/
rename.sh in=ReedbuckSn/mito/ReedbuckSn_mt.damaged_s1.fq.gz in2=ReedbuckSn/mito/ReedbuckSn_mt.damaged_s2.fq.gz out=ReedbuckSn/mito/ReedbuckSn_mt.damaged.rename_s1.fq.gz out2=ReedbuckSn/mito/ReedbuckSn_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
rm ReedbuckSn/mito/ReedbuckSn_mt.damaged_s1.fq.gz ReedbuckSn/mito/ReedbuckSn_mt.damaged_s2.fq.gz
date

## Grey rhebok
echo "Processing grey rhebok..."
gargammel -c 20 --comp 0,0,1 -f $DAMAGE/Fragment_lengths.txt -mapdamage $DAMAGE/misincorporation.txt single -rl 150 -o RhebokGrey/mito/RhebokGrey_mt.damaged RhebokGrey/mito/
rename.sh in=RhebokGrey/mito/RhebokGrey_mt.damaged_s1.fq.gz in2=RhebokGrey/mito/RhebokGrey_mt.damaged_s2.fq.gz out=RhebokGrey/mito/RhebokGrey_mt.damaged.rename_s1.fq.gz out2=RhebokGrey/mito/RhebokGrey_mt.damaged.rename_s2.fq.gz prefix=MT addprefix=t ow=t
rm RhebokGrey/mito/RhebokGrey_mt.damaged_s1.fq.gz RhebokGrey/mito/RhebokGrey_mt.damaged_s2.fq.gz
date
echo "Done."
