#!/bin/bash

# Print compute node name and the date
hostname; date

# This line prints how many CPUs are being used
# It's a sanity check against your SLURM and software-specific settings
echo "Running job on $SLURM_CPUS_ON_NODE CPU cores"

# Load modules/activate conda environments
## Initiate conda environment on the compute node first
source /home/pzx702/.bashrc
### Activate package
conda deactivate
conda activate paleomix # This is version 1.3.9

# Navigate to working directory where yaml file is
cd /projects/lorenzen/people/pzx702/palaeobovids/paleomix/compMap_nc

# Execute paleomix
echo "$(date) Mapping sample RmA001..."
paleomix bam run paleomix_RmA001_compMap_nc.yaml --max-threads 10 --adapterremoval-max-threads 10 --bwa-max-threads 10 --temp-root ./temp.RmA001
echo "$(date) Done mapping sample RmA001."

