# Generate simulated ancient mitochondrial reads
 - We used the damage patterns from sample RmA009 (Cape grysbok) mapped to the Cape grysbok nuclear genome to generate the ancient DNA reads for all samples.
 - The script `1.1.1-prep_mapDamage.sh` generates the `Fragment_lengths.txt` file needed by Gargammel.
 - The script `1.1.2-palaeobovids_mito_gargammel.slurm` then uses gargammel to generate the damaged reads from the mitogenome of each target species.