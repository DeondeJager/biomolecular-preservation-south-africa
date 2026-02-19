# Masking microbial-like regions

## Eland and Cape grysbok reference genomes
 - Because we altered the eland reference genome (see 1.2-MakeCapeGrysbokRef) and made the Cape grysbok reference genome from scratch, we could not use the microbial-like BED files from [Oskolkov et al. (2025)](https://doi.org/10.1093/gigascience/giaf108) to mask these reference genomes (or rather mask BAM files once reads were mapped to them).
 - Thus, we followed the MCWorkflow from Oskolkov et al. (2025) found [here](https://github.com/NikolayOskolkov/MCWorkflow) to identify microbial-like regions and produce the BED files for each of these two reference genomes.
 - This is implemented in the two `1.4.1-*` scripts.

## Masking microbial-like regions in BAM files
 - See the folder `2-Mapping/2.2-CompMapNuclear`.