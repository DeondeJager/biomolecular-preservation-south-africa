# Prep mapDamage output files for aDNA damage simulation with gargammel.

# We need to do is create a txt file with the proportion of each fragment length based on the mapdamage read lengths from the empirical data (the file lgdistribution.txt in the mapdamage output directory):.
# See: https://github.com/Mvwestbury/AncientGenomicsWorkshop25/tree/main

## Mito
cd /projects/lorenzen/people/userid/palaeobovids/gargammel/paleomix/RmA009_SCR086/mito/SCR_086.GrysbokCape.mapDamage/SCR_086
awk '/\+/{sum+=$3; count[$2]+=$3} END{for (i in count) print i"\t"count[i]/sum}' lgdistribution.txt > Fragment_lengths.txt

