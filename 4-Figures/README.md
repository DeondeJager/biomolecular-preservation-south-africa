# Code for producing figures

## Figure 1
 - The R code `4.1-plot_Fig1.R` uses the files `Table S1 - Specimen metadata.csv` and `Sites_GPS.csv` to produce the individual components of Figure 1.
 - The components were then combined in [Inkscape](https://inkscape.org/), where further cosmetic tweaks were made.

## Figure 2
 - The R code `4.2-plot_preservation_Fig2.R` uses the file `Table S3 - paleomix nuclear & collagen.csv` to produce Figure 2.
 - Further cosmetic tweaks were made in Inkscape, with the most important listed below:
	- Add a non-transparent border to all green triangles (using "select same fill" by right-clicking on a triangle) by setting the stroke paint to `2a9e72a6`.
	- Change the transparency of all grey triangles by changing the fill colour to `a9a9a91a` and adding a non-transparent border by setting the stroke paint to `969696ff`.

## Figure 3
 - Fig. 3 is produced by the script `compare_SCR_BEST.R` in folder `2-Mapping\2.4-SCR_BEST_Comparison`; specifically section 3.3 of the script.

## Figure S1
 - Produced by the R code `4.3-plot_mito_panel_FigS1.R`, which uses the file `Table S3 - paleomix nuclear & collagen.csv`.

## Figure S2
 - The R code `4.4-plot_damageProfiler_FigS2.R` uses the `misincorporation.txt` and JSON output files from DamageProfiler in the `damageProfiler` folder here and the file `Table S3 - paleomix nuclear & collagen.csv` to produce Figure S2.

## Figure S3-S5
 - The R code `4.5-plot_preservation_FigS3toS5.R` produces Figures S3 to S5, using the file `Table S3 - paleomix nuclear & collagen.csv`.
 - The same cosmetic tweaks were made to these plots as for Figure 2, in Inkscape.