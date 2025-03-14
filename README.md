This repository contains scripts is associated with the paper:  

## Effects of repeated freeze and thaw cycles on the stability of faecal microbiome composition.
Published in [*** XYZ ***](http...), on ...

[DOI: XYZ](http...)

#### Authors
**Matteo Sangermani**<sup>1,2</sup>, Indri Desiati<sup>3</sup>, Nicole Quattrini<sup>1</sup>, Guro F. Giskeødegård<sup>1,2</sup>

<h6>
<sup>1</sup> Dep. of Public Health and Nursing NTNU, Trondheim. 
<sup>2</sup> Dep of Surgery, St. Olavs University Hospital, Trondheim. 
<sup>3</sup> Dep. of Circulation and Medical Imaging, NTNU, Trondheim. 
</h6>

## Abstract of the paper
<h5>Background:
	Over the past two decades large population studies exploring the gut microbiome have resulted in the accumulation of collections of faecal samples in the freezers of laboratories. Data generated from previous microbiome measurements can be unsuitable for today’s research, due to the rapid improvement in quality and accuracy of sequencing technology, as well as more advanced bioinformatic tools. The amassed past faecal collections that lay unused could be remeasured to address new research questions, bypassing often lengthy procedure of recruiting and sampling a cohort. However, before reutilising past collections for newer microbiome analysis, it is necessary to dispel uncertainties on potential alterations the faecal microbiome composition may experience upon repeated freezing and thawing, thereby ensuring that sequencing faithfully represent the sample’s original gut microbiome composition.
</h5> 
<h5>Results:
	To this end, we investigated the effects that freezing and thawing (FT) and how this may alter the microbiome compositions. Specifically, we explored the effect FT cycle had on 16S rRNA sequencing and its ability to reproduce the gut microbiome composition consistently. Our study demonstrates that repeated FT cycles have minimal impact on the stability of the gut microbiome. Overall, changes due to FT cycles appear to be minor and overshadowed by interindividual effect. Repeated freezing and thawing did not compromise the ability to uniquely identify individuals using the gut microbiome, even after 6 FT cycles. Differential abundance analyses revealed that the first freezing and thawing had a few significantly altered microbial genera, when compared to the fresh, unfrozen material. Further, the two FT cycles after the first thawing showed no significant changes. 
</h5> 
<h5>Conclusions:
	Minor changes are detected during repeated freeze and thaw cycle, restricted to the first time freezing of fresh samples and after 4 FT cycles. These findings provide confidence for the potential to reuse stored faecal samples that have been thawed only once, which opens the possibility to reutilise past collections to conduct new microbiome studies.
</h5> 

## Description of the Scripts
This repository includes several scripts used in the study. Below is a brief overview of each:

- **Figure_1A_Metrics.R**: &nbsp;&nbsp;&nbsp;*Generates box plots of alpha and beta diversity metrics.*
- **Figure_1B_CoV_ICC.R**: &nbsp;&nbsp;&nbsp;*Analyzes the Coefficient of Variation (CoV) and Interclass correlation (ICC).*
- **Figure_2_Multivariate.R**: &nbsp;&nbsp;&nbsp;*Multivariate analisys of sample groups, via tSNE, PCA, PCoA and NMDS; generate many different subplot.*
- **Figure_3_MaAslin2.R**: &nbsp;&nbsp;&nbsp;*Perform differential abundance analysis using the MaAsLin2 method, and plot the results for significant taxa.*
- **SupFig_1_Temperature_vs_Time.R**: &nbsp;&nbsp;&nbsp;*Plot the temeperature versus time during the thawing of fecal samples.*
- **SupTable_1_Wilcox_v2.R**: &nbsp;&nbsp;&nbsp;*Perform feature wise Wilcoxon Rank sum tests (Mann-Whitney test).*
- **SupTable_ALDEx_glm.R**: &nbsp;&nbsp;&nbsp;*Perform differential abundance analysis using the ALDEx2, and plot the results for significant taxa.*


Each script is designed to be run independently and they can generate several subfigures of tables presented in the published paper. Moreover, we have two repository-like scripts:

- **Plotting_Repository.R**: &nbsp;&nbsp;&nbsp;*This script gather all the function related to plotting and layouts, to allow to apply a unified theme across all scripts*
- **PruneWrangle_Data.R**: &nbsp;&nbsp;&nbsp;*Mix of functions that re-format the microbiome feature tables or perform complex, repeated calculation.* 


## Description of Data Tables:
Table used for analysis are included in this repository. Folder *Data_Seq* contains the count data and taxonomic information of microbiome features obtained with 16S rRNA sequencing. Below is an overview of the files:
- *Data_Seq*
	- *ASV_freeze.tsv*: read count table of the cohort, where each column is a microbial ASV.
	- *TAXO_freeze.tsv*: taxonomic information for each ASV present in the cohort.
	- *META_freeze*: The metadata information of the cohort, importantly the variables "Sample_ID", "Exp_FTcycle",	"Exp_Name".
	- *RefList_Class.txt* and *RefList_Phyla.txt*: a simple list of all the classes and phyla present in the cohort listed by abundance and with a reference color to use to coherently highlight features in plots.
	- *Drop_Features_Named.txt*: taxonomy in 16S rRNS sequencing is often incomplete, specifically at species rank. These is a list of names (such as "unkown", "unclassified genome", etc) that should ignore as "species" and such ASV will be re-labelled in the scripts.
	- *LUT_Microbes*: LUT color table; each class in the cohort is assigned a range of tonalities from the same base color.
	- *Temperature_Record_FT_fecal.xlsx* is the temperature measurements of the stool thawing
	
## Main Figures
The following figures showcase the results:

| ![Figure 1](Final_Figures/Figure_1.png) | ![Figure 2](Final_Figures/Figure_2.png) | ![Figure 3](Final_Figures/Figure_3.png) |
|-----------------------------------------|-----------------------------------------|-----------------------------------------|
| Figure 1                                | Figure 2                                | Figure 3                                |

## Usage
To run these scripts, the following R packages are required:
#### Visualization:
- `ggplot2`
- `cowplot`
- `tiff`
- `ggrepel`
- `ggtext`
- `showtext`
#### Data Manipulation:
- `readr`
- `stringr`
- `tidyr`
- `tibble`
- `dplyr`
#### Statistical Analysis:
- `Rtsne`
- `mixOmics`
- `nlme`
- `nortest`
- `rstatisitx`
- `vegan`
- `Maaslin2`
- `ALDEx2`

#### Setup Instructions:
In each script, change the variable `cDir_Rscript` to the "absolute path" of the repository on the machine running the code. All other paths are relative to this and follow Unix/MacOS format.












