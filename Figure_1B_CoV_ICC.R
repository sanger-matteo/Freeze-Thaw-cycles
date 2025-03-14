# +++++ DESCRIPTION - CoV and ICC analysis +++++++++++++++++++++++++++++++++++++
#
# This script summarise the results from the analysis of the Coefficient of
# Variation (CoV) of microbiome as well as the InterClass Correlation (ICC). 
# This can be done at any specific the taxonomic rank (variable "aggr_at_rank").
# 
# DEF: ICC quantifies the degree to which individuals with a fixed degree of 
# relatedness (assumption) resemble each other in terms of a quantitative trait.
#
# When run the script generate two charts:
# - Fig A - Coefficient of variation (CoV)          --> hnd_A
# - Fig B - Interclass correlation (ICC)            --> hnd_B
#
# For the published article we run the script to generate three sets of plots, 
# each analysing one specific taxonomic rank:
# Genus:     for Figure 1E to 1F 
# Family:    for Suppl. Figure 3C to 3D   
# Class:     for Suppl. Figure 3A to 3B   
#
# The taxa are colored consistently across the paper, using a file that ranks 
# and assign specific color tone based on taxonomic "Class". This is a reference
# for the Class ranking order as well.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# --- Load Libraries and functions
library(ggplot2)
library(cowplot)
library(stringr)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(nlme)
library(nortest)

# Define the ICC function: very simple, no package needed
icc.lme <- function(lmeobj){
  varcorr <- VarCorr(lmeobj)
  return( as.numeric(varcorr[[1]]) / (as.numeric(varcorr[[1]])+as.numeric(varcorr[[2]]) ) )
}  
# Define function to perform division by zero and return zero (not NaN)
"%/0%" <- function(x,y) { res <- x/y;   res[is.na(res)] <-0;  return(res) }

cDir_Rscript <- setwd("/Users/sm/Documents/R/FT_Cycle")
path_2_data  <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/" )
path_saveFig <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/Figures/" )
save_subFigs <- F

# Import custom functions, specifically contains those to calculated microbiome diversity metrics
source( paste0(path_2_data, "Plotting_Repository.R") )
source( paste0(path_2_data, "PruneWrangle_Data.R") )



# ----- CoV --------------------------------------------------------------------------------------

# +++ Load DATA tables +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load ASV table(s) and Metadata
ASV  <- read.table( paste0( path_2_data,"Data_Seq/ASV_freeze.tsv"), header=T, sep="\t")
TAXO <- read.table( paste0( path_2_data,"Data_Seq/TAXO_freeze.tsv"), header=T, sep="\t")
MetaData <- read.table( paste0(path_2_data, "Data_Seq/META_freeze.txt"), header=T, sep="\t")

# Subset Metadata to only the key grouping variables
sub_MetaData <- MetaData[ , c("Exp_Name", "Exp_FTcycle")]


# +++ Prune Microbiome Dataset +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++ Aggregate ASVs at rank X:
# Import user defined aggregation function
show_topFeat <- 25
filter_features <- T
aggr_at_rank <- "Family"
returned_list <- aggregate_at_Rank( ASV, TAXO, aggr_at_rank, F, F, F )
reduc_ASVs <- returned_list[[1]]
reduc_Taxo <- returned_list[[2]]


# +++ Filter features (based on RA):
if( filter_features==TRUE ){
  # - thres_abundance:     min relative abundance of a feature
  # - thres_freqCohort:    min frequency of a feature in the cohort
  # Generally, use Rel.Abu. less than 0.03%, and freq. greater than 10% of the cohort
  # Also, adjust the TAXO table to contain only features present in the cohort.
  thres_abundance  <- 0.003 /100
  thres_freqCohort <- 0.01  /100
  taxa_RA_Freq <- data.frame( Taxa  =colnames(reduc_ASVs[,-1]) ,
                              N_Obs_above_RA_thres=colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) ,
                              Percent_of_cohort   =colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] *100,
                              Keep_taxa           =( colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] ) >= thres_freqCohort )
  rownames(taxa_RA_Freq) <- seq(1,nrow(taxa_RA_Freq))
  mask_filterASV <- as.vector( ( colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] ) >= thres_freqCohort )
  sprintf(" --> # taxon to be remove=%d", sum(!mask_filterASV) )
  reduc_ASVs <- reduc_ASVs[ , c(T,mask_filterASV) ]
  # Mask and reduce Taxo to the essential taxonomy, that used in reduc_ASVs table
  reduc_Taxo <- reduc_Taxo[ mask_filterASV , ]
}

# > > > Normalize reads to Relative abundance
reduc_ASVs[,-1] <- (reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) *100


# Make ASV table fully numeric and with only ASV features
rownames(reduc_ASVs) <- reduc_ASVs$Sample_ID
reduc_ASVs <- reduc_ASVs[,-1]
# Replace ASV_ID colnames with taxonomic names
colnames(reduc_ASVs) <- reduc_Taxo[[aggr_at_rank]]

# Drop the "unkown" features, which are taxonomically unclassified.
# There is only 1 in the top 50 most abundant Genera. Removing it is not a big inaccuracy
idx <- grep( "unknown_", colnames(reduc_ASVs) )
if ( length(idx)!=0 ){
  reduc_ASVs <- reduc_ASVs[ , -idx ]
  reduc_Taxo <- reduc_Taxo[ -idx, ]
}



# +++ Calculate CoV ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate using user defined different pairs of Freeze-Thaw timepoint (TP_2)
# and a compared to specific reference timepoint (TP_ref)
TP_ref <- "C0"
TP_2   <- c("C1","C2","C3","C6")
cov_val <- data.frame()
for (ii in seq(1,length(TP_2))) {
  selected_Samples_ASVs     <- reduc_ASVs[sub_MetaData$Exp_FTcycle %in% c(TP_ref,TP_2[ii]) , ]
  selected_Samples_MetaData <- sub_MetaData[sub_MetaData$Exp_FTcycle %in% c(TP_ref,TP_2[ii]) , ]
  
  temp <- cbind(selected_Samples_MetaData,  selected_Samples_ASVs) %>% 
              group_by(Exp_Name) %>%
              reframe( across( colnames(selected_Samples_ASVs), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
              as.data.frame()
  temp <- replace(temp, is.na(temp), 0)
  temp <- temp %>% mutate( ., FT_pair=paste0( TP_ref, "_", TP_2[ii]), .after=Exp_Name )
  cov_val <- rbind( cov_val, temp )
}

dsPlot <- cov_val %>% pivot_longer( cols     =-c(Exp_Name, FT_pair),  # exclude specific column
                                    names_to =aggr_at_rank,           # selected column for a new variable ("groupings")
                                    values_to="CoV")                  # values all into one new variables


# +++ Add color information ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import color LUT for microbiome ranks and set order to match color palette
upper_rank <- "Class"
ref_ranks <- read.table( paste0( path_2_data, "Data_Seq/RefList_",upper_rank,".txt"), header=T, sep="\t")
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")

# Merge taxonomic and color information to dsPlot
dsPlot <- left_join( dsPlot, reduc_Taxo[, c( "ASV_ID","Phyla","Class","Family","Genus","Species")], by=aggr_at_rank )
dsPlot <- left_join( dsPlot, ref_ranks, by=upper_rank )

# Order the features (rank = aggr_at_rank) by abundance
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
# Factorize the aggr_at_rank (e.g., Genus); then subset to take the "top" taxa
dsPlot[[aggr_at_rank]] <- factor( dsPlot[[aggr_at_rank]], levels=rev( ordered_Features ) )
dsPlot <- dsPlot[ dsPlot[[aggr_at_rank]] %in% ordered_Features[1:show_topFeat], ]
# Set up the palette colors, based on upper_rank
uniq_aggr  <- distinct(dsPlot, .data[[aggr_at_rank]], .keep_all=T)
my_palette <- setNames( uniq_aggr$Colour, uniq_aggr[[upper_rank]] )

# Calculate the average CoV to display on top of charts
average_CoV  <- dsPlot %>%
                group_by(FT_pair) %>%
                summarize(avg_CoV=mean(CoV))


# +++ Plotting CoV +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hnd_A <- ggplot( dsPlot, aes(x=.data[[aggr_at_rank]], y=CoV)) +  
          geom_boxplot( aes(fill=.data[[upper_rank]]), color="#777777", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
          geom_jitter(color="#555555", alpha=0.5, size=0.6, width=0.2) +
          coord_flip() + 
          ylim(-0.1, 2.01) + 
          labs(fill="Class:") +
          xlab("") + ylab("CoV") +
          Theme_BoxPlot +
          scale_colour_manual(values=my_palette) +
          scale_fill_manual(  values=my_palette) +
          facet_wrap( ~FT_pair, nrow=1, scales="free_x", strip.position="top",
                      labeller=labeller(FT_pair=function(x) {
                        avg_CoV <- round(average_CoV$avg_CoV[average_CoV$FT_pair==x], digits=3)
                        label <- sprintf( "%s \n(mean: %.3f)", x, avg_CoV ) 
                        return(label)
                      })
          )
hnd_A



# ----- Control CoV -------------------------------------------------------------------------------

# Load and include the REPEATED MEASUREMENTS control in the CoV plot
# +++ Prepare QC data ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load ASV table(s) and Metadata
ASV  <- read_tsv( "/Users/mattesa/ZenBook/R/16SMT_v1/DataTables/Reg_All/comb_ASV_Count.txt",
                  skip_empty_rows=TRUE ,show_col_types=FALSE,
                  guess_max =min( 10 , 300 ) )
TAXO <- read_tsv( "/Users/mattesa/ZenBook/R/16SMT_v1/DataTables/Reg_All/comb_TAXO.txt",
                  skip_empty_rows=TRUE ,show_col_types=FALSE,
                  guess_max =min( 10 , 300 ) ) 
MetaData <- read_tsv( "/Users/mattesa/ZenBook/R/16SMT_v1/DataTables/comb_Metadata.txt",
                      skip_empty_rows=TRUE ,show_col_types=FALSE,
                      guess_max =min( 10 , 300 ) ) 

# +++ Prune QC Dataset +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove control samples, with too low read numbers
idx <- ASV$Sample_ID[ ASV$Sample_TotReads < 1000 ]
ASV <- ASV[ ! ASV$Sample_ID %in% idx , ]
MetaData <- MetaData[ ! MetaData$Sample_ID %in% idx , ]

# > > > NOTE: Normalize reads to Relative abundance as we did for the major dataset
relative_features <- T
source( paste0( path_2_data, "Prune_Ctrl_Data.R") ) 


# Take only features EXP2 datasets
mask <- str_detect( MetaData$External_ID , "Exp2_DNA_" )
mask <- mask + str_detect( MetaData$External_ID , "A09" )    # A01, A06, A09, A16, A20
mask <- mask + str_detect( MetaData$External_ID , "C[123]" )     # C1 to C6    
mask <- mask == 3
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]

sub_MetaData <- MetaData

# +++ Calculate CoV for Ctrlt ++++++++++++++++++++++++++++++++++++++++++++++++++
# Group by individual and calculate the coefficient of variation
cov_val <- cbind(sub_MetaData,  reduc_ASVs) %>% 
            group_by(Exp_Name) %>%
            reframe(  across( colnames(reduc_ASVs), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
            as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)

dsCTRL <- cov_val %>% pivot_longer( cols     =-Exp_Name,     # exclude specific column
                                    names_to ="Feature",     # selected column for a new variable ("groupings")
                                    values_to="CoV")         # values all into one new variables

# Order the Features to be plotted by the relative frequency in which they appear
# by using the previously plotted order: ordered_Features
dsCTRL <- dsCTRL[ !is.na(dsCTRL$Feature) , ]
dsCTRL <- dsCTRL[ match(ordered_Features,dsCTRL$Feature) , ]
# Factorize the "Feature" and take only the # of feature to show that matches original plot
dsCTRL$Feature <- factor( dsCTRL$Feature, levels=rev( ordered_Features) )         
dsCTRL <- dsCTRL[ dsCTRL$Feature %in% ordered_Features[1:show_topFeat], ]

# +++ OVERLAY QC on hnd_A plot +++++++++++++++++++++++++++++++++++++++++++++++++
hnd_A <- hnd_A + #ggplot(NULL) + 
          geom_jitter(  data=dsCTRL, aes(x=Feature, y=CoV) ,color="#8a4949", fill="#963939", shape=23, size=1.5, width=0) +
          coord_flip() + 
          ylim( -0.1, 2.01)+
          Theme_BoxPlot +
          scale_colour_manual( values=my_palette) +
          scale_fill_manual( values=my_palette)

hnd_A



# ----- ICC -------------------------------------------------------------------------------

# +++ Prepare QC data ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load ASV table(s) and Metadata
ASV  <- read.table( paste0( path_2_data,"Data_Seq/ASV_freeze.tsv"), header=T, sep="\t")
TAXO <- read.table( paste0( path_2_data,"Data_Seq/TAXO_freeze.tsv"), header=T, sep="\t")
MetaData <- read.table( paste0(path_2_data, "Data_Seq/META_freeze.txt"), header=T, sep="\t")

# +++ Aggregate ASVs at rank X:
returned_list <- aggregate_at_Rank( ASV, TAXO, aggr_at_rank, F, F, F )
reduc_ASVs <- returned_list[[1]]
reduc_Taxo <- returned_list[[2]]

# +++ Filter features (based on RA):
if(filter_features==TRUE){
  # - thres_abundance:     min relative abundance of a feature
  # - thres_freqCohort:    min frequency of a feature in the cohort
  # Generally, use Rel.Abu. less than 0.03%, and freq. greater than 10% of the cohort
  # Also, adjust the TAXO table to contain only features present in the cohort.
  thres_abundance  <- 0.003 /100
  thres_freqCohort <- 0.01  /100
  taxa_RA_Freq <- data.frame( Taxa  =colnames(reduc_ASVs[,-1]) ,
                              N_Obs_above_RA_thres=colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) ,
                              Percent_of_cohort   =colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] *100,
                              Keep_taxa           =( colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] ) >= thres_freqCohort )
  rownames(taxa_RA_Freq) <- seq(1,nrow(taxa_RA_Freq))
  mask_filterASV <- as.vector( ( colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] ) >= thres_freqCohort )
  sprintf(" --> # taxon to be remove=%d", sum(!mask_filterASV) )
  reduc_ASVs <- reduc_ASVs[ , c(T,mask_filterASV) ]
  # Mask and reduce Taxo to the essential taxonomy, that used in reduc_ASVs table
  reduc_Taxo <- reduc_Taxo[ mask_filterASV , ]
}

# > > > Normalize reads to Relative abundance
reduc_ASVs[,-1] <- (reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) *100

# Make ASV table fully numeric and with only ASV features
rownames(reduc_ASVs) <- reduc_ASVs$Sample_ID
reduc_ASVs <- reduc_ASVs[,-1]
# Replace ASV_ID colnames with taxonomic names
colnames(reduc_ASVs) <- reduc_Taxo[[aggr_at_rank]]

# Drop the "unkown" features, which have unclassified taxonomy.
# There is only 1 in the top 50 most abundant Genera. Removing it is not a big inaccuracy
idx <- grep( "unknown_", colnames(reduc_ASVs) )
if ( length(idx)!=0 ){
  reduc_ASVs <- reduc_ASVs[ , -idx ]
  reduc_Taxo <- reduc_Taxo[ -idx, ]
}


# +++ Calcaulate ICC +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate the ICC, including specified number of timepoints
mask <- str_detect( MetaData$External_ID , "C[0123456]" )     # C1 to C6
MetaData <- MetaData[ mask, ]
Data <- reduc_ASVs[ mask, ]
Data <- Data[ , colSums(Data) != 0 ]    # Remove empty features

Individuals <- MetaData$Exp_Name
n_features  <- dim(Data)[2]
ICC_res     <- vector('numeric',n_features)
P_res       <- vector('numeric',n_features)


# Run the ICC aalgorithm via LME
for (ii in seq(1:dim(Data)[2])) {
  y_Feature <- Data[ , ii]
  resultat  <- lme( y_Feature ~1, random=~1 | as.factor(Individuals) , na.action=na.omit ) 
  resultatY <- summary(resultat)
  temp <- as.data.frame(resultatY[["tTable"]])
  ICC_res[ii] <- icc.lme(resultat)
  P_res[ii]   <- temp[["p-value"]]
}
dsPlot_ICC=data.frame( Feature_Name=colnames(Data), ICC=ICC_res, P.value=P_res)
avg_icc=mean(dsPlot_ICC$ICC)
colnames(dsPlot_ICC)[1] <- aggr_at_rank

# Check if LME was significance and add "*" symbol that we can add in the plot
dsPlot_ICC$Significance <- ""
dsPlot_ICC[ dsPlot_ICC$P.value < 0.05, "Significance" ] <- "*"



# +++ Adding color and formatting ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Join the taxonomic and color information to dsPlot
dsPlot_ICC <- left_join( dsPlot_ICC, reduc_Taxo[, c( "ASV_ID","Phyla","Class","Family","Genus","Species")], by=aggr_at_rank )
dsPlot_ICC <- left_join( dsPlot_ICC, ref_ranks, by=upper_rank )
# Factorize the aggr_at_rank (e.g., Genus); then subset to take the "top" taxa
dsPlot_ICC[[aggr_at_rank]] <- factor( dsPlot_ICC[[aggr_at_rank]], levels=rev( ordered_Features ) )
dsPlot_ICC <- dsPlot_ICC[ dsPlot_ICC[[aggr_at_rank]] %in% ordered_Features[1:show_topFeat], ]
# NOTE: We use same color palette as define above for CoV, so no need to redefine it


# +++ Plotting ICC +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
hnd_B <- ggplot( dsPlot_ICC, aes( x=.data[[aggr_at_rank]], y=ICC) )  +
          geom_point( aes(fill=.data[[upper_rank]]), color="#666666", size=3, alpha=0.8, shape=23) +  
          geom_hline( yintercept=c(0.0,0.8), linetype="dashed", alpha=0.0) +
          scale_fill_manual(values=my_palette) +
          
          # geom_text( aes( x=.data[[aggr_at_rank]], y=1.07, label=Significance),
          #            color="#400000", size=6, vjust=0.8) +
  
          coord_flip() + ylim(0,1.08) + 
          labs(x="", colour="Metabolite Group:") +
          Theme_BoxPlot + 
          theme( axis.text.y=element_blank() ) +
          scale_y_continuous(n.breaks=6) 
hnd_B



# --- COMBINE sub-Figures --------------------------------------------------------------------------
# Place sub-plots together into one figure
ggdraw() +
  draw_plot(hnd_A, x=0, y=0, width=0.81, height =.96) +
  draw_plot(hnd_B, x=0.82, y=0, width=0.18, height =.910) +      # .910  .885  .845
  draw_plot_label(label=c("C","D"), size=17, x=c(0.07,0.84), y=c(.98,.98)) 

ggsave( paste0(path_saveFig, "Figure_1EF_",aggr_at_rank,".tiff"), width=10, height=6, dpi=300)
  

# Figure size and positions to use when plotting different taxonomic rank
#           height:       Fig label:             Save to file name:
# Genus:      9           label=c("C"),          Figure_1EF_          Genus
# Family:     6           label=c("C","D"),      Suppl_Figure_3CD_    Family
# Class:      3.5         label=c("A","B"),      Suppl_Figure_3AB_    Class














