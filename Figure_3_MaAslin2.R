# +++++ DESCRIPTION --- MaAsLin2.R +++++++++++++++++++++++++++++++++++++++++++++
#
# Differential Abundance Analysis using the algorithm:
#  MaAsLin2 - Microbiome Multivariable Association with Linear Models.
#             https://github.com/biobakery/biobakery/wiki/maaslin2
#
# We run in two mode for the article, changing the variable:
# - ref_TP :    either with value  "C[0123456]"  or  "C[123456]"  
# The MaAsLin algorithm uses the first level group as the reference to which it 
# compares all others. Using "C[0123456]" sets C0 as the reference time-point, 
# whereas using "C[123456]" set C1 as the reference time-point.
#
# The results of the MaAsLin2 are save in a .txt table to which we added 
# important cohort variables and model information. The same table is then used
# by the script to generate 3 subfigures for the main Figure 3 of the article:
# - hnd_2Ys  -  plot the amount of significant changes in taxa across FT cycles
#               The amount is plotted on with two reference axis: one reporting
#               the number of taxa found significantly changed at eac FT cycle ;
#               the proportion of all taxa found significantly changed.
# - hnd_Log2  - Chart the effect size in log2 scale of significantly changed 
#               taxa splitting the findings by FT cycle and placing them side
#               by side.
#               Significantly changed taxa (p<0.05) are plotted with full color
# - hnd_Percent - Same as above, but instead of log2 scale, data is show as 
#                 Percept change in respect to the reference time-point
#
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# --- Load Libraries, functions and data tables
library(ggplot2)
library(ggpattern)
library(cowplot)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(Maaslin2)

cDir_Rscript <- setwd("/Users/sm/Documents/R/FT_Cycle")
path_2_data  <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/" )
path_saveFig <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/Figures/" )
save_subFigs <- F

# Import custom functions, specifically contains those to calculated microbiome diversity metrics
source( paste0(path_2_data, "Plotting_Repository.R") )
source( paste0(path_2_data, "PruneWrangle_Data.R") )

# Load ASV table(s) and Metadata
ASV  <- read.table( paste0( path_2_data,"Data_Seq/ASV_freeze.tsv"), header=T, sep="\t")
TAXO <- read.table( paste0( path_2_data,"Data_Seq/TAXO_freeze.tsv"), header=T, sep="\t")
MetaData <- read.table( paste0(path_2_data, "Data_Seq/META_freeze.txt"), header=T, sep="\t")



# --- Prune Microbiome Dataset -------------------------------------------------
# +++ Aggregate ASVs at rank X:
# Import user defined aggregation function
aggr_at_rank <- "Genus"
returned_list <- aggregate_at_Rank( ASV, TAXO, aggr_at_rank, F, F, F )
reduc_ASVs <- returned_list[[1]]
reduc_Taxo <- returned_list[[2]]

# Make ASV table fully numeric and with only ASV features
rownames(reduc_ASVs) <- reduc_ASVs$Sample_ID
reduc_ASVs <- reduc_ASVs[,-1]
# Replace ASV_ID colnames with taxonomic names
colnames(reduc_ASVs) <- reduc_Taxo[[aggr_at_rank]]

# Include only specified number of timepoints
# NOTE: later we factorize, and thus make the "earliest" FT cycle the reference 
#       point for the Maaslin GLM; for our interest it will be either C0 or C1.
ref_TP <- "C[0123456]"       # "C[123456]"  or "C[0123456]"
mask <- stringr::str_detect( MetaData$External_ID , ref_TP )
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]
reduc_ASVs <- reduc_ASVs[ , colSums(reduc_ASVs) != 0 ]    # Remove empty features

# > > > Normalize reads to Relative abundance
temp <- (reduc_ASVs/rowSums(reduc_ASVs)) *100
rel_abundance <- data.frame( "RA_mean" =sapply( temp, function(x) mean(x) ) ,
                             "RA_stdev"=sapply( temp, function(x) sd(x) )   )
rel_abundance$Genus <- rownames(rel_abundance)
rownames(rel_abundance) <-  seq(1,nrow(rel_abundance))


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MaAsLin2 requires a strict rownames match between the feature table and metadata
sub_Meta <- MetaData[ , c("Sample_ID","Exp_FTcycle","Exp_Name")]
# MaAsLin2 requires that we factorize the categorical variables
sub_Meta$Exp_FTcycle <- factor( sub_Meta$Exp_FTcycle, levels=unique(sub_Meta$Exp_FTcycle) )
sub_Meta$Exp_Name    <- factor( sub_Meta$Exp_Name   , levels=unique(sub_Meta$Exp_Name) )
rownames(sub_Meta) <- sub_Meta$Sample_ID
fix_eff  <- "Exp_FTcycle"
rand_eff <- "Exp_Name"

# --- Run MaAsLin2 algorithm - using the recommended parameters that fit small cohort sizes, like ours
fit_data <- Maaslin2( reduc_ASVs, sub_Meta, paste0( path_2_data, "Output_MaAsLin2"),
                      fixed_effects  = fix_eff,
                      random_effects = rand_eff,
                      analysis_method= "LM",          # LM , CPLM (takes a lot of time) or NEGBIN
                      normalization  = "TSS",         # TSS, CLR or CSS
                      transform      = "LOG",         # LOG
                      correction     = "BH",
                      min_abundance  = 0.0,           # No filtering of features, we do it ...
                      min_prevalence = 0.0,           # ... in upper part of code
                      cores=5, 
                      standardize=FALSE, plot_scatter=FALSE ) 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# NOTE: In MaAsLin2 outupt 
# P.value - is the un-adjusted p-value
# Q.value - is the adjusted p.value
# coef    - "coefficient" (coef) is the effect size, which quantifies the
#           relationship between the predictor variable (FT-cycle) and the  
#           response variable for each feature.
#           A coefficient of 0.1 (log10 scale) means that a unit increase in 
#           the predictor (FT Cycle) is associated with a 10% increase in the 
#           relative abundance of the feature.
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# Import the just calculated results of MaAsLin2, format and save as result table
res_Maaslin <- as.data.frame( read_tsv( paste0( path_2_data, "Output_MaAsLin2/all_results.tsv" )) )

colnames(res_Maaslin) <- c( aggr_at_rank, "Fixed_effect", "FixEff_Value", "Coeff", "StDevErr", "N", "N.not.0", "P.value", "Q.value")

if(is.null(rand_eff)){   res_Maaslin <- res_Maaslin %>% mutate( "Random_effect"=NA,       .before="Coeff" )
}else{                   res_Maaslin <- res_Maaslin %>% mutate( "Random_effect"=rand_eff, .before="Coeff" ) }
res_Maaslin <- res_Maaslin %>% mutate( "Formula"=paste(fix_eff,collapse=" + "), .after="Random_effect" )
# Add taxonomic information
res_Maaslin <- merge( reduc_Taxo[,c("Phyla","Class",aggr_at_rank)], res_Maaslin, by=aggr_at_rank  )
res_Maaslin <- res_Maaslin %>% relocate( all_of(aggr_at_rank), .after=3 )

res_Maaslin <- res_Maaslin[ order(res_Maaslin[[aggr_at_rank]], res_Maaslin$FixEff_Value ) , ]
res_Maaslin$Significance <- ""
res_Maaslin[ res_Maaslin$Q.value < 0.05, "Significance" ] <- "*"
res_Maaslin[ res_Maaslin$Q.value < 0.01, "Significance" ] <- "**"

write.table( res_Maaslin, paste0(path_2_data, "Res_Tables/Sig_Table_MaAsLin2_",aggr_at_rank,".txt"),
             row.names=FALSE, col.names=TRUE, sep='\t' )



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 

# --- Plot Effect Size ----------------------------------------------------------------
# Import data for plotting significant changes in taxa over time
aggr_at_rank <- "Genus"
display_rank <- "Class"
res_Maaslin  <- read.table( paste0(path_2_data, "Res_Tables/Sig_Table_MaAsLin2_",aggr_at_rank,".txt"), header=T, sep="\t")

# Subset long format table
df_Long <- res_Maaslin[ c(display_rank,aggr_at_rank,"FixEff_Value","Coeff","Q.value") ]
df_Long$Time <- sapply( df_Long$FixEff_Value, function(x) as.numeric(strsplit(x,"C")[[1]][2]) )
df_Long <- relocate( df_Long, "Time", .after=3 )

# Find significant
sign_thres <- 0.05
sig_levels <- unique(df_Long[ df_Long$Q.value < sign_thres, aggr_at_rank ])
dsPlot <- df_Long[ df_Long[[aggr_at_rank]] %in% sig_levels,  ]
# Filter out Classes that have significance
dsPlot$Significance <- ""
dsPlot[ dsPlot$Q.value < sign_thres, "Significance" ] <- "*"

if( ref_TP == "C[123456]") {
  dsPlot <- dsPlot[ dsPlot$FixEff_Value %in% c("C4","C5","C6"), ]
}

# # Transform coefficients to proportional changes
dsPlot$Percent   <- (2^dsPlot$Coeff -1) *100        # Convert log2 to percentage change
# Add RA mean and SD values for each genera (aggr_at_rank)
dsPlot <- left_join( dsPlot, rel_abundance, by=aggr_at_rank )


# +++ Add color information +++++++++
ref_ranks <- read.table( paste0( path_2_data, "Data_Seq/RefList_",display_rank,".txt"), header=T, sep="\t")
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(display_rank, "Colour")
# Merge taxonomic and color information to dsPlot
dsPlot <- left_join( dsPlot, ref_ranks, by=display_rank )
# Order the features by abundance and factorize aggr_at_rank (e.g., Genus).
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
dsPlot[[aggr_at_rank]] <- factor( dsPlot[[aggr_at_rank]], levels= rev(ordered_Features) )
# Set up the palette colors, based on display_rank
my_palette <- setNames( ref_ranks$Colour, ref_ranks[[display_rank]] )


# --- Plot effect size, dividing into panels for the different FT cycles
# --- Plot the Log2 change
res_scale <- "Coeff"
yMax <- abs( ceiling( max( dsPlot[[res_scale]]) ))
yMin <- -12   # floor(   min( dsPlot[[res_scale]]) )
hnd_Log2 <- ggplot( dsPlot, aes(x=.data[[aggr_at_rank]], y=.data[[res_scale]], 
                                fill=.data[[display_rank]], pattern=Q.value <sign_thres ) ) +
          geom_bar_pattern(stat="identity",  position="dodge", alpha=0.7, linewidth=0.35,
                           color="#777777", pattern_color="#DDDDDD", pattern_fill="#EEEEEE",
                           pattern_angle=45, pattern_density=0.2,
                           pattern_spacing=0.04, pattern_key_scale_factor=0.6) + 
          geom_hline( aes(yintercept=0), size=1, linetype="solid", color="#544444") +
          coord_flip() +
          facet_wrap( ~FixEff_Value, scales="fixed", nrow=1 ) + 
          # scale_y_continuous( breaks=seq(yMin,yMax,by=2), labels=seq(yMin,yMax,by=2) ) +
          scale_fill_manual( values = my_palette) +
          scale_pattern_manual(  values=c("TRUE"="none", "FALSE"="stripe") ) +
          scale_pattern_fill_manual( values=my_palette) +
          Theme_BarHoriz +
          labs( x=NULL, y="Log2 (Proportional change)" )
hnd_Log2


# --- Plot percent change
res_scale <- "Percent"
cap_value <-  125    #5
dsPlot$display_outlier <- ifelse(dsPlot[[res_scale]] > cap_value, cap_value, dsPlot[[res_scale]]) # Cap values
yMax <- cap_value  # abs( ceiling( max( dsPlot[[res_scale]]) ) )
yMin <- floor(   min( dsPlot[[res_scale]]) )
hnd_Percent <- ggplot( dsPlot, aes(x=.data[[aggr_at_rank]], y=display_outlier,
                                   fill=.data[[display_rank]], pattern=Q.value <sign_thres) ) +
              geom_bar_pattern(stat="identity",  position="dodge", alpha=0.7, linewidth=0.35, 
                               color="#777777", pattern_color="#DDDDDD", pattern_fill="#EEEEEE",
                               pattern_angle=45, pattern_density=0.2,
                               pattern_spacing=0.04, pattern_key_scale_factor=0.6) + 
              geom_hline( aes(yintercept=0), size=1, linetype="solid", color="#544444") +
              # Display outliers as text values
              geom_text( aes(label=ifelse(.data[[res_scale]] > cap_value, round(.data[[res_scale]]), "")), hjust=1, size=3, color="#333333") +
              coord_flip() +
              facet_wrap( ~FixEff_Value, scales="free_x", nrow=1 ) +
              scale_fill_manual(      values = my_palette) +
              scale_pattern_manual( values=c("TRUE"="none", "FALSE"="stripe") ) +
              scale_pattern_fill_manual( values=my_palette) +
              ylim( yMin, yMax ) +
              Theme_BarHoriz +
              labs( x=NULL, y="Percent change" ) 
hnd_Percent



# --- Plot Significance_vs_Time ----------------------------------------------------------------
# Subset and create summary table
df_Long <- res_Maaslin[ c(display_rank,"FixEff_Value","Q.value") ]       # Q.value  N.not.0
Sig_to_Time <- df_Long %>%
                group_by(FixEff_Value, .data[[display_rank]]) %>%
                summarize( Sig_Count=sum(Q.value < sign_thres, na.rm=TRUE), .groups="drop") 

# Create the class "Total" which store the counts for all "Genera" that had a 
# significantly changed over the FT cycles
total <- Sig_to_Time %>% group_by( FixEff_Value ) %>% summarize( Sig_Count=sum(Sig_Count), .groups = "drop") 
total$Class <- "Total"
Sig_to_Time <- rbind( Sig_to_Time, total )
Sig_to_Time$Time <- sapply( Sig_to_Time$FixEff_Value, function(x) as.numeric(strsplit(x,"C")[[1]][2]) )
Sig_to_Time <- relocate( Sig_to_Time, "Time", .after=1 )

# Filter out Classes where there was no significance at any FT cycle (Sig_Count)
Sig_to_Time <- Sig_to_Time %>% group_by( .data[[display_rank]] ) %>%
                filter( any(Sig_Count > 0) )  %>%  ungroup()

# Convert counts into proportions (%). We will count how many total genera 
# are found at each FT cycle, so that the proportions are "normalized" to the 
# FT cycle
Sig_to_Time$Sig_Prop <- 0
list_TP <- unique(Sig_to_Time$FixEff_Value)
taxonCount.FT_cycles <- data.frame()
for( ii in seq(1,length(list_TP))){
  # Find specific FT cycle and calculate total number of non-zero taxons
  row_mask <- grep( list_TP[ii], rownames(reduc_ASVs) )
  sub_reduc_ASVs <- reduc_ASVs[ row_mask, ]
  sub_reduc_ASVs <- sub_reduc_ASVs[ , colSums(sub_reduc_ASVs) != 0 ]
  # for the ii-th FT cycle, transform counts into proportions
  row_mask  <- Sig_to_Time$FixEff_Value %in% list_TP[ii]
  Sig_to_Time[ row_mask, "Sig_Prop" ] <- Sig_to_Time[ row_mask, "Sig_Count" ] / ncol(sub_reduc_ASVs)
  taxonCount.FT_cycles <- rbind( taxonCount.FT_cycles, 
                                 data.frame( Time=list_TP[ii], Tot_lower=ncol(sub_reduc_ASVs) ) )
}


# +++ Add color information +++++++++
# Import color LUT for microbiome ranks and set order to match color palette
# +++ Add color information +++++++++
ref_ranks <- read.table( paste0( path_2_data, "Data_Seq/RefList_",display_rank,".txt"), header=T, sep="\t")
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(display_rank, "Colour")
# Merge taxonomic and color information to dsPlot
Sig_to_Time <- left_join( Sig_to_Time, ref_ranks, by=display_rank )
# Factorize the aggr_at_rank (e.g., Genus); then subset to take the "top" taxa
Sig_to_Time[[display_rank]] <- factor( Sig_to_Time[[display_rank]], levels=c(ref_ranks$Class, "Total") )
my_palette <- setNames( c(ref_ranks$Colour, "#333333"), c(ref_ranks$Class, "Total") )

# Define the scale factor between the left y-axis (Sig_Count) and right y-axis 
# (Percentage), so that we can plot matched scales. Thus the line are overlap
# exactly and we can show two y-axis scales.
max_count    <- max(Sig_to_Time$Sig_Count, na.rm = TRUE)
max_percent  <- max(Sig_to_Time$Sig_Prop , na.rm = TRUE)
scale_factor <- max_percent / max_count

# NOTE: The Erysipelotrichia class is always at Sig_Count = 1 and is overshadow 
# by the  other plot line. Here we add a small shift to make it popo ou in the chart.
Sig_to_Time$Sig_Count[Sig_to_Time$Class=="Erysipelotrichia"] <- 1.15

# Plot the data
hnd_2Ys <- ggplot(Sig_to_Time, aes(x=Time) ) +
              # Primary axis: counts
              geom_line( aes( y=Sig_Count, color=.data[[display_rank]], group=.data[[display_rank]]), size=1) +
              geom_point(aes( y=Sig_Count, color=.data[[display_rank]], group=.data[[display_rank]]), size=2, shape=21) +
              # Secondary axis: percentage
              geom_line( aes( y=Sig_Prop /scale_factor, color=.data[[display_rank]], group=.data[[display_rank]]), 
                         size=0.0, linetype="dashed") +
              
              scale_color_manual(values=my_palette, name=NULL) +
              scale_fill_manual( values=my_palette, name=NULL) +

              scale_x_continuous(breaks=unique(Sig_to_Time$Time), labels=unique(Sig_to_Time$FixEff_Value) ) +
              scale_y_continuous( name="Genera sign. changed (n)",
                                  breaks=seq(0, 150, 2), 
                                  sec.axis=sec_axis(~ . * scale_factor*100, name="Genera sign. changed (%)")  # Secondary y-axis
                                ) +
              ## Add the labels first
              # geom_text_repel( data = Sig_to_Time %>% filter( Time == 6 ),
              #                  aes(y = Sig_Count, label = Class, color = .data[[display_rank]]),  # Use y=Sig_Count for label positioning
              #                  max.overlaps=20, size=3,
              #                  direction="y", xlim=c(8, 12), hjust=1,  box.padding=.6,
              #                  segment.size=.4, segment.alpha=.5, segment.linetype="dotted",
              #                  segment.ncp=1,   segment.angle=90, segment.curvature=0.9
              #                  ) +
              
              # Ensure clipping is off
              coord_cartesian(clip = "off") +
  
              xlab("FT cycles") +  # Axis labels and theme
              main_layout + theme(
                panel.border = element_rect(size=0.8, colour = "#666666", linetype = "solid", fill=NA) ,
                panel.grid.major.x=element_line(size=0.25, colour="#cbcbcb", linetype="dashed"),
                panel.grid.major.y=element_line(size=0.25, colour="#cbcbcb", linetype="dashed"),
                # legend.position = "left",
                # plot.margin = margin(1, 110, 1, 1) ,
                axis.title.y.left  = element_text(color = "#005E76",size=figFontSize-2, vjust=2,   hjust=.5),
                axis.title.y.right = element_text(color = "#A23232",size=figFontSize-2, vjust=1.5, hjust=.5)
              )
hnd_2Ys


# --- Save plots ---------------------------------------------------------------
if (ref_TP == "C[0123456]"){    # using C0 as reference point
  
ggdraw() +
  draw_plot(hnd_2Ys, x=.04, y=.01, width=.38, height =.95) +
  draw_plot_label(label = c("A"), size=17, x = c(.02), y = c(.99) )
ggsave( paste0(path_saveFig, "Figure_3_A__reference_C0.tiff"), width=10, height=3.2, dpi=300)

ggdraw() +
  draw_plot(hnd_Log2, x=.05, y=.01, width=.90, height =.95) +
  draw_plot_label(label = c("C"), size=17, x = c(.05), y = c(.99) )
ggsave( paste0(path_saveFig, "Figure_3_C__reference_C0.tiff"), width=10, height=4.5, dpi=300) 

ggdraw() +
  draw_plot(hnd_Percent, x=.05, y=.01, width=.90, height =.95) +
  draw_plot_label(label = c("C"), size=17, x = c(.05), y = c(.99) )
ggsave( paste0(path_saveFig, "Figure_3_P__reference_C0.tiff"), width=10, height=4.5, dpi=300) 


}else {   # "C[123456]"         # using C1 as reference point
  
  ggdraw() +
    draw_plot(hnd_2Ys, x=.54, y=.01, width=.30, height =.95) +
    draw_plot_label(label = c("B"), size=17, x = c(.52), y = c(.99) )
  ggsave( paste0(path_saveFig, "Figure_3_B__reference_C1.tiff"), width=10, height=3.2, dpi=300)

  ggdraw() +
    draw_plot(hnd_Log2, x=.09, y=.01, width=.48, height =.95) +
    draw_plot_label(label = c("D"), size=17, x = c(.05), y = c(.99) )
  ggsave( paste0(path_saveFig, "Figure_3_D__reference_C1.tiff"), width=10, height=2, dpi=300) 
  
  ggdraw() +
    draw_plot(hnd_Percent, x=.09, y=.01, width=.48, height =.95) +
    draw_plot_label(label = c("D"), size=17, x = c(.05), y = c(.99) )
  ggsave( paste0(path_saveFig, "Figure_3_P__reference_C1.tiff"), width=10, height=2, dpi=300) 
  
}




