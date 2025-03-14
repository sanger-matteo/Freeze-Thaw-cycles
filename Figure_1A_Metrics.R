# +++++ DESCRIPTION - Diversity Metrics ++++++++++++++++++++++++++++++++++++++++
#
# The script calculate several Alpha and Beta diversity metric of the microbiome 
# composition to analyse the differences between freeze-thaw cycle. 
#
# ANOVA test, Kruskal Wallis test, Wilcoxon and T-test are performed to assess
# if there are significant differences between F-T cycle groups. We also create
# a Linear Mix Model and perform test.
# (NOTE: These are run in the script TestStatistics_Metric)
#
# All statistical tests are stored in a summary table (saved as txt). 
#
# The sub-figures created are:
# ---> for Figure 1:
# - F 1 A - DNA concentration      --> hnd_DNA
# - F 1 B - Shannon index          --> hnd_Shan
# - F 1 C - BC similarity          --> hnd_BC
# - F 1 D - Aitchison Similarity   --> hnd_Aitch
#
# ---> for Supplementary Figure 2:
# - SF 2 A - Richness              --> hnd_Rich
# - SF 2 B - Jaccard Similarity    --> hnd_Jacc
# - SF 2 C - Chao Similarity       --> hnd_Chao
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# --- Load Libraries, functions and data tables
library(ggplot2)
library(cowplot)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(rstatix)
library(vegan)

# Round to the nearest upper or lower "round" with a specific precision factor
round_up_to <- function(num, round_fact=0.01){    ceiling( num / round_fact ) *round_fact  }
round_dw_to <- function(num, round_fact=0.01){    floor(   num / round_fact ) *round_fact  }

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

# --- Aggregate ASVs at rank X:
aggr_at_rank <- "Genus"
if (aggr_at_rank != "ASV" ){
  returned_list <- aggregate_at_Rank( ASV, TAXO, aggr_at_rank, F, F, F )
  reduc_ASVs <- returned_list[[1]]
  reduc_Taxo <- returned_list[[2]]
}else{
  reduc_ASVs <- ASV
  reduc_Taxo <- TAXO
}

# --- Filter features (based on RA):
relative_features <- T
if(relative_features==TRUE){
  # - thres_abundance:     min relative abundance of a feature
  # - thres_freqCohort:    min frequency of a feature in the cohort
  # Generally, use Rel.Abu. less than 0.03%, and freq. greater than 10% of the cohort
  # Also, adjust the TAXO table to contain only features present in the cohort.
  thres_abundance  <- 0.003 /100
  thres_freqCohort <- 0.01  /100
  taxa_RA_Freq <- data.frame( Taxa   = colnames(reduc_ASVs[,-1]) ,
                              N_Obs_above_RA_thres = colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) ,
                              Percent_of_cohort    = colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] *100,
                              Keep_taxa            = ( colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] ) >= thres_freqCohort )
  rownames(taxa_RA_Freq) <- seq(1,nrow(taxa_RA_Freq))
  mask_filterASV <- as.vector( ( colSums((reduc_ASVs[,-1]/rowSums(reduc_ASVs[,-1])) > thres_abundance, na.rm=T) / dim(reduc_ASVs[,-1])[1] ) >= thres_freqCohort )
  sprintf(" --> # taxon to be remove = %d", sum(!mask_filterASV) )
  reduc_ASVs <- reduc_ASVs[ , c(T,mask_filterASV) ]
  # Mask and reduce Taxo to the essential taxonomy, that used in reduc_ASVs table
  reduc_Taxo <- reduc_Taxo[ mask_filterASV , ]
}

# Make ASV table fully numeric and with only ASV features
rownames(reduc_ASVs) <- reduc_ASVs$Sample_ID
reduc_ASVs <- reduc_ASVs[,-1]

# Empty variable to store the summary of test statistics (e.g., p-value and F-ratio)
summary_stats <- data.frame()




# ******************************************************************************
# ************************* ALPHA DIVERSITY METRICS ****************************
# ******************************************************************************


# ---- DNA Concentrations ------------------------------------------------------
dfPlot <- MetaData %>% dplyr::select( ., DNA_Extract_nguL, Exp_Name, Exp_FTcycle)
colnames(dfPlot) <- c( "Values","RowID","FT_cycle" )

hnd_DNA <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=Exp_FTcycle)) + 
            geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
            geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                          lwd=0.2, fatten=2.5 ) +  
            scale_y_continuous( expand = c(0, NA), limits = c(100, 600) ) +
            scale_fill_manual( values = c("#00C385","#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF")  ) + 
                                      # c("#FF7052","#A12844","#C93255","#DF4B72","#F06991","#FB8CB1","#FFB8D1")
            xlab("")  +  ylab("DNA (ng/µL)") +
            Theme_BoxPlot 
if (save_subFigs){
  hnd_DNA
  ggsave( paste0(path_saveFig, "DNA_Conc_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]

# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "DNA_Conc" )
summary_stats <- rbind ( summary_stats , return_list )



# ---- Shannon diversity -------------------------------------------------------
# Calculate Shannon diversity index
M_Shan <- data.frame( "Values"=diversity(reduc_ASVs, index="shannon") )
M_Shan$RowID    <- unlist( sapply( rownames(M_Shan), function(x) strsplit(x,"_")[[1]][2], simplify="array" ) )
M_Shan$FT_cycle <- unlist( sapply( rownames(M_Shan), function(x) strsplit(x,"_")[[1]][1] ) )
dfPlot <- M_Shan
# Transform to wider format, to better check results in a table format
M_Shan <- pivot_wider( M_Shan, names_from=FT_cycle, values_from=Values )
# Automatically find y-axis limits 
y_Max <- round_up_to( max(dfPlot$Values), 2 )
y_Min <- 0      # round_dw_to( min(dfPlot$Values), 3 )

hnd_Shan <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) +  
          scale_y_continuous( expand = c(0, NA), limits=c(y_Min,y_Max) ) +
          scale_fill_manual( values = c("#00C385","#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Shannon index") +
          Theme_BoxPlot 
if (save_subFigs){
  hnd_Shan 
  ggsave( paste0(path_saveFig, "Metrics_Shannon_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Shannon" )
summary_stats <- rbind ( summary_stats , return_list )



# ---- Simpson diversity -------------------------------------------------------
# Calculate Simpson diversity index
M_Simp <- data.frame( "Values" = diversity(reduc_ASVs, index = "simpson") )
M_Simp$RowID    <- unlist( sapply( rownames(M_Simp), function(x) strsplit(x,"_")[[1]][2], simplify="array" ) )
M_Simp$FT_cycle <- unlist( sapply( rownames(M_Simp), function(x) strsplit(x,"_")[[1]][1] ) )
dfPlot <- M_Simp
# Transform to wider format, to better check results in a table format
M_Simp <- pivot_wider( M_Simp, names_from=FT_cycle, values_from=Values )
# Automatically find y-axis limits 
y_Max <- round_up_to( max(dfPlot$Values), 0.1 )
y_Min <- round_dw_to( min(dfPlot$Values), 0.3 )

hnd_Simp <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) + 
          scale_y_continuous( expand = c(0, NA), limits=c(y_Min,y_Max) ) +
          scale_fill_manual( values = c("#00C385","#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Simpson index") +
          Theme_BoxPlot 
if (save_subFigs){ 
  hnd_Simp
  ggsave( paste0(path_saveFig, "Metrics_Simpson_BoxPlot.tiff"), width=3, height=3, dpi=300) 
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Simpson" )
summary_stats <- rbind ( summary_stats , return_list )



# ---- Richness diversity -------------------------------------------------------
# Calculate Simpson diversity index
M_Rich <- data.frame( "Values" = specnumber(reduc_ASVs) )
M_Rich$RowID    <- unlist( sapply( rownames(M_Rich), function(x) strsplit(x,"_")[[1]][2], simplify="array" ) )
M_Rich$FT_cycle <- unlist( sapply( rownames(M_Rich), function(x) strsplit(x,"_")[[1]][1] ) )
dfPlot <- M_Rich
# Transform to wider format, to better check results in a table format
M_Rich <- pivot_wider( M_Rich, names_from=FT_cycle, values_from=Values )
# Automatically find y-axis limits 
y_Max <- round_up_to( max(dfPlot$Values), 50 )
y_Min <- 0       # round_dw_to( min(dfPlot$Values), 0.01 )

hnd_Rich <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) + 
          scale_y_continuous( expand = c(0, NA), limits=c(y_Min,y_Max) ) +
          scale_fill_manual( values = c("#00C385","#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Richness") +
          Theme_BoxPlot 
if (save_subFigs){
  hnd_Rich
  ggsave( paste0(path_saveFig, "Metrics_Richness_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Richness" )
summary_stats <- rbind ( summary_stats , return_list )




# ******************************************************************************
# ************************* BETA DIVERSITY METRICS *****************************
# ******************************************************************************


# ---- Jaccard diversity -------------------------------------------------------
# Calculate the metric. However, the vegan calculates all possible combinations
# of samples pairs. We are interested in all of them: namely only C0, as a 
# reference, against all other CX, for each individual.
# We will gather only the desired combinations for plotting and test statistics.
unique_FT_ID <- paste0( MetaData$Exp_FTcycle, "_", MetaData$Exp_Name ) 
MetricMatrix <- as.data.frame( as.matrix( vegan::vegdist( reduc_ASVs , method="jaccard" ) ))
idx <- data.frame( iRow=rep(seq(1,5),6), iCol=seq(6,35) )
# Extract the desired calculated metrics from the matrix
ds <- list()
for (ii in seq(1,dim(idx)[1])){
  ds <- c(ds, MetricMatrix[ idx[ii,1], idx[ii,2] ] )
}
select_paired_ID <- paste0( unique_FT_ID [idx[,1]],"-",unique_FT_ID [idx[,2]] )
M_Jacc           <- data.frame( "Values" = unlist(ds) ) 
M_Jacc$FT_cycle  <- paste0( MetaData$Exp_FTcycle[idx[,1]],"-",MetaData$Exp_FTcycle[idx[,2]] )
M_Jacc$RowID     <- MetaData$Exp_Name[idx[,2]]
dfPlot <- M_Jacc
M_Jacc <- pivot_wider( M_Jacc, names_from=FT_cycle, values_from=Values )

hnd_Jacc <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) + 
          scale_y_continuous( expand = c(0, NA), limits=c(0, 1)) +
          scale_fill_manual( values = c("#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Jaccard similarity") +
          Theme_BoxPlot
if (save_subFigs){
  hnd_Jacc
  ggsave( paste0(path_saveFig, "Metrics_Jaccard_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Jaccard" )
summary_stats <- rbind ( summary_stats , return_list )



# ---- Bray-Curtis diversity -------------------------------------------------------
# Calculate the metric. However, the vegan calculates all possible combinations
# of samples pairs. We are interested in all of them: namely only C0, as a 
# reference, against all other CX, for each individual.
# We will gather only the desired combinations for plotting and test statistics.
unique_FT_ID <- paste0( MetaData$Exp_FTcycle, "_", MetaData$Exp_Name ) 
MetricMatrix <- as.data.frame( as.matrix( vegan::vegdist( reduc_ASVs , method="bray" ) ))
idx <- data.frame( iRow=rep(seq(1,5),6), iCol=seq(6,35) )
# Extract the desired calculated metrics from the matrix
ds <- list()
for (ii in seq(1,dim(idx)[1])){
  ds <- c(ds, MetricMatrix[ idx[ii,1], idx[ii,2] ] )
}
select_paired_ID <- paste0( unique_FT_ID [idx[,1]],"-",unique_FT_ID [idx[,2]] )
M_Bray           <- data.frame( "Values" = unlist(ds) ) 
M_Bray$FT_cycle  <- paste0( MetaData$Exp_FTcycle[idx[,1]],"-",MetaData$Exp_FTcycle[idx[,2]] )
M_Bray$RowID     <- MetaData$Exp_Name[idx[,2]]
dfPlot <- M_Bray
M_Bray <- pivot_wider( M_Bray, names_from=FT_cycle, values_from=Values )

hnd_BC <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) + 
          scale_y_continuous( expand = c(0, NA), limits=c(0, 1)) +
          scale_fill_manual( values = c("#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Bray-Curtis similarity") +
          Theme_BoxPlot
if (save_subFigs){
  hnd_BC
  ggsave( paste0(path_saveFig, "Metrics_BrayCurtis_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Bray-Curtis" )
summary_stats <- rbind ( summary_stats , return_list )



# ---- Aitchison diversity -------------------------------------------------------
# Calculate the metric. However, the vegan calculates all possible combinations
# of samples pairs. We are interested in all of them: namely only C0, as a 
# reference, against all other CX, for each individual.
# We will gather only the desired combinations for plotting and test statistics.
unique_FT_ID <- paste0( MetaData$Exp_FTcycle, "_", MetaData$Exp_Name ) 
MetricMatrix <- as.data.frame( as.matrix( vegan::vegdist( reduc_ASVs , method="robust.aitchison" ) ))
idx <- data.frame( iRow=rep(seq(1,5),6), iCol=seq(6,35) )
# Extract the desired calculated metrics from the matrix
ds <- list()
for (ii in seq(1,dim(idx)[1])){
  ds <- c(ds, MetricMatrix[ idx[ii,1], idx[ii,2] ] )
}
select_paired_ID <- paste0( unique_FT_ID [idx[,1]],"-",unique_FT_ID [idx[,2]] )
M_Aich           <- data.frame( "Values" = unlist(ds) ) 
M_Aich$FT_cycle  <- paste0( MetaData$Exp_FTcycle[idx[,1]],"-",MetaData$Exp_FTcycle[idx[,2]] )
M_Aich$RowID     <- MetaData$Exp_Name[idx[,2]]
dfPlot <- M_Aich
M_Aich <- pivot_wider( M_Aich, names_from=FT_cycle, values_from=Values )
# Automatically find y-axis limits 
y_Max <- round_up_to( max(dfPlot$Values), 10 )
y_Min <- 0       # round_dw_to( min(dfPlot$Values), 0.01 )

hnd_Aitch <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) + 
          scale_y_continuous( expand = c(0, NA) , limits=c(y_Min,y_Max) ) +
          scale_fill_manual( values = c("#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Aitchison similarity") +
          Theme_BoxPlot
if (save_subFigs){
  hnd_Aitch
  ggsave( paste0(path_saveFig, "Metrics_Aitchison_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Aitchison" )
summary_stats <- rbind ( summary_stats , return_list )



# ---- CHAO diversity -------------------------------------------------------
# Calculate the metric. However, the vegan calculates all possible combinations
# of samples pairs. We are interested in all of them: namely only C0, as a 
# reference, against all other CX, for each individual.
# We will gather only the desired combinations for plotting and test statistics.
unique_FT_ID <- paste0( MetaData$Exp_FTcycle, "_", MetaData$Exp_Name ) 
MetricMatrix <- as.data.frame( as.matrix( vegan::vegdist( reduc_ASVs , method="chao" ) ))
idx <- data.frame( iRow=rep(seq(1,5),6), iCol=seq(6,35) )
# Extract the desired calculated metrics from the matrix
ds <- list()
for (ii in seq(1,dim(idx)[1])){
  ds <- c(ds, MetricMatrix[ idx[ii,1], idx[ii,2] ] )
}
select_paired_ID <- paste0( unique_FT_ID [idx[,1]],"-",unique_FT_ID [idx[,2]] )
M_Chao           <- data.frame( "Values" = unlist(ds) ) 
M_Chao$FT_cycle  <- paste0( MetaData$Exp_FTcycle[idx[,1]],"-",MetaData$Exp_FTcycle[idx[,2]] )
M_Chao$RowID     <- MetaData$Exp_Name[idx[,2]]
dfPlot <- M_Chao
M_Chao <- pivot_wider( M_Chao, names_from=FT_cycle, values_from=Values )
# Automatically find y-axis limits 
y_Max <- round_up_to( max(dfPlot$Values), 0.01 )
y_Min <- 0       # round_dw_to( min(dfPlot$Values), 0.01 )

hnd_Chao <- ggplot( dfPlot, aes(x=FT_cycle, y=Values, fill=FT_cycle)) + 
          geom_jitter(  aes(fill=FT_cycle), color="#666666", alpha=0.5, size=1.2, width=0.2) +  
          geom_boxplot( aes(fill=FT_cycle), color="#333333", alpha=0.8, outlier.shape=NA, 
                        lwd=0.2, fatten=2.5 ) + 
          scale_y_continuous( expand = c(0, NA), limits=c(y_Min,y_Max) ) +
          scale_fill_manual( values = c("#00C385","#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF") ) +
          xlab("")  +  ylab("Chao similarity") +
          Theme_BoxPlot
if (save_subFigs){
  hnd_Chao
  ggsave( paste0(path_saveFig, "Metrics_Chao_BoxPlot.tiff"), width=3, height=3, dpi=300)
}
# Perform LMM test
return_list <- Test_LMM_Metrics( dfPlot ) 
LMM_sig  <- return_list[[1]]
LMM_pair <- return_list[[2]]
# Perform a range of statistical tests (param. and non param.) and store results in summary table
return_list <- Test_Significance_Metrics( dfPlot, "Values", "FT_cycle", "RowID", "Chao" )
summary_stats <- rbind ( summary_stats , return_list )




# --- Save Figure and Table -----------------------------------------------------------------------------
# Save summary table of statistical tests
write_tsv( summary_stats, paste0( path_2_data, "Res_Tables/Table_ANOVA_DivMetrics_",aggr_at_rank,".txt") )

# Place sub-plots together into one figure
ggdraw() +
  draw_plot(hnd_DNA ,  x = .00, y = .0, width = .24, height = .91) +
  draw_plot(hnd_Shan,  x = .25, y = .0, width = .24, height = .91) +
  draw_plot(hnd_BC,    x = .50, y = .0, width = .24, height = .91) +
  draw_plot(hnd_Aitch, x = .75, y = .0, width = .24, height = .91) +
  draw_plot_label(label = c("A","B","C","D"), size = 17, 
                  x = c(.00,.25,.50,.75), y = c(.99,.99,.99,.99) )
path_saveFig <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/Figures/" )
ggsave( paste0(path_saveFig, "Figure_1AD_",aggr_at_rank,".tiff"), width=10, height=3, dpi=300)


ggdraw() +
  draw_plot(hnd_Rich,  x = .12, y = .0, width = .24, height = .9) +
  draw_plot(hnd_Jacc,  x = .37, y = .0, width = .24, height = .9) +
  draw_plot(hnd_Chao,  x = .62, y = .0, width = .24, height = .9) +
  draw_plot_label(label = c("A","B","C"), size = 17, 
                  x = c(.12,.37,.62), y = c(.99,.99,.99) )
path_saveFig <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/Figures/" )
ggsave( paste0(path_saveFig, "Suppl_Figure_2_AB_",aggr_at_rank,".tiff"), width=10, height=3, dpi=300)

 












