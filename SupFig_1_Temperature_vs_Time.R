# +++++ DESCRIPTION --- Temperature_vs_Time.R ++++++++++++++++++++++++++++++++++
#
# Analyse and plot the results of the thawing experiments, which is where we 
# followed the temperature changes in 3 different types of tubes containing 
# different amounts of fecal matter. The temperature was generally taken every 
# 5 min during the entire thawing process.
# The temperature limit for the IR thermometer was -35°C.
#
# The script will plot the trend line of the temperature with a solid color, 
# which is calculated from the average of repeated experiments; the script also 
# adds a ribbon area that represent the standard deviation.
# Lastly, on top of the plots we chart blue boxes that indicate at which 
# temperature condition we placed the samples.
#
# In the INPUT excel file we have:
#  - sheet 1, is the slow thawing experiments
#  - sheet 2, is the fast thawing experiments
# in both, each row is an separate set of tube and experiment
#
# This script generate the Supplementary Figure 1 A and 1 B.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Load Libraries and functions
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2) 
library(readxl)


cDir_Rscript <- setwd("/Users/sm/Documents/R/FT_Cycle")
path_2_data  <- paste0( cDir_Rscript, "/Freeze-Thaw_Results/" )
# Import custom functions, specifically contains those to calculated microbiome diversity metrics
source( paste0(path_2_data, "Plotting_Repository.R") )


# ---- Sheet 1 - Slow Thaw ---------------------------------------------------------------------------------------
# Load and plot the temperature measurements for the SLOW thawing process
dsTemp <- read_excel( paste0(path_2_data, "Data_Seq/Temperature_Record_FT_fecal.xlsx"), sheet=1 )

dsTemp <- as.data.frame(dsTemp)
dsTemp <- dsTemp[-nrow(dsTemp), ]
dsTemp <- dsTemp[dsTemp$Date != "24.04.24", ]
dsTemp <- dsTemp[!sapply( dsTemp,function(x) all(is.na(x)) )]

# Remove becker sample with very little mass
dsTemp <- dsTemp[ !(dsTemp$BioRep=="FT0" & dsTemp$TechRep=="R3" & dsTemp$TubeType=="15mL"), ]

# Format all measured values as numeric and transform to long format data
dsTemp[-c(1,2,3,4)] <- sapply( dsTemp[-c(1,2,3,4)], function(x) as.numeric(x))
dsTemp <- dsTemp[ , -5]
# Calculate mean and standard deviation for each group (Tube Type), then 
# transfrom to long format and join them togeteher as a data frame for plotting
dsMean <- dsTemp %>% group_by(TubeType) %>%
  summarise(across(where(is.numeric), \(x) mean(x,na.rm = TRUE)))
dsStdev <- dsTemp %>% group_by(TubeType) %>%
  summarise(across(where(is.numeric), \(x) sd(x,na.rm = TRUE)))

LongMean <- dsMean %>% pivot_longer( cols=-TubeType,
                                     names_to="Time",
                                     values_to="MeanTemp") 
LongStdev <- dsStdev %>% pivot_longer( cols=-TubeType,
                                       names_to="Time",
                                       values_to="StdevTemp") 

LongTemp <- left_join(LongMean,LongStdev)
LongTemp$Time <- as.numeric(LongTemp$Time)
# Add information about the weight of faeces in the tubes, to display in the legend
LongTemp$TubeType <- ifelse( LongTemp$TubeType == "2mL", "2mL: 1.5g faeces", LongTemp$TubeType )
LongTemp$TubeType <- ifelse( LongTemp$TubeType == "5mL", "5mL: 4g faeces", LongTemp$TubeType )
LongTemp$TubeType <- ifelse( LongTemp$TubeType == "15mL", "15mL: 8g faeces", LongTemp$TubeType )
LongTemp$TubeType <- factor( LongTemp$TubeType, levels= c("2mL: 1.5g faeces","5mL: 4g faeces","15mL: 8g faeces"))


# ---- PLOT Slow Thawing Approach
# Plot temperature change over time: line=Mean + shadow_area=StDev
hnd_Slow <- ggplot( LongTemp, aes(x=Time, y=Temp, fill=TubeType, color=TubeType)) + 
          # Reference gridlines at temperature changes and 0 degrees
          # NOTE: plot it first so that it does not overlay and remains in background
          geom_hline( yintercept=seq(-40,30,10), size=0.5, colour="#DDDDDD", linetype="dashed") +
          geom_vline( xintercept=c(0,10,20,25,30,35,40,45,50,55,60,65,70), size=0.5,  colour="#DDDDDD", linetype="dashed") +
          geom_vline( xintercept=c(0,20,40), size=0.7, colour="#BBBBBB") +
          geom_hline( yintercept=0, size=0.7, color="#888888", ) +
          # Plot main data, temperature change over time: line=Mean + shadow_area=StDev
          geom_line( aes(y=MeanTemp), size=1)+
          geom_ribbon(aes(y=MeanTemp, ymin=MeanTemp-StdevTemp, ymax=MeanTemp+StdevTemp), alpha=0.2, linetype="blank")+
          geom_point(aes(y=MeanTemp), size=1.4, shape=21)+
          # Format plot with nice theme and add information of the temperature changes
          labs( x="Time (min)", y="Temperature (˚C)" ) +
          scale_color_manual(values=c("#AA1A36","#D25E3A","#E9AD59"), name=NULL) +     # add colors and remove legend title
          scale_fill_manual( values=c("#AA1A36","#D25E3A","#E9AD59"), name=NULL) +
          scale_x_continuous( breaks=unique(LongTemp$Time) ) +
          Theme_Temp2Time +
          # Annotate plot with additional information: text number to indicate the 
          # temp. change, arrow and rectangle to show the time period at set temp.
          annotate("rect", xmin=1,  xmax=19, ymin=28, ymax=33, fill="#20ABC7", alpha=0.6, color="#666666", size=0.5) +
          annotate("rect", xmin=21, xmax=39, ymin=28, ymax=33, fill="#20ABC7", alpha=0.3, color="#666666", size=0.5) +
          annotate("rect", xmin=41, xmax=69, ymin=28, ymax=33, fill="#20ABC7", alpha=0.1, color="#666666", size=0.5) +
          annotate("text", x=c(6,26,46), y=c(28,28,28), label=c("–25˚C","+4˚C","RT"),vjust=-0.5, hjust=0.5 ) +
          # Manually position legend inside chart area
          theme(  legend.position=c(0.95, 0.35), legend.justification=c(1, 1) )
hnd_Slow




# ---- Sheet 2 - Fast Thaw --------------------------------------------------------------------------------------
# Load and plot the temperature measurements for the FAST thawing process
dsTemp <- read_excel( paste0(path_2_data, "Data_Seq/Temperature_Record_FT_fecal.xlsx"), sheet=2 )

dsTemp <- as.data.frame(dsTemp)
dsTemp <- dsTemp[-nrow(dsTemp), ]
dsTemp <- dsTemp[!sapply( dsTemp,function(x) all(is.na(x)) )]

# Remove becker sample with very little mass
dsTemp <- dsTemp[ !(dsTemp$BioRep=="FT2" & dsTemp$TechRep=="R2" & dsTemp$TubeType=="15mL"), ]

# Format all measured values as numeric and transform to long format data
dsTemp[-c(1,2,3,4)] <- sapply( dsTemp[-c(1,2,3,4)], function(x) as.numeric(x))
dsTemp <- dsTemp[ , -5]
# Calculate mean and standard deviation for each group (Tube Type), then 
# transfrom to long format and join them togeteher as a data frame for plotting
dsMean <- dsTemp %>% group_by(TubeType) %>%
          summarise(across(where(is.numeric), \(x) mean(x,na.rm = TRUE)))
dsStdev <- dsTemp %>% group_by(TubeType) %>%
           summarise(across(where(is.numeric), \(x) sd(x,na.rm = TRUE)))

LongMean <- dsMean %>% pivot_longer( cols=-TubeType,
                                     names_to="Time",
                                     values_to="MeanTemp") 
LongStdev <- dsStdev %>% pivot_longer( cols=-TubeType,
                                     names_to="Time",
                                     values_to="StdevTemp") 

LongTemp <- left_join(LongMean,LongStdev)
LongTemp$Time <- as.numeric(LongTemp$Time)
# Add information about the weight of faeces in the tubes, to display in the legend
LongTemp$TubeType <- ifelse( LongTemp$TubeType == "2mL", "2mL: 1.5g faeces", LongTemp$TubeType )
LongTemp$TubeType <- ifelse( LongTemp$TubeType == "5mL", "5mL: 4g faeces", LongTemp$TubeType )
LongTemp$TubeType <- ifelse( LongTemp$TubeType == "15mL", "15mL: 8g faeces", LongTemp$TubeType )
LongTemp$TubeType <- factor( LongTemp$TubeType, levels= c("2mL: 1.5g faeces","5mL: 4g faeces","15mL: 8g faeces"))

  
# ---- PLOT Fast Thawing approach
# Plot temperature change over time: line=Mean + shadow_area=StDev
hnd_Fast <- ggplot( LongTemp, aes(x=Time, group=TubeType, fill=TubeType, color=TubeType)) + #y=Temp, 
          # Reference gridlines at temperature changes and 0 degrees
          # NOTE: plot it first so that it does not overlay and remains in background
          geom_hline( yintercept=seq(-40,30,10), size=0.5, colour="#DDDDDD", linetype="dashed") +
          geom_vline( xintercept=c(0,5,10,15,20,25,30,35,40,45,50,55,60), size=0.5,  colour="#DDDDDD", linetype="dashed") +
          geom_vline( xintercept=0, size=0.7, colour="#BBBBBB") +
          geom_hline( yintercept=0, size=0.7, colour="#888888", ) +
          # Plot main data, temperature change over time: line=Mean + shadow_area=StDev
          geom_line( aes(y=MeanTemp), size=1)+
          geom_ribbon(aes(y=MeanTemp, ymin=MeanTemp-StdevTemp, ymax=MeanTemp+StdevTemp), alpha=0.2, linetype="blank")+
          geom_point(aes(y=MeanTemp), size=1.4, shape=21)+
          # Format plot with nice theme and add information of the temperature changes
          labs( x="Time (min)", y="Temperature (˚C)" ) +
          scale_color_manual(values=c("#AA1A36","#D25E3A","#E9AD59"), name=NULL) +     # add colors and remove legend title
          scale_fill_manual( values=c("#AA1A36","#D25E3A","#E9AD59"), name=NULL) +
          scale_x_continuous( breaks=unique(LongTemp$Time) ) +
          Theme_Temp2Time +
          # Annotate plot with additional information: text number to indicate the 
          # temp. change, arrow and rectangle to show the time period at set temp.
          annotate("rect", xmin=1,  xmax=59, ymin=28, ymax=33, fill="#20ABC7", alpha=0.1, color="#666666", size=0.5) +
          annotate("text", x=6, y=28, label="RT",vjust=-0.5, hjust=0.5 ) +
          # Manually position legend inside chart area
          theme(  legend.position=c(0.99, 0.35), legend.justification=c(1, 1) )
hnd_Fast



# --- Combine sub-figures -------------------------------------------------------------------------------------
# Place sub-plots together into one figure 
library("cowplot")

ggdraw() +
  draw_plot(hnd_Fast, x=.03, y=.0, width=.47, height=.96) +
  draw_plot(hnd_Slow, x=.53, y=.0, width=.47, height=.96) +
  draw_plot_label(label=c("A","B"), size=17, x=c(.02,.52), y=c(1,1) )

ggsave( paste0(path_2_data, "Figures/Suppl_Figure_1.tiff"), width=10, height=4, dpi=300)





 







