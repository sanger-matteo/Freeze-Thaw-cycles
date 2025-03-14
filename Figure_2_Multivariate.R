# +++++ DESCRIPTION - Multivariate Analysis ++++++++++++++++++++++++++++++++++++
#
# THis script perform several multivariate analysis of the dataset:
# - Principal component analysis (PCA), 
# - t-distributed stochastic neighbor embedding (tSNE)
# - PCoA and NMDS - for these we also run a Permanova test to see if there are 
#                   statistical difference in clustering of groups.
#
# The subfigures generated are divided amont Figure 2 and Suppl. Figure 4:
# - Fig A - PCA scores 1,2         --> hnd_PC1
# - Fig B - PCA scores 5,6         --> hnd_PC3
# - Fig C - PCA loads 5,6          --> hnd_L2
# - Fig D - NMDS - Bray            --> hnd_NMDS
#
# - SupFig A - tSNE                --> hnd_tSNE
# - SupFig B - PCA scores 3,4      --> hnd_PC2
# - SupFig C - PCA loads 1,2       --> hnd_L1
# - SupFig D - PCoA - Aitchison    --> hnd_PCoA
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
library(cowplot)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(Rtsne)
library(mixOmics)
library(ggrepel) 
library(vegan)

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
filter_features <- T
aggr_at_rank <- "Genus"
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

# Make ASV table fully numeric (place Sample_ID as rownames) and replace colnames with taxonomic names
rownames(reduc_ASVs) <- reduc_ASVs$Sample_ID
reduc_ASVs <- reduc_ASVs[,-1]
colnames(reduc_ASVs) <- reduc_Taxo[[aggr_at_rank]]    

# > > > Normalize reads to Relative abundance
reduc_ASVs <- reduc_ASVs/rowSums(reduc_ASVs)



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
# ------ tSNE ---------------------------------------------------------------------------
# Perform t-SNE and scatter plot to see clustering
perplex   <- 8                                 # approximately   ~ floor((nrow(reduc_ASVs)-1)/3)
ScaleData <- scale(reduc_ASVs)                  # Scale data
tSNE_fit  <- Rtsne( X = data.matrix( ScaleData ), perplexity = perplex,
                    pca_center = TRUE, pca_scale  = TRUE, 
                    dims = 2 , max_iter = 1000, check_duplicates = FALSE )
tSNE_ds<- tSNE_fit$Y %>% as.data.frame() 
tSNE_ds <- cbind( tSNE_ds, MetaData )        # Add metadata back to tSNE results

# X and Y are parametrized to allow to plot any combination of two PCs
X_score <- "V1" ;       X_n <- as.numeric( substr( X_score, 2, 2) )
Y_score <- "V2" ;       Y_n <- as.numeric( substr( Y_score, 2, 2) )
axLimits <- square_plot_limits( tSNE_ds[ ,X_score], tSNE_ds[ ,Y_score] )

hnd_tSNE <- plot_Scores_scatter( tSNE_ds, X_score, Y_score, "Exp_FTcycle", "Exp_Name", NULL,
                                 paste0("tSNE ",X_n), paste0("tSNE ",Y_n), axLimits )
hnd_tSNE
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_tSNE_Variance.tiff"), width=3, height=3, dpi=300)
}



# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
# --- Perform PCA -----------------------------------------------------------------------
# # Remove columns of NaN values, and select response variable (-> Obs_grouping)
# mask     <- colSums( reduc_ASVs ) == 0
# plotData <- reduc_ASVs[ , !mask]

plotData <- reduc_ASVs
grouped_at <- "Exp_FTcycle"        # choose from "Exp_Name", "Exp_FTcycle" 
Obs_grouping <- factor(MetaData[[grouped_at]])     

# Perform PCA and then create data.frames for plotting the results of PCA
FT.pca <- mixOmics::pca( plotData, ncomp = 10, center = TRUE, scale = TRUE)

scorePCA <- data.frame( PC1 = as.data.frame(FT.pca$variates)$X.PC1,
                        PC2 = as.data.frame(FT.pca$variates)$X.PC2,
                        PC3 = as.data.frame(FT.pca$variates)$X.PC3,
                        PC4 = as.data.frame(FT.pca$variates)$X.PC4,
                        PC5 = as.data.frame(FT.pca$variates)$X.PC5,
                        PC6 = as.data.frame(FT.pca$variates)$X.PC6,
                        PC7 = as.data.frame(FT.pca$variates)$X.PC7,
                        PC8 = as.data.frame(FT.pca$variates)$X.PC8,
                        Grouping = Obs_grouping,
                        Exp_Name = MetaData$Exp_Name )
loadsPCA <- data.frame( Load_1 = as.data.frame(FT.pca$load)$X.PC1,
                        Load_2 = as.data.frame(FT.pca$load)$X.PC2,
                        Load_3 = as.data.frame(FT.pca$load)$X.PC3,
                        Load_4 = as.data.frame(FT.pca$load)$X.PC4,
                        Load_5 = as.data.frame(FT.pca$load)$X.PC5,
                        Load_6 = as.data.frame(FT.pca$load)$X.PC6,
                        Load_7 = as.data.frame(FT.pca$load)$X.PC7,
                        Load_8 = as.data.frame(FT.pca$load)$X.PC8 )
row.names(loadsPCA) <- colnames(FT.pca$X)
expl_var <- FT.pca$prop_expl_var


# +++ Plot explained VARIANCE (Scree Plot) +++++++++++++++++++++++++++++++++++++++++++++++++++++++
varPlot <- as.data.frame(expl_var)
varPlot$PCs <- factor( rownames(varPlot), levels=rownames(varPlot) )
varPlot$CumSum <- cumsum(varPlot$X) *100
varPlot$X <- varPlot$X *100
hnd <- ggplot( varPlot ,aes(x=PCs, y=X, group=2) ) +
        geom_col(  size=0.2, color="#999999", fill="#299CD2") +
        geom_point( aes(x=PCs, y=CumSum/3), color="#b82744") +
        geom_line(  aes(x=PCs, y=CumSum/3), color="#b82744") +
        scale_y_continuous( sec.axis = sec_axis(~.*3, name = expression(paste("Cum. ", sigma^2, " explained [%]",sep = "")))) +
        ggtitle( "" ) + xlab( "" ) + 
        ylab( expression(paste(sigma^2, " explained [%]", sep = "")) ) +
        theme( axis.text.x = element_text(angle = 45, hjust=1)) +
        Theme_PCA
hnd
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCA_Variance.tiff"), width=3, height=3, dpi=300)
}


# +++ Plot PRINCIPAL COMPONENTS score +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X and Y are parametrized to allow to plot any combination of two PCs
X_score <- "PC1" ;       X_n <- as.numeric( substr( X_score, 3, 3) )
Y_score <- "PC2" ;       Y_n <- as.numeric( substr( Y_score, 3, 3) )
axLimits <- square_plot_limits( scorePCA[ ,X_score], scorePCA[ ,Y_score] )
hnd_PC1 <- plot_Scores_scatter( scorePCA, X_score, Y_score, "Grouping", "Exp_Name", NULL, axisLimits=axLimits,
                                paste0( "PC ",X_n," (", round(expl_var[[1]][[X_n]]*100, 1) ,"%)"), 
                                paste0( "PC ",Y_n," (", round(expl_var[[1]][[Y_n]]*100, 1) ,"%)") )
hnd_PC1
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCAScores_",grouped_at,"_",X_n,"-",Y_n,".tiff"), width=4, height=4, dpi=300)
}


X_score <- "PC3" ;       X_n <- as.numeric( substr( X_score, 3, 3) )
Y_score <- "PC4" ;       Y_n <- as.numeric( substr( Y_score, 3, 3) )
axLimits <- square_plot_limits( scorePCA[ ,X_score], scorePCA[ ,Y_score] )
hnd_PC2 <- plot_Scores_scatter( scorePCA, X_score, Y_score, "Grouping", "Exp_Name", NULL, axisLimits=axLimits,
                                paste0( "PC ",X_n," (", round(expl_var[[1]][[X_n]]*100, 1) ,"%)"), 
                                paste0( "PC ",Y_n," (", round(expl_var[[1]][[Y_n]]*100, 1) ,"%)")  )
hnd_PC2
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCAScores_",grouped_at,"_",X_n,"-",Y_n,".tiff"), width=4, height=4, dpi=300)
}


X_score <- "PC5" ;       X_n <- as.numeric( substr( X_score, 3, 3) )
Y_score <- "PC6" ;       Y_n <- as.numeric( substr( Y_score, 3, 3) )
axLimits <- square_plot_limits( scorePCA[ ,X_score], scorePCA[ ,Y_score] )
hnd_PC3 <- plot_Scores_scatter( scorePCA, X_score, Y_score, "Grouping", "Exp_Name", NULL, axisLimits=axLimits,
                                paste0( "PC ",X_n," (", round(expl_var[[1]][[X_n]]*100, 1) ,"%)"), 
                                paste0( "PC ",Y_n," (", round(expl_var[[1]][[Y_n]]*100, 1) ,"%)")  )
hnd_PC3
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCAScores_",grouped_at,"_",X_n,"-",Y_n,".tiff"), width=4, height=4, dpi=300)
}




# +++ Plot LOADINGS values +++ SCATTER +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X and Y are parametrized to allow to plot any combination of two Loads
plot_loadsPCA <- loadsPCA
plot_loadsPCA <- plot_loadsPCA %>% mutate( "Taxa"=rownames(plot_loadsPCA), .before=1)
rownames(plot_loadsPCA) <- NULL

# Import microbiome ranking and relative colors. We can then factorize and plot 
# such that each feature is colored by its taxa upper_rank
upper_rank <- "Class"
ref_ranks  <- read_tsv(paste0(path_2_data, "Data_Seq/RefList_",upper_rank,".txt"),
                       col_names = FALSE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                       guess_max  = min( 10 , 300 ) )
# Take only the columns "aggr_at_rank" and the upper hierarchy to use for coloring
sub2_rank <- reduc_Taxo[ , c(upper_rank, aggr_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")

# First, merge with selection from reduc_TAXO, to march "aggr_at_rank" (e.g. Families) to colors and upper rank.
sub2_rank <- merge( sub2_rank, ref_ranks, by=upper_rank, all=TRUE)
# Second, select only unique elements of "aggr_at_rank" (since it is in long-form)
# and then merge with the plot_loadsPCA, by the "aggr_at_rank" names, to match loads and colors
temp <- unique( sub2_rank[, c(aggr_at_rank,"Colour",upper_rank)])
colnames(temp)[1] <- "Taxa"
plot_loadsPCA <- merge( temp, plot_loadsPCA, by="Taxa", all=TRUE)
plot_loadsPCA <- plot_loadsPCA[ !is.na(plot_loadsPCA$Taxa), ]

# Matching colors order with both "aggr_at_rank" and upper_rank columns is tricky;
# despite that in DF all looks good. One way is to factorize and order 
# to lock the order in which observations are plotted and the legend is organized
plot_loadsPCA[[upper_rank]] <- factor(plot_loadsPCA[[upper_rank]], levels =ref_ranks[[upper_rank]] )
plot_loadsPCA <- plot_loadsPCA[ order(plot_loadsPCA[,1]) , ]

plot_loadsPCA <- plot_loadsPCA[ complete.cases(plot_loadsPCA$Taxa), ]


X_load <- "Load_1" ;
Y_load <- "Load_2" ;
axLimits <- square_plot_limits( plot_loadsPCA[ ,X_load], plot_loadsPCA[ ,Y_load] )
hnd_L1 <- plot_Loadings_2D_arrows( plot_loadsPCA, X_load, Y_load, upper_rank, "Taxa",
                                   axisLimits=axLimits, thres_sig=0.1 )
hnd_L1 
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCA_Loads_",X_n,"-",Y_n,".tiff"), width=6, height=4, dpi=300)
}


X_load <- "Load_5" ;
Y_load <- "Load_6" ;
axLimits <- square_plot_limits( plot_loadsPCA[ ,X_load], plot_loadsPCA[ ,Y_load] )
hnd_L2 <- plot_Loadings_2D_arrows( plot_loadsPCA, X_load, Y_load, upper_rank, "Taxa",
                                   axisLimits=axLimits, thres_sig=0.1 )
hnd_L2 
if (save_subFigs){
  ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCA_Loads_",X_n,"-",Y_n,".tiff"), width=6, height=4, dpi=300)
}




# *** SUPPLEMENTARY INFORMATION ***************************************************

# +++ Plot LOADINGS as LIST +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Here we plot the ordered list of loading scores, one Component at a time.
list_Loads <- c("Load_1")                    # "Load_2", "Load_5", "Load_6"
threshold_load <- 0.0
for( jj in seq(1,length(list_Loads))) {
  sele_Load <- list_Loads[jj]
  loadsPCA[[aggr_at_rank]] <- rownames(loadsPCA)
  reference_Factors <- distinct( sub2_rank, .data[[aggr_at_rank]], .keep_all=TRUE )
  dsPlot <- merge( loadsPCA, reference_Factors, by=aggr_at_rank, .all=TRUE)
  dsPlot <- dsPlot[ order(dsPlot[[sele_Load]]), ]
  # Filter by the value of the Load
  #dsPlot <- dsPlot[ dsPlot[[sele_Load]] < -0.05 | dsPlot[[sele_Load]] > +0.05 , ]
  dsPlot <- dsPlot  %>%  filter(abs(.data[[sele_Load]]) > threshold_load |
                                abs(.data[[sele_Load]]) > threshold_load )
  dsPlot$Colour[ is.na(dsPlot$Colour) ] <- "#666666"
  # Factorize "Features"
  dsPlot[[aggr_at_rank]]  <- factor( dsPlot[[aggr_at_rank]], levels= unique(dsPlot[[aggr_at_rank]]) ) 
  
  hnd_SA <- plot_Loadings_1D_ranking( dsPlot, sele_Load, aggr_at_rank, 
                            thres_sig=threshold_load, remove_NA=TRUE )
  hnd_SA
  if (save_subFigs){
    ggsave( paste0(path_saveFig, "uB_Aggr-",aggr_at_rank,"_PCA_",sele_Load,"_ordered.tiff"), 
            width=3, height=6, dpi=300)
  }
}






# --- Ordination Analysis -------------------------------------------------------------------

# ******************************************************************************
# We need to reload the data because NMDS and PERMANOVA use Counts data, and 
# not relative data
# ******************************************************************************

path_2_data  <- "/Users/mattesa/ZenBook/R/Freeze-Thaw_Results/"
path_saveFig <- "/Users/mattesa/ZenBook/R/Freeze-Thaw_Results/Figures/" 
save_subFigs <- F

# Load ASV table(s) and Metadata
ASV  <- read.table( paste0( path_2_data,"Data_Seq/ASV_freeze.tsv"), header=T, sep="\t")
TAXO <- read.table( paste0( path_2_data,"Data_Seq/TAXO_freeze.tsv"), header=T, sep="\t")
MetaData <- read.table( paste0(path_2_data, "Data_Seq/META_freeze.txt"), header=T, sep="\t")

# +++ Aggregate ASVs at rank X:
# Import user defined aggregation function
filter_features <- T
aggr_at_rank <- "Genus"
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

# Make ASV table fully numeric (place Sample_ID as rownames) and replace colnames with taxonomic names
rownames(reduc_ASVs) <- reduc_ASVs$Sample_ID
reduc_ASVs <- reduc_ASVs[,-1]
colnames(reduc_ASVs) <- reduc_Taxo[[aggr_at_rank]]   

# NO NORMALIZATION ----- we need count data for ordination 



# --- PERMANOVA ----------------------------------------------------------------
# Store the permanova test statistics F-ration and Pr(<F)
permanova_df <- data.frame()

# Choose distance metrics to use for comparison, amongst others: 
# "bray", "jaccard", "binomial", "chao", "aitchison", "robust.aitchison"
distance_parameter <- "robust.aitchison"

# +++ DISTANCE, calculate the metric and create the distance matrix
distances <- vegan::vegdist( reduc_ASVs, method=distance_parameter )
distance_Mat <- as.matrix(distances)
colnames(distance_Mat) <- MetaData$Sample_ID
rownames(distance_Mat) <- MetaData$Sample_ID

# +++ PERMANOVA, perform test and print the test statistic results
Meta_factors <- MetaData[ , c("Exp_Name", "Exp_FTcycle")]

permanova <- vegan::adonis2( reduc_ASVs ~ Exp_FTcycle, data=Meta_factors, 
                             permutations=10000, method=distance_parameter)
permanova
# Gather the PERMANOVA test statistics results
permanova_df <- rbind ( permanova_df, data.frame( Metric=distance_parameter, Method="adoni2" ,
                                                  F_ratio=permanova[[4]][1], Pr_F=permanova[[5]][1] ) )

# +++ HOMOGENEITY, test if the groups have different variations
# if p < 0.05, then we can use distance matrix to display clustering in a  by performing a PcoA
disp_type <- "centroid"
dispersion   <- vegan::betadisper( distances, type=disp_type, group=Meta_factors$Exp_FTcycle )
permute_disp <- vegan::permutest( dispersion )
expl_var  <- as.data.frame( rbind( SD = sqrt(dispersion$eig), 
                                   Proportion = dispersion$eig/sum(dispersion$eig),
                                   Cumulative = cumsum(dispersion$eig)/sum(dispersion$eig) ) )
# Gather the homogeneity and dispersion test results
permanova_df <- rbind ( permanova_df , data.frame( Metric=distance_parameter, Method=paste0("Homogeneity_",disp_type ) , 
                                                   F_ratio=permute_disp$tab[[4]][1], Pr_F=permute_disp$tab[[6]][1] ) )

# +++ Plot PCoA scores 
scorePCoA <- as.data.frame( dispersion$vectors )
scorePCoA$Grouping <- Meta_factors$Exp_FTcycle 
scorePCoA$ID_Name  <- Meta_factors$Exp_Name 
centroidPCoA <- as.data.frame( dispersion$centroids )
centroidPCoA$Grouping <- rownames(centroidPCoA)  

# X and Y are parametrized to allow to plot any combination of two PCs
X_score <- "PCoA1" ;       X_n <- as.numeric( substr( X_score, 5, 5) )
Y_score <- "PCoA2" ;       Y_n <- as.numeric( substr( Y_score, 5, 5) )
limits  <- 25

hnd_PCoA <- ggplot( scorePCoA, aes(x=.data[[X_score]], y=.data[[Y_score]]) ) +
        geom_hline( aes(yintercept = 0), alpha = 0.4, linetype="dotted") +
        geom_vline( aes(xintercept = 0), alpha = 0.4, linetype="dotted") +
        
        # Plot Observation points and ellipse
        geom_point( aes(fill=Grouping, shape=ID_Name), alpha = 1, size = 3.2 , stroke=.4, color="#AAAAAA" ) +
        stat_ellipse(aes(color=Grouping), geom="path", level=0.95,  alpha=0.8, linewidth=0.8 ) + 
        scale_fill_manual( values = palette_FT_cycle ) + 
        scale_color_manual(values = palette_FT_cycle ) + 
        scale_shape_manual(values = c(21,22,23,24,25)) +  # shape vector for the bio ID
        
        # Plot centroid points
        geom_point( data=centroidPCoA, aes(x=.data[[X_score]], y=.data[[Y_score]]), 
                    color=palette_FT_cycle, shape=4, alpha=1, size=2.5 , stroke=1 ) +
        
        coord_fixed(ratio = 1) +
        xlim( -limits, limits) +
        ylim( -limits, limits) +
        xlab( paste0( X_score, " (", round( expl_var["Proportion",X_score]*100, 1) ,"%)") ) + 
        ylab( paste0( Y_score, " (", round( expl_var["Proportion",Y_score]*100, 1) ,"%)") ) + 
        Theme_PCA
hnd_PCoA 
if (save_subFigs){
  ggsave( paste0(path_saveFig, "PCoA-",distance_parameter,"_",aggr_at_rank,"_PC_",X_n,"-",Y_n,".tiff"), width=4, height=4, dpi=300)
}

write_tsv( permanova_df, paste0( path_2_data, "Res_Tables/Table_",aggr_at_rank,"_PERMANOVA.txt") )




# ---- NMDS --------------------------------------------------------------------
Meta_factors <- MetaData[ , c("Sample_ID","Exp_Name","Exp_FTcycle")]
# Store the permanova test statistics F-ration and Pr(<F)
nmds_df <- data.frame()

# Choose distance metrics to use for comparison, amongst others: 
# "bray", "jaccard", "binomial", "chao", "aitchison", "robust.aitchison"
distance_parameter <-"bray"

# +++ DISTANCE, calculate the metric and create the distance matrix
# The vegan::avgdist function computes the dissimilarity matrix multiple times
# using random subsampling of the the dataset each time. Then, averages 
# calculated iterations to yield the final distance matrix. 
distances <- vegan::avgdist( reduc_ASVs, dmethod=distance_parameter, sample=2500, iteration=250 )
distance_Mat <- as.matrix(distances)
colnames(distance_Mat) <- MetaData$Sample_ID
rownames(distance_Mat) <- MetaData$Sample_ID


# +++ Plot the  STRESS - Calculate the NMDS stress and plot it
stressPlot <- data.frame()
for (ii in seq(2,6)){ 
  stressPlot <- rbind( stressPlot, data.frame( Dimension = ii, Stress = unlist(
    metaMDS(distances, k=ii, model="global", try=50)["stress"]) ))
}
hnd <- ggplot( stressPlot ,aes(x=Dimension, y=Stress) ) +
  geom_point( aes(x=Dimension, y=Stress), color="#b82744", size=2) +
  geom_line(  aes(x=Dimension, y=Stress), color="#b82744", linewidth=1) +
  xlab( "NMDS Dimensions" ) + ylab( "Stress" ) + ylim(0,.12) +
  Theme_PCA
hnd
if (save_subFigs){
  ggsave( paste0(path_saveFig, "NMDS_Stress.tiff"), width=3, height=3, dpi=300)
}

# +++ Performs Non-metric Multidimensional Scaling (NMDS),
# Perform the NMDS until stable solution is found, using several random starts.
# ("try" is the number of starts before giving up)
nmds <- metaMDS( distances, k=2, model = "global", try=50) %>% scores() %>% as_tibble()
nmds$Sample_ID <- MetaData$Sample_ID
score_NMDS <- inner_join( Meta_factors, nmds )
centroid   <- score_NMDS %>% group_by(Exp_FTcycle) %>% summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

# +++ Plot NMDS scores 
limits <- 0.8
hnd_NMDS <- ggplot( score_NMDS, aes(x=NMDS1, y=NMDS2) ) +
        geom_hline( aes(yintercept = 0), alpha = 0.4, linetype="dotted") +
        geom_vline( aes(xintercept = 0), alpha = 0.4, linetype="dotted") +
        
        # Plot Observation points and ellipse
        geom_point( aes(fill=Exp_FTcycle, shape=Exp_Name), alpha = 1, size = 3.2 , stroke=.4, color="#AAAAAA" ) +
        stat_ellipse(aes(color=Exp_FTcycle), geom="path", level=0.95,  alpha=0.8, linewidth=0.8) + 
        scale_fill_manual( values = palette_FT_cycle ) + 
        scale_color_manual(values = palette_FT_cycle ) + 
        scale_shape_manual(values = c(21,22,23,24,25)) +  # shape vector for the bio ID
  
        # Plot centroid points
        geom_point( data=centroid, aes(x=NMDS1, y=NMDS2), 
                    color=palette_FT_cycle, shape=4, alpha=1, size=2.5 , stroke=1 ) +
        
        coord_fixed(ratio = 1) +
        xlim( -limits, limits) +
        ylim( -limits, limits) +
        Theme_PCA
hnd_NMDS 
if (save_subFigs){
  ggsave( paste0(path_saveFig, "NMDS_-",distance_parameter,"_",aggr_at_rank,".tiff"), width=4, height=4, dpi=300)
}


# +++ PERMANOVA, perform test and print the test statistic results
permanova <- vegan::adonis2( distances ~ Meta_factors$Exp_FTcycle, permutations=1000, 
                             method=distance_parameter )   # strata=Meta_factors$Exp_Name
permanova
# Gather the PERMANOVA test statistics results
nmds_df <- rbind ( nmds_df, data.frame( Metric=distance_parameter, Method="adoni2" ,
                                        F_ratio=permanova[[4]][1], Pr_F=permanova[[5]][1] ) )


# +++ HOMOGENEITY, test homogeneity of multivariate dispersions
# Test if there is a significant variation in beta dispersion
disp_type  <- "centroid"
dispersion <- vegan::betadisper( distances, type=disp_type, group=Meta_factors$Exp_FTcycle )
expl_var   <- as.data.frame( rbind( SD = sqrt(dispersion$eig), 
                                    Proportion = dispersion$eig/sum(dispersion$eig),
                                    Cumulative = cumsum(dispersion$eig)/sum(dispersion$eig) ) )
# Option 1 - use ANOVA
disp_anova   <- anova( dispersion )
# Option 2 - use permutation test
disp_permute <- permutest( dispersion )
# Gather the homogeneity and dispersion test results
nmds_df <- rbind ( nmds_df , data.frame( Metric=distance_parameter, Method=paste0("Homogeneity_",disp_type ) , 
                                         F_ratio=disp_permute$tab[[4]][1], Pr_F=disp_permute$tab[[6]][1] ) )

write_tsv( nmds_df, paste0( path_2_data, "Res_Tables/Table_",aggr_at_rank,"_",distance_parameter,"_NMDS.txt") )



# --- Save combined Plots ----------------------------------------------------------------------
# Place sub-plots together into one figure

ggdraw() +
  draw_plot(hnd_PC1,  x = .08, y = .50, width = .40, height = .48) +
  draw_plot(hnd_PC3,  x = .52, y = .50, width = .40, height = .48) +
  draw_plot(hnd_L2,   x = .08, y = .0,  width = .40, height = .48) +
  draw_plot(hnd_NMDS, x = .52, y = .0,  width = .40, height = .48) +
  draw_plot_label(label = c("A","B","C","D"), size = 17, 
                  x = c(.07,.51,.07,.51), y = c(.99,.99,.49,.49) )  
path_saveFig <- "/Users/mattesa/ZenBook/R/Freeze-Thaw_Results/Figures/" 
ggsave( paste0(path_saveFig, "Figure_2_Multivariate.tiff"), width=10, height=8, dpi=300)



ggdraw() +
  draw_plot(hnd_tSNE, x = .08, y = .50, width = .40, height = .48) +
  draw_plot(hnd_PC2,  x = .52, y = .50, width = .40, height = .48) +
  draw_plot(hnd_L1,   x = .08, y = .0,  width = .40, height = .48) +
  draw_plot(hnd_PCoA, x = .52, y = .0,  width = .40, height = .48) +
  draw_plot_label(label = c("A","B","C","D"), size = 17, 
                  x = c(.07,.51,.07,.51), y = c(.99,.99,.49,.49) ) 
path_saveFig <- "/Users/mattesa/ZenBook/R/Freeze-Thaw_Results/Figures/" 
ggsave( paste0(path_saveFig, "Suppl_Figure_4_Multivariate.tiff"), width=10, height=8, dpi=300)















