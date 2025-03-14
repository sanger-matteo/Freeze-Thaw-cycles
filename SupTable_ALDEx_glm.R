# +++++ DESCRIPTION --- DAA_ALDEx2.R +++++++++++++++++++++++++++++
#
# Differential Abundance Analysis using the algorithm:
# ALDEx2 - ANOVA-Like Differential Expression.
#          https://bioconductor.org/packages/release/bioc/html/ALDEx2.html
#          https://www.bioconductor.org/packages/devel/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html
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
library(readr)
library(tidyr)
library(dplyr)
library(ggtext)     # for markdown writing
library(cowplot)
library(ALDEx2)

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
timepoint1 <- "C0"
timepoint2 <- "C6"

ref_TP <- "C[0123456]"    # "C[123456]"  or "C[0123456]"
mask <- stringr::str_detect( MetaData$External_ID , ref_TP )     # C1 to C6
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]
reduc_ASVs <- reduc_ASVs[ , colSums(reduc_ASVs) != 0 ]    # Remove empty features



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Perform ALDEx2 - All approaches
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# You can run ALDEx2 in three ways. The results will be essentially the same.
# Approch 1 is more straightforward and do not return intermediate results.
# OPT_1 - aldex() calls aldex.clr, aldex.ttest, aldex.effect in turn and then
#         merges the data into one dataframe called x.all.
# x.aldex <- aldex( count_ds, condition_test, mc.samples=128, test="t", effect=TRUE, 
#                   include.sample.summary=FALSE, denom="all", verbose=T, paired.test=FALSE, gamma=NULL)
# NOTE: One of the most powerful aspect of ALDEx2 is that everything is
#       calculated on posterior distributions. The underlying posterior
#       distribution can be viewed if you use the individual modules (see next).
# 
# OPT_2 - Run aldex.clr(), aldex.ttest(), aldex.effect() individually and
#         and calculate. This allows more in-depth exploration of the data.
#         In addition, we obtain all the intermediate values.
# 
# x.clr    <- aldex.clr( count_ds, condition_test, mc.samples=2500, denom="all", verbose=T)
# x.tt     <- aldex.ttest( x.clr, hist.plot=T, paired.test=F, verbose=T )
# x.effect <- aldex.effect(x.clr, CI=T, verbose=F, include.sample.summary=F, 
#                          paired.test=FALSE, glm.conds=NULL, useMC=F)
# x.all <- data.frame(x.tt,x.effect)
# 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# # -- PLOT A -- Create default plots
# par(mfrow=c(1,3))
# aldex.plot(x.aldex, type="MA", test="welch", xlab="Log-ratio abundance",
#            ylab="Difference", main='Bland-Altman plot')
# aldex.plot(x.aldex, type="MW", test="welch", xlab="Dispersion",
#            ylab="Difference", main='Effect plot')
# aldex.plot(x.aldex, type="volcano", test="welch", xlab="Difference",
#            ylab="-1(log10(q))", main='Volcano plot') 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# OPT 3 - to model complex study design we can use GLM. 
# First we must create a dynamic design matrix using independnet variables.
count_ds <- as.matrix( t( reduc_ASVs) )        # NOTE: Transpose because ALDEx2 wants observations as columns
test_var <- "Exp_FTcycle"
formula <- as.formula( paste0("~ ",test_var," + (1|Exp_Name)") )
MetaData$Exp_Name <- as.numeric( gsub( "A", "", MetaData$Exp_Name ) )
design_matrix <- model.matrix(formula, data=MetaData )

x.glm    <- aldex.clr( count_ds, design_matrix, mc.samples=300, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, design_matrix, fdr.method='BH') 
glm.eff  <- aldex.glm.effect(x.glm)
glm.all  <- data.frame(glm.test,glm.eff)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# GLM results from ALDEx2 automatically creates column names of effect
# size, p.values, error within, etc; for each fixed and random variable
# used in the model. With grepl() we can find the columns of interest
# (essentially test_var) and save the data in a summary table
significance_col <- paste0("Exp_FTcycle",timepoint2,".pval")
glm.all <- glm.all[ order( glm.all[ ,significance_col ] ),  ]

glm_names <- colnames(glm.all)
col_selection <- c( "Feature", 
                    glm_names[ grep( ".*\\.pval$", glm_names) ],
                    glm_names[ grep( ".*\\.pval.padj$", glm_names) ], 
                    glm_names[ grep( ".*\\.effect$", glm_names) ] )

summary_SigTable <- filter( glm.all, .data[[significance_col]]<0.05 )
summary_SigTable$Feature <- rownames(summary_SigTable)
summary_SigTable <- summary_SigTable[ col_selection ]

if ( nrow(summary_SigTable)!=0 ){
  summary_SigTable <- mutate( summary_SigTable, Rank=aggr_at_rank, .before=1 )
  summary_SigTable <- mutate( summary_SigTable, Formula=as.character(formula[2]), .after=2 )
  summary_SigTable <- mutate( summary_SigTable, TestVar=test_var, .after=3 )
}
# Add taxonomic information
colnames(summary_SigTable)[2] <- aggr_at_rank
summary_SigTable <- merge( reduc_Taxo[,c("Phyla","Class",aggr_at_rank)], summary_SigTable, by=aggr_at_rank  )
summary_SigTable <- summary_SigTable %>% relocate( all_of(c("Phyla","Class",aggr_at_rank)), .before=3 )

# Save the table with significant features
write_tsv(summary_SigTable, paste0( path_2_data, "Res_Tables/Sig_Table_ALDEx2_",aggr_at_rank,".txt") )



# -- PLOT B -- Create custom plots with ggplot, which allows us to color the 
#              features by rank "Class" and better visualize the results
# Log2 transform the dataset and replace the "Inf" (count =0) with 0
post_hoc_dataset <- log2( count_ds[ rownames(glm.all), ] )   # hist( post_hoc_dataset[ , 1], main="raw counts" )
post_hoc_dataset[post_hoc_dataset == -Inf] <- 0
post_hoc_dataset[post_hoc_dataset == Inf ] <- 0

mu_A <- apply( post_hoc_dataset[ , grepl( timepoint1, colnames(post_hoc_dataset)) ] , 1, function(x) mean(x) )
mu_B <- apply( post_hoc_dataset[ , grepl( timepoint2, colnames(post_hoc_dataset)) ] , 1, function(x) mean(x) )

var_A <- apply( post_hoc_dataset , 1, function(x) var(x)/5 )
var_B <- apply( post_hoc_dataset , 1, function(x) var(x)/5 )

t_statistic <- (mu_A-mu_B) / (sqrt(var_A+var_B))

glm.all$t_statistic <- t_statistic
glm.all$Feature <- rownames(glm.all)

dsPlot <- merge( reduc_Taxo[,-c(1,2,3)], glm.all, by.x=aggr_at_rank, by.y="Feature", all=TRUE)



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Set plot colors
# Define palette and factorization taxa nomenclature to the reference_Taxa table
upper_rank <- "Class"
reference_Taxa <- as.data.frame( read_tsv( paste0(path_2_data, "Data_Seq/RefList_",upper_rank,".txt") ))
dsPlot <- merge( dsPlot, reference_Taxa[c("Feature","Colour")], by.x=upper_rank, by.y="Feature")

dsPlot$Class <- factor( dsPlot$Class, levels=reference_Taxa$Class )
dsPlot$Genus <- factor( dsPlot$Genus, levels=dsPlot[ order(dsPlot$Class), "Genus"] )
# Set palette as list of combinations colour==name
my_palette <- setNames( dsPlot$Colour, dsPlot$Genus )


# We can only show one timepoint at a time. Select which to show in the MA and MW plots
select_TP <- timepoint2
col_rab.all  <- glm_names[  grep( paste0("^",test_var,".*\\.rab.all$"), glm_names) ]
col_rab.all  <- col_rab.all[ grepl( select_TP, col_rab.all) ]
col_diff.btw <- glm_names[   grep( paste0("^",test_var,".*\\.diff.btw$"), glm_names) ]
col_diff.btw <- col_diff.btw[ grepl( select_TP, col_diff.btw) ]
col_diff.win <- glm_names[   grep( paste0("^",test_var,".*\\.diff.win$"), glm_names) ]
col_diff.win <- col_diff.win[ grepl( select_TP, col_diff.win) ]
col_we.eBH   <- glm_names[   grep( paste0("^",test_var,".*\\.we.eBH$"), glm_names) ]
col_we.eBH   <- col_we.eBH[ grepl( select_TP, col_we.eBH) ]
col_effect   <- glm_names[   grep( paste0("^",test_var,".*\\.effect$"), glm_names) ]
col_effect   <- col_effect[ grepl( select_TP, col_effect) ]

# Find features that had a significant change in distribution: look for we.eBH < 0.05
dsPlot$significants <- "no"
mask_sig <- which(dsPlot[col_we.eBH] <0.05 )
if ( length(mask_sig)!=0 ){ dsPlot$significants[ mask_sig] <- "yes" }

# Use rel. abu. to determine those with values above or below threshold 0
# Then set the accordingly a specific line-colour
dsPlot$relAbu_thres <- "grey"
dsPlot$relAbu_thres[ which(dsPlot[[col_rab.all]] <0) ] <- "black"
dsPlot$relAbu_thres[ mask_sig ] <- "darkred"
dsPlot$relAbu_thres <- factor(dsPlot$relAbu_thres, levels=c("grey","black","darkred"))



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# +++ Bland-Altman plot (or MA plot) +++++++
hnd_A <-ggplot( dsPlot, aes( x=.data[[col_rab.all]], y=.data[[col_diff.btw]],
                             fill=.data[[aggr_at_rank]], shape=relAbu_thres, color=relAbu_thres))+
        # Draw threshold line
        geom_vline( aes(xintercept=0), color="#888888", alpha=0.7, linetype="dashed") +   # ratio phi =2
        geom_point( size=2.7, alpha = 0.8, stroke=0.75) +  # MW plot
        scale_shape_manual( values = c(21,23,24) ) +
        scale_color_manual( values = c("grey","black","darkred") ) +
        scale_fill_manual(  values = my_palette ) +
        # ggrepel::geom_label_repel( aes(label=Feature), size=2, box.padding=unit(0.25, "lines"), point.padding=unit(0.2, "lines")) +
        # Use markdown axis labels to add subscription (see theme_Plot)
        labs( x="Median Log<sub>2</sub> Rel. Abundance" ,            # between-group mean of feature
              y="Median Log<sub>2</sub> Difference [&delta;]" ) +    # between-group mean difference of feature
        Theme_ALDEx2
hnd_A

# +++ Effect plot (or MW plot) +++++++++++++
# Select subset of features based on effect value, and only for these we show the name in the chart next to the points
Taxa_Subset <- dsPlot %>% filter( col_effect < -0.50 | col_effect > +0.50 )
hnd_B <-ggplot( dsPlot, aes(  x=.data[[col_diff.win]], y=.data[[col_diff.btw]],
                              fill=.data[[aggr_at_rank]], shape=relAbu_thres, color=relAbu_thres)) +
        # Draw diagonal line, to set ratios of difference/dispersion (effect size)
        geom_abline(aes(intercept=0), slope= +1, color="#666666", alpha=0.4, linetype="dotted") +     # ratio phi =1
        geom_abline(aes(intercept=0), slope= -1, color="#666666", alpha=0.4, linetype="dotted") +
        geom_abline(aes(intercept=0), slope= +2, color="#666666", alpha=0.5, linetype="dashed") +     # ratio phi =2
        geom_abline(aes(intercept=0), slope= -2, color="#666666", alpha=0.5, linetype="dashed") +
        geom_abline(aes(intercept=0), slope= +3, color="#666666", alpha=0.6, linetype="longdash") +   # ratio phi =3
        geom_abline(aes(intercept=0), slope= -3, color="#666666", alpha=0.6, linetype="longdash") +

        geom_point( size=2.7, alpha = 0.8, stroke=0.75) +  # MW plot
        scale_shape_manual( values = c(21,23,24) ) +
        scale_color_manual( values = c("grey","black","darkred") ) +
        scale_fill_manual(  values = my_palette) +
        ggrepel::geom_label_repel( data=Taxa_Subset, aes(label=.data[[aggr_at_rank]]), color="#444444", alpha=0.7, fill="#FFFFFF", segment.color ="#444444",
                                   size=2.5, force=80, min.segment.length=0, box.padding=0.1, max.overlaps=15 ) +
        coord_cartesian(xlim = c(0.65,NA)) +
        # Use markdown axis labels to add subscription (see theme_Plot)
        labs( x="Median Log<sub>2</sub> Dispersion [&sigma;]",     # dispersion, as pooled standard deviation
              y="Median Log<sub>2</sub> Difference [&delta;]" ) +  # difference between the means
        Theme_ALDEx2
hnd_B



# --- Combine chats in one figure ---------------------------------------------------------------
ggdraw() +
  draw_plot(hnd_A, x = .0, y = .0,  width = .45, height = .95) +
  draw_plot(hnd_B, x = .5, y = .0,  width = .45, height = .95) +
  draw_plot_label(label = c("A","B"), size = 17, x = c(.02,.52), y = c(1,1) )
path_saveFig <- "/Users/mattesa/ZenBook/R/Freeze-Thaw_Results/Figures/" 
ggsave( paste0(path_saveFig, "Suppl_Figure_3_GLM",select_TP,".tiff"), width=10, height=5, dpi=300)


hnd_B
ggsave( paste0(path_saveFig, "Suppl_Figure_3_GLM_MW",select_TP,".tiff"), width=5, height=5, dpi=300)




