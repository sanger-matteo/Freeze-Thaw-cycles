# +++++ DESCRIPTION --- DAA_MaAsLin2.R +++++++++++++++++++++++++++++
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
library(broom)

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
ref_TP <- "C[0123456]"         # "C[123456]"  or "C[0123456]"
mask <- stringr::str_detect( MetaData$External_ID , ref_TP )     # C1 to C6
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]
reduc_ASVs <- reduc_ASVs[ , colSums(reduc_ASVs) != 0 ]    # Remove empty features

# Normalize
reduc_ASVs <- as.data.frame( compositions::clr(reduc_ASVs) )
# reduc_ASVs <- (reduc_ASVs / rowSums(reduc_ASVs)) 

test_var <- "Exp_FTcycle"

# Find significantly different taxa with simple T-test and Wilcoxon test
wide_df <- cbind( as.vector(MetaData[test_var]), reduc_ASVs )

# Select only two time points, wince Wilcoxon test only handles two groups for comparison.
mask <- wide_df[[test_var]] %in% c("C0","C6")
wide_df <- wide_df[ mask, ]
list_features <- colnames(reduc_ASVs)

sig_Table <- data.frame()
for( rr in seq(1,length(list_features)) ) {
  returned_list <- Test_Significance( wide_df, list_features[rr], test_var, test_paired=FALSE )
  sig_Table <- rbind( sig_Table, returned_list[[1]] )
}
# Adjust p-values, accounting for the repeated test statistics
sig_Table$adj_P.Ttest <- p.adjust( sig_Table$P.Ttest, method="BH")
sig_Table$adj_P.Wilcx <- p.adjust( sig_Table$P.Wilcx, method="BH")
sig_Table$Sign     <- (sig_Table$P.Ttest<0.05) + (sig_Table$P.Wilcx<0.05)
sig_Table$adj_Sign <- (sig_Table$adj_P.Ttest<0.05) + (sig_Table$adj_P.Wilcx<0.05)


sig_Table <- merge( reduc_Taxo[,c("Phyla","Class",aggr_at_rank)], sig_Table, by.x=aggr_at_rank, by.y="Dep_var" )
sig_Table <- sig_Table %>% relocate( all_of(aggr_at_rank), .after=3 )

# Order by p-values (unadjusted) and filter significant features
sig_Table <- sig_Table[ order(sig_Table$P.Wilcx), ]

sig_features <- filter( sig_Table, P.Ttest<0.05 )
# NOTE 2: If the table comes back "empty" it is because no significantly
#         different feature was detected

write.table( sig_Table, paste0(path_2_data, "Res_Tables/Sig_Table_Wilcox_",aggr_at_rank,".txt"),
             row.names=FALSE, col.names=TRUE, sep='\t')





# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# # --- ALTERNATIVE solution:
# # Select only two time points, wince Wilcoxon test only handles two groups for
# # comparison. Then drop features that have no counts in any sample left in the cohoort
# mask <- MetaData$Exp_FTcycle %in% c("C0","C6")
# reduc_ASVs <- reduc_ASVs[ mask, ]
# reduc_ASVs <- reduc_ASVs[ , colSums(reduc_ASVs) != 0 ]
# 
# # Add metadata to ASV table and transform in long format; necessary to run
# # test statistic
# long_df <- cbind( MetaData[ mask ,c("Exp_Name","Exp_FTcycle","Sample_ID")], reduc_ASVs ) %>%
#   pivot_longer( cols = -c("Exp_Name", "Exp_FTcycle", "Sample_ID"),
#                 names_to = aggr_at_rank, 
#                 values_to = "Counts")
# selec_TAXOcolumns <- c("Phyla","Class","Order","Family","Genus")
# long_df <- merge( long_df, reduc_Taxo[ ,selec_TAXOcolumns], by=aggr_at_rank)
# 
# # Find significantly different taxa with simple Wilcoxon test
# # NOTE: We must adjust p-values, accounting for the repeated test statistics
# # NOTE 2: If the table comes back "empty" it is because no significantly 
# #         different feature was detected
# sig_features <- long_df %>% 
#                 nest(data = -all_of(selec_TAXOcolumns)) %>%
#                 mutate(test = map(data, ~ wilcox.test( Counts ~ Exp_FTcycle, 
#                                                        data=.x, 
#                                                        exact = FALSE,
#                                                        paired = TRUE, # Include pairing
#                                                        subset = !is.na(Sample_ID) 
#                                                        ) %>% 
#                                     broom::tidy() ) ) %>%
#                 unnest(test) %>%
#                 mutate(p_adjust = p.adjust(p.value, method="BH"))  # %>%
# # filter( p_adjust <0.05 )  %>%
# # select(all_of(selec_TAXOcolumns[-1]), p.value, p_adjust)
# # Order by p-values (unadjusted)
# sig_features <- sig_features[ order(sig_features$p.value), ]
# 
# # Save table
# write_tsv( sig_features, paste0( path_2_data, "Res_Tables/Table_Wilcox_significance_",aggr_at_rank,"_ConcNorm.txt") )
# 
# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




