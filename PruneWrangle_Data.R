# +++++ DESCRIPTION - Prune and Wrangle Datasets +++++++++++++++++++++++++++++++
#
# Mix of functions that re-format the microbiome feature tables or perform 
# advanced calculation using microbiome dataset.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# +++++ Aggregate ASVs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Aggregate ASVs at a specified taxonomic rank; i.e., summing up column features
#  that share the same name at the chosen taxonomic rank.
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

aggregate_at_Rank <- function( ds_ASV, ds_taxo, aggregate_at_rank, reindex_ASV=TRUE, 
                               Taxon_as_Colname=FALSE, add_RegionName=FALSE ){
  # Aggregate all the ASVs (rowSums) that belong to the same unique taxon,
  # thus, reducing the features to only the unique taxa of the rank of choice 
  # ("aggreg_at_rank").
  # The script also reduce the TAXO table to contain only the taxonomy that is 
  # relevant for the reduced ASV table. 
  #
  # Because the script keeps in reduce TAXO the first valid row nomenclature that 
  # is common to all the rows sharing the nomenclature until the aggreg_at_rank:
  # - a - "ASV_XYZ" names of feature are still unique, but are somewhat useless
  #        because the features are now aggregates of multiple ASVs
  # - b -  the taxonomic nomenclature after "aggreg_at_rank" are to be ignored;   
  #        the names are carried over and kept just to maintain TAXO table with the 
  #        same dataframe dimensions
  #
  # Intput: 
  # - ASVs:   First column are "Sample_ID"; remaining columns are features and 
  #           named by "ASV_ID"; rows are observations.
  # - TAXO:   each row is the full taxonomy of an ASV. First, is the "ASV_ID"
  #           in the form of "ASV_XXXXX", then the taxonomy ranking:
  #           "Kingdom","Phyla","Class","Order","Family","Genus","Species"
  # Output: 
  # - reduc_ASVs: Same as ASVs
  # - reduc_TAXO: Same as TAXO, but smaller, containing only the first feature
  #               features
  #
  # +++ Run as SOURCECODE within other scripts:
  # +++  > source(paste0( path_workdir, "pruning_MicrobiomeData.R"))
  
  uniq_names <- unique( ds_taxo[[aggregate_at_rank]] )
  
  for (kk in seq(1,length(uniq_names)) ){
    # Find the features whose name is the same at the specified rank
    idx_feats <- ds_taxo[ ds_taxo[[aggregate_at_rank]]==uniq_names[kk], "ASV_ID"]
    new_feat_name <- idx_feats[1]
    
    # 1st IF --- at first iteration, initialize storing variable reduc_taxo and 
    # reduc_feat; otherwise row.bind (for taxo) or col.bind (for ASV)
    # 2nd IF --- if there is only one ASV matching the list_SPP[kk], add it; 
    # otherwise there many ASVs and you need to rowSums them all into one column.
    #
    # In any case, take the first ASV ID found to match uniq_names[kk] (i.e. 
    # idx_feats[1]) as reference taxonomy to create new taxonomic table and to
    # rename the newly added feature column of reduc_feat (use new_feat_name)
    if (kk == 1){  
      aggreg_taxo <- ds_taxo[ ds_taxo$ASV_ID==new_feat_name, ]
      if(length( which(colnames(ds_ASV)%in%idx_feats) )==1){   # On first iteration include Sample_ID names in ASV table
        aggreg_feat <- data.frame( Sample_ID=ds_ASV[,1] ,
                                   ds_ASV[ , colnames(ds_ASV)%in%idx_feats] )
      }else{
        aggreg_feat <- data.frame( Sample_ID=ds_ASV[,1] ,
                                   rowSums( ds_ASV[ , which(colnames(ds_ASV)%in%idx_feats)], na.rm=T ) )
      }
      colnames(aggreg_feat)[2] <- new_feat_name
      
    } else { 
      aggreg_taxo <- rbind( aggreg_taxo, ds_taxo[ ds_taxo$ASV_ID==new_feat_name, ] )
      if(length( which(colnames(ds_ASV)%in%idx_feats) )==1){
        temp <- data.frame( ds_ASV[ , which(colnames(ds_ASV)%in%idx_feats)] )
      }else{
        temp <- data.frame( rowSums( ds_ASV[ , which(colnames(ds_ASV)%in%idx_feats)], na.rm=T ) )
      }
      colnames(temp) <- new_feat_name
      aggreg_feat <- cbind( aggreg_feat, temp )
    }
  }
  
  # Re-index the ASV_ID with or without region VXVY tag
  if(reindex_ASV==T){
    if(add_RegionName==F | add_RegionName==""){ 
      aggreg_taxo$ASV_ID <- sprintf( "ASV_%05d", seq(1,nrow(aggreg_taxo)) ) 
    }else{             
      aggreg_taxo$ASV_ID <- sprintf( "ASV_%s_%05d", add_RegionName, seq(1,nrow(aggreg_taxo)) )
    }
    colnames(aggreg_feat)[-1] <- aggreg_taxo$ASV_ID
  }
  # Re-index row numbers
  rownames(aggreg_taxo) <- seq(1,nrow(aggreg_taxo))
  
  # Assign to column names of ASV table the taxon name at specified rank
  if(Taxon_as_Colname==T){
    colnames(aggreg_feat)[-1] <- aggreg_taxo[[aggregate_at_rank]]
  }
  
  return( list(aggreg_feat, aggreg_taxo) )
}





# +++++ Calculate Significance +++++++++++++++++++++++++++++++++++++++++++++++++
#
# Perform several significant tests to find if there are statistical differences 
# between compared groups
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Test_Significance <- function( dataset, dep_var, indep_var, test_paired=TRUE ) {
  # Test statistics provided dependent variable (dep_var) using observations 
  # groups from indepedendent variable (indep_var). Then extract significance 
  # of the tests (i.e., p-values):
  #
  # 1) Ind_var group == 2    - t-test
  #                          - Wilcoxon test
  # 2) Ind_var group  > 2    - ANOVA  (+ normality and uniformity test)
  #                          - Kruskal-Wallis
  #
  # Return a ready made list of p-value table, and the complete results of each 
  # as a list
  
  res_List <- list()
  res_Pval <- data.frame( "Dep_var"=dep_var, "Ind_var"=indep_var, 
                          "N-groups"=length(unique(dataset[[indep_var]])) )
  
  # We have to ensure that column names have no " " or "-" characters, otherwise
  # It is not possible to create a formula
  ori_colnames <- colnames(dataset)
  new_colnames <- gsub( " ", "__", ori_colnames)
  new_colnames <- gsub( "-", "..", new_colnames)
  dep_var   <- gsub( " ", "__", dep_var)
  dep_var   <- gsub( "-", "..", dep_var)
  indep_var <- gsub( " ", "__", indep_var)
  indep_var <- gsub( "-", "..", indep_var)
  colnames(dataset) <- new_colnames
  
  # Define the main formula for the dependent and independent variables
  main_formula <- paste0(dep_var," ~ ",indep_var)
  
  # Case == 2 groups 
  # performT-test and Wilcoxon test
  if (length(unique(dataset[[indep_var]])) == 2) {
    if (test_paired==T){
      group <- unique(dataset[[indep_var]])
      ttest_result <- t.test( dataset[ dataset[[indep_var]]==group[1], dep_var][[1]],
                              dataset[ dataset[[indep_var]]==group[2], dep_var][[1]],  paired=T )
      wilx_result  <- wilcox.test( dataset[ dataset[[indep_var]]==group[1], dep_var][[1]],
                                   dataset[ dataset[[indep_var]]==group[2], dep_var][[1]],  paired=T )
    }else{
      ttest_result <- t.test( as.formula(main_formula), data=dataset )
      wilx_result  <- wilcox.test( as.formula(main_formula), data=dataset )
    }
    res_List$t_test   <- ttest_result
    res_List$wilcoxon <- wilx_result
    res_Pval <- cbind( res_Pval, data.frame( P.Ttest=ttest_result$p.value,
                                             P.Wilcx=wilx_result$p.value ))
    
  }else {
    # Case  > 2 groups
    # --- ANOVA test (if p-value < 0.05, then reject H0: µ are equal)
    aov_fit <- aov( as.formula(main_formula), data=dataset)
    # -> Assumption 1 - CHECK NORMALITY of variable(S) - Shapiro test (if p-value < 0.05, then reject H0: distribution is normal) 
    shapiro_table <- shapiro.test( aov_fit$residual )
    Shapiro=shapiro_table$p.value
    # -> Assumption2 - CHECK VARIANCE uniformity - Levene test (if p-value < 0.05, then reject H0: variances are equal) 
    levene_table <- car::leveneTest( as.formula(main_formula) , data=dataset )
    
    # --- Kruskal-Wallis test
    # ALTERNATIVE, for small cohort and if normality or variance assumptions are not met
    # if p-value < 0.05, then reject H0: the means are equal
    kruskal_result  <- kruskal.test( as.formula(main_formula), data=dataset)
    
    res_List$anova  <- summary(aov_fit)
    res_List$kruskal_wallis <- kruskal_result
    res_Pval <- cbind( res_Pval, data.frame( P.Anova=res_List$anova[[1]][1,5],
                                             Shapiro=shapiro_table$p.value,
                                             Levene =levene_table[1,3],
                                             P.KruskalWallis=res_List$kruskal_wallis$p.value ))
  }
  
  res_Pval <- as.data.frame( lapply( res_Pval, function(x) if(is.numeric(x)) round(x, 6) else x) )
  
  return( list(res_Pval, res_List) )
}




Test_Significance_Metrics <- function( dataset, dep_var, indep_var, pair_category, row_name ) {
  # Test statistics provided dependent variable (dep_var) using observations 
  # groups from indepedendent variable (indep_var). This function is to be used 
  # specifically with "metrics" data (script Figure_2_Metrics); it will perform 
  # many tests, such as:
  # - ANOVA,
  # - Kruskal Wallis test
  # - Friedman test
  # - Pairwise Wilcoxon
  # - Pairwise T-test
  # Then all the tests results are save a one-row dataframe, which can be 
  # appended (rbind())) when repeated tests are performed in Figure_2_Metrics.R 
  # script, thus, creating one summary table of all test conditions and variables.
  
  #  +++ ANOVA test ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # (if p-value < 0.05, then reject H0: µ are equal)
  fit_aov <- aov( dataset[[dep_var]] ~ dataset[[indep_var]] + dataset[[pair_category]] )
  anova_table <- summary(fit_aov) 
  
  # --> Assumption 1 - CHECK NORMALITY of variable(S) - Shapiro test (if p-value < 0.05, then reject H0: distribution is normal) 
  residuals <- data.frame( Res=fit_aov$residuals)   
  shapiro_table <- shapiro.test( residuals$Res )
  
  # --> Assumption2 - CHECK VARIANCE uniformity - Levene test (if p-value < 0.05, then reject H0: variances are equal) 
  levene_table <- car::leveneTest( dataset[[dep_var]] ~ dataset[[indep_var]] )
  
  
  # +++ Kruskal Wallis test ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # ALTERNATIVE, for small cohort and if normality or variance assumptions are not met
  # if p-value < 0.05, then reject H0: the means are equal
  kruskal_t1 <- kruskal.test(  dataset[[dep_var]] ~ dataset[[indep_var]] )
  kruskal_t2 <- kruskal.test(  dataset[[dep_var]] ~ dataset[[pair_category]] )
  FSA::dunnTest( dataset[[dep_var]] ~ dataset[[pair_category]], data = dataset,  method="bh" )
  
  
  # +++ Friedman test ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # A non-parametric alternative to one-way repeated measures  ANOVA test.
  # Use to assess whether there are any stat. sig. differences between 
  # distributions of three or more paired groups. Recommended when the normality
  # assumptions is not met or when the dependent variable is measured on ordinal scale.
  my_formula <- as.formula( paste(dep_var," ~ ",indep_var," | ",pair_category ) )
  fried_test <- friedman_test( dataset, my_formula )
  # Effect size: The Kendall’s W coefficient assumes the value from 0 
  #     (indicating no relationship) to 1 (indicating a perfect relationship).
  #     Kendall’s W uses the Cohen’s interpretation guidelines of:
  #     -  0.1 - < 0.3 (small effect), 
  #     -  0.3 - < 0.5 (moderate effect) and 
  #     -  >= 0.5 (large effect). 
  #     Confidence intervals are calculated by bootstap.
  my_formula <- as.formula( paste(dep_var," ~ ",indep_var," | ",pair_category ) )
  fried_eff <- friedman_effsize( dataset, my_formula )    # ci=TRUE )
  
  
  # +++ Pairwise Wilcoxon signed-rank test +++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  # P-values are adjusted using the BH multiple testing correction method.
  my_formula <- as.formula( paste(dep_var," ~ ",indep_var ) )
  Wilx <- wilcox_test( dataset, my_formula, paired=TRUE, p.adjust.method="BH")
  Wx_colnames <- paste0(Wilx$group1, "_", Wilx$group2, "__Tt_P")
  
  
  # +++ Pairwise t-test signed-rank test +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  # P-values are adjusted using the BH multiple testing correction method.
  my_formula <- as.formula( paste(dep_var," ~ ",indep_var ) )
  Tt_pair <- pairwise_t_test(dataset, my_formula, paired=TRUE, p.adjust.method="BH")
  Tt_colnames <- paste0(Tt_pair$group1, "_", Tt_pair$group2, "__Tt_P")
  
  
  # +++ Save and return data
  test_summary <- data.frame( Method =row_name,                         Df     =anova_table[[1]][1,1], 
                              Sum_Sq =anova_table[[1]][1,2],            Mean_Sq=anova_table[[1]][1,3], 
                              F_value=anova_table[[1]][1,4],            Pr_F   =anova_table[[1]][1,5],
                              Shapiro=shapiro_table$p.value,            Levene =levene_table[1,3], 
                              KW_test_F1=kruskal_t1$p.value ,           KW_test_F2=kruskal_t2$p.value,
                              Friedman_p=fried_test$p       ,           Fried_EffSize=fried_eff$effsize,
                              # Fried_Eff_conf.low=fried_eff$conf.low,    Fried_Eff_conf.high=fried_eff$conf.high,
                              Wilcox_any_05=any(Wilx$p <0.05),          Wilcox_any_padj_05=any(Wilx$p.adj <0.05),
                              T_test_any_05=any(Tt_pair$p <0.05),       T_test_any_padj_05=any(Tt_pair$p.adj <0.05))
  return(test_summary)
}



Test_LMM_Metrics <- function( dataset ) {
  # Run a Linear Mixed-Effect Model, with the specific variables:
  #  - Values  :   measurements (dependent)
  #  - FT_cycle:   time points (independent)
  #  - RowId   :   random effect
  # Then do post-hoc pairwise comparisons using estimated marginal means (EMMs) 
  # with FDR adjustment to reveal significant differences between each
  # possible pair of FT cycles
  
  LMM_model <- lmerTest::lmer(Values ~ FT_cycle + (1 | RowID), data = dataset)
  # anova(LMM_model)   # display p-value of LMM
  LMM_Pval <- anova(LMM_model)[[6]]

  LMM_sig   <- emmeans::emmeans( LMM_model, pairwise ~ FT_cycle, adjust="fdr")
  LMM_pairwise <- as.data.frame( LMM_sig$contrasts )
  LMM_pairwise <- LMM_pairwise %>% mutate( Sig= case_when(  p.value>=0.05 ~ "", p.value<0.05 ~ "*", p.value<0.01 ~ "**" ) )
  
  return( list( LMM_Pval, LMM_pairwise) )
}






