# +++++ DESCRIPTION - Plotting functions, repository +++++++++++++++++++++++++++
#
# This script gathers functions and common themes used for plotting the 
# results of the Freeze-Thaw paper. Thus, it is easier to maintain a unified  
# theme across all figures and analysis performed, as well as avoiding many
# lengthy repetition for plotting the many sub-figures, which essentially share
# the same plot layout.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
require(ggpattern)
require(ggplot2)
require(ggrepel)
require(ggtext)

# General plotting setting
# --- Plotting preferences
# Common options and THEME for plotting
palette_FT_cycle <- c("#00C385","#103354","#19527E","#2175A8","#299CD2","#45CEF2","#77F6FF" )
figFontSize <- 12

scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = palette_FT_cycle)
}
scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = palette_FT_cycle)
} 


main_layout <- theme(
      plot.margin = margin(t=0.25, r=0.25, b=0.25, l=0.25, "cm"),
      panel.background = element_blank() ,
      panel.grid.minor = element_blank() ,
      panel.grid.major = element_blank() , 
      plot.title   = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
      axis.title.x = element_text( color="#444444", size=figFontSize  , vjust = -.6, hjust = .5) ,
      axis.title.y = element_text( color="#444444", size=figFontSize  , vjust =  .6, hjust = .5) ,
      axis.text.x  = element_text( color="#444444", size=figFontSize-2 ) ,   
      axis.text.y  = element_text( color="#444444", size=figFontSize-2 ) ,
      # Legend
      legend.position = "none",
      legend.background = element_blank() ,
      legend.title = element_text(size = figFontSize-1), 
      legend.text  = element_text(size = figFontSize-2),
      legend.key.size = unit(0.20, 'cm'),
)


Theme_PCA <- main_layout + 
      theme(
        plot.title   = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
        panel.border = element_rect(size = 1.5, colour = "#666666", linetype = "solid", fill=NA) ,
        strip.background = element_blank() 
      )


Theme_BoxPlot <- main_layout + 
      theme(
        panel.grid.major.x=element_line( size=0.5, colour="#cbcbcb", linetype="dashed") ,
        panel.grid.major.y=element_line( size=0.5, colour="#cbcbcb", linetype="solid") ,
        axis.ticks  =element_blank(),
        axis.text.x =element_text(size=figFontSize-2, angle=45, hjust=1) ,
        axis.line.x =element_line(size=1, colour="#444444"),
        # Facet wrap options
        strip.background=element_blank(), 
        strip.placement ="outside",
        strip.text.x =element_text(size=figFontSize),  # Increase font size of facet titles
        strip.text.y =element_text(size=figFontSize),
        panel.spacing=unit(1, "lines")
      )

Theme_BarHoriz <- main_layout + 
      theme(
        panel.grid.major.x=element_line( size=0.5, colour="#cbcbcb", linetype="dashed") ,
        panel.grid.major.y=element_line( size=0.5, colour="#cbcbcb", linetype="solid") ,
        axis.ticks  =element_blank(),
        axis.text.x =element_text(size=figFontSize-5, angle=45, hjust=1) ,
        axis.line.x =element_line(size=0.6, colour="#555555"),
        # Facet wrap options
        strip.background=element_blank(), 
        strip.placement ="outside",
        strip.text.x =element_text(size=figFontSize),  # Increase font size of facet titles
        strip.text.y =element_text(size=figFontSize),
        panel.spacing=unit(1, "lines")
      )

Theme_Temp2Time <- main_layout + 
      theme(
      panel.border = element_rect(size = 1.5, colour = "#666666", linetype = "solid", fill=NA) ,
      )

Theme_ALDEx2 <- main_layout + 
      theme(
        panel.grid.major.x = element_line( size = 0.5, colour = "#cbcbcb", linetype="dashed") ,
        panel.grid.major.y = element_line( size = 0.5, colour = "#cbcbcb", linetype="dashed") ,
        plot.title  = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
        # Make title Markdown to add superscript and greek letters
        axis.title  = element_markdown( color="#444444", size=figFontSize ) ,
        axis.ticks  = element_blank(),
        axis.line.x = element_line(size = 1, colour = "#444444"),
        axis.text.x = element_text(angle = 45, hjust=1) 
      )


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

square_plot_limits <- function( vec_X, vec_Y ){
    # Given two vectors (of two variables to plot in a Cartesian coordinates)
    # find the min/max limits and then "square them. I.e., it corrects the 
    # difference (delta) between the two axes limits to ensure the plot will 
    # display a perfect square area that includes all the points.
    if ( any(abs(vec_X)<1) | any(abs(vec_Y)<1) ) {  # plotting loading, they are always in range -1 to +1
      x_MaxL <- max(vec_X) +0.1
      x_MinL <- min(vec_X) -0.1
      y_MaxL <- max(vec_Y) +0.1
      y_MinL <- min(vec_Y) -0.1
    } else{   # any other case the range will likely be several units (e.g, 3, 14, or 25)
      x_MaxL <- ceiling(max(vec_X)) +1
      x_MinL <- floor(  min(vec_X)) -1
      y_MaxL <- ceiling(max(vec_Y)) +1
      y_MinL <- floor(  min(vec_Y)) -1
    }
    delta  <- (x_MaxL+abs(x_MinL)) - (y_MaxL+abs(y_MinL))
    if( delta !=0 ){
      if( delta > 0){         
        y_MaxL <- y_MaxL+delta/2;        
        y_MinL <- y_MinL-delta/2
      }else if( delta < 1 ){  
        x_MaxL <- x_MaxL+abs(delta/2);    
        x_MinL <- x_MinL-abs(delta/2)
      }}
    return( c(x_MaxL,x_MinL, y_MaxL,y_MinL) )
}


plot_Scores_scatter <- function( dataset, xAxis, yAxis, Category_1, Category_2,
                                 title, xLabel, yLabel, axisLimits ){
  # Graph the PC scores in a scatter plot, for component X and Y. 
  # Use grouping variables: "Category_1" (as colmane) to color points and 
  #                         "Category_2" to shape the points.
  return( 
    ggplot( dataset, aes(x=.data[[xAxis]], y=.data[[yAxis]], fill=.data[[Category_1]] ) ) +
      geom_hline( aes(yintercept = 0), alpha = 0.5, linetype="dotted") +
      geom_vline( aes(xintercept = 0), alpha = 0.5, linetype="dotted") +
      geom_point( aes(shape=.data[[Category_2]]), alpha = 1, size = 3.2 , stroke=.4, color="#AAAAAA" ) +
      xlim( axisLimits[1], axisLimits[2]) +
      ylim( axisLimits[3], axisLimits[4]) +
      scale_shape_manual(values = c(21,22,23,24,25)) +  # shape vector for the bio ID
      coord_fixed(ratio = 1) +
      ggtitle( title ) + xlab( xLabel ) + ylab( yLabel ) + 
      Theme_PCA
  )
}


plot_Loadings_2D_arrows <- function( ds_loads, xAxis, yAxis, categoryVar, 
                                     feat_name, axisLimits, thres_sig=0.1 ){
  # Plot the Loading as arrows in 2D plot, for component X and Y. 
  
  X_n <- as.numeric( strsplit( X_load,"_" )[[1]][2] )   # Extract numeric value of chosen components
  Y_n <- as.numeric( strsplit( Y_load,"_" )[[1]][2] )

  # Filter data based on threshold condition
  filtered_data <- ds_loads %>% filter(abs(.data[[X_load]]) > thres_sig | abs(.data[[Y_load]]) > thres_sig )
  
  return( ggplot( ) +
            geom_hline(aes(yintercept = 0), alpha = 0.4, linetype="dashed") +
            geom_vline(aes(xintercept = 0), alpha = 0.4, linetype="dashed") +
            geom_segment( 
                data = filtered_data, 
                aes(x=0, y=0, xend=.data[[X_load]], yend=.data[[Y_load]], color=.data[[categoryVar]] ),              
                arrow = arrow(length = unit(1/2, "picas")), size=0.5, color=filtered_data$Colour  ) + 
            ggrepel::geom_text_repel( 
                data=filtered_data,
                aes(x=.data[[X_load]], y=.data[[Y_load]], label=.data[[feat_name]]), size=2, alpha=0.8,
                segment.color="#666666", segment.linetype="dashed", segment.size=0.2) +    # Set segment.color=NA
            xlim( axisLimits[2], axisLimits[1]) +
            ylim( axisLimits[4], axisLimits[3]) +
            coord_fixed(ratio = 1) + 
            xlab( paste0("Loads ", X_n) ) + 
            ylab( paste0("Loads ", Y_n) ) + 
            Theme_PCA
  )
}


plot_Loadings_1D_ranking <- function( ds_loads, selected_Load, categoryVar, 
                                      thres_sig=0.05, remove_NA=TRUE ){
  # Plot Loading for a single component, as either ordered by the implicit 
  # (knowledge based) ranking given by feat_name, or ordered by the loading
  # values (from highest to lowest)

  # Remove NA (i.e. unmeasured) metabolites if chose and filter loads by threshold value 
  # (if the load is true at any value, the load is kept and displayed)
  if ( remove_NA==TRUE ){    
    ds_loads <- ds_loads[ !is.na(ds_loads$Load_1) , ]  
  }
  thres_pass <- rowSums( ds_loads[ ,grep("Load_",colnames(ds_loads))] >= thres_sig ) != 0
  ds_loads   <- ds_loads[ thres_pass, ]
  
  # Factorize Features by loading value
  order_by_load <- order(ds_loads[,selected_Load], decreasing=T) 
  ds_loads[[categoryVar]] <- factor( ds_loads[[categoryVar]], levels=rev(ds_loads[[categoryVar]][order_by_load]) )
  
  # Factorize Features by implicit ranking found in "metabo_Ranking"
  return( ggplot( ds_loads, aes(x=.data[[categoryVar]], y=.data[[selected_Load]]) ) + 
            geom_hline(  aes(yintercept = 0), alpha = 0.5, colour="#4e4e4e") +
            geom_jitter( color="#888888", size=2.5, alpha=0.9, shape=23, stroke=0.7, fill=ds_loads$Colour) +
            coord_flip() +
            xlab("") + ylab( selected_Load ) +
            ggtitle("") +
            Theme_BoxPlot
  )
}


plot_HorizBoxPlots <- function(dataset, xAxis, yAxis, Category,  xLabel, yLabel, 
                               axLimits, axStep, colour_palette)  {
  set_axis_scale <- seq(axLimits[1], axLimits[2], by=axStep)
  return(
    ggplot(dataset, aes(x=.data[[xAxis]], y=.data[[yAxis]], fill=.data[[Category]]) ) +  
      geom_boxplot( color="#333333", alpha=0.8, outlier.shape=NA, lwd=0.1, fatten=2.5) +  
      scale_x_discrete( labels = function(x) str_wrap(as.character(x), width=33) ) +
      scale_fill_manual( values=colour_palette ) +
      xlab(xLabel)   +   ylab(yLabel) +
      scale_y_continuous( breaks=set_axis_scale, labels=set_axis_scale )  +
      coord_flip( ylim=c(axLimits[1], axLimits[2]) ) + 
      Theme_BoxPlot + 
      theme(
        axis.text.y =element_text(color="#444444", size=figFontSize-5, hjust=0),
        axis.text.x =element_text(size=figFontSize-4, angle=45, hjust=1) ,
        strip.text.x=element_text(size=7, hjust=0)
      )
  )
}





