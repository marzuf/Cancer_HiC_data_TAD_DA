# zscoreDT = NULL  
# ranked_branches =T
# plotMap = "square" 
# low_values_col = "#FFFFFF" 
# low_limit_col = NULL
# high_limit_col = NULL
# high_values_col ="#FF0000" 
# na_values_col = "white"
# fill_legName = "" 
# dendroLabSize = 4
# # dendroLabPosition = 0.2
# # dendroLeftWidth = 0.5
# # dendroRightWidth = 0.5
# # ylab_size = 1
# annotateMat = TRUE
# annotateMean = TRUE
# comparisonName = "comparison"
# lab_color_vect = NULL
# ylab_symbols = NULL
# lowerTri = TRUE
# addClusterDot = FALSE
# legCategoryCols = NULL
# labFontFace="bold"
# 
# 
# 
# x=as.matrix(corMat)
# ranked_branches =T
# plotMap = "square" 
# low_limit_col = 0
# high_limit_col = 1
# fill_legName = "MoC" 
# dendroLabSize = 4
# addClusterDot = TRUE
# annotateMat = TRUE
# annotateMean = TRUE
# comparisonName = "caller"
# legCategoryCols = NULL
# lab_color_vect = NULL



plot_ggheatmap_with_left_rowdendro <- function(x,
                      zscoreDT = NULL,  
                      ranked_branches =T,
                      plotMap = "square", 
                      low_values_col = "#FFFFFF", 
                      low_limit_col = NULL,
                      high_limit_col = NULL,
                      high_values_col ="#FF0000", 
                      na_values_col = "white",
                      fill_legName = "", 
                      dendroLabSize = 4,
                      # dendroLabPosition = 0.2,
                      # dendroLeftWidth = 0.5,
                      # dendroRightWidth = 0.5,
                      # ylab_size = 1,
                      annotateMat = TRUE,
                      annotateMean = TRUE,
                      comparisonName = "comparison",
                      lab_color_vect = NULL,
                      ylab_symbols = NULL,
                      lowerTri = TRUE,
                      addClusterDot = FALSE,
                      legCategoryCols = NULL,
                      labFontFace="bold",
                      perso_rightLab = NULL,
                      perso_rightLeg = NULL
                      ) {

suppressPackageStartupMessages(library(dendextend, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(cowplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


axisLabelSize <- 4
legTextSize <- 4
legTitSize <- 4
cellSize <- 4

  if(!is.null(zscoreDT)) {
    stopifnot("MoC" %in% colnames(zscoreDT))
    stopifnot("comparison" %in% colnames(zscoreDT))
    stopifnot("zscores" %in% colnames(zscoreDT))
  }
  if(is.null(colnames(x)))
    colnames(x) <- sprintf("col%s",1:ncol(x))
  if(is.null(rownames(x)))
    rownames(x) <- sprintf("row%s",1:nrow(x))
  ## plot a heatmap
  ## x is an expression matrix
  
  ##########################################################################################################
  ##########################################################################################################  DENDROPLOT: DENDROPLOT + YLAB
  ##########################################################################################################
  
  row.hc <- hclust(dist(x), "ward")
  col.hc <- hclust(dist(t(x)), "ward")
  row.dend <- as.dendrogram(row.hc)
  row.dend <- sort(row.dend, decreasing = FALSE)
  if(! exists("outfile_dendro"))
    outfile_dendro <- file.path(outFold, "row.dend_check.png")
  png(outfile_dendro)
  plot(row.dend)
  foo <- dev.off()
  
  if(ranked_branches) row.dend <- rank_branches(row.dend)

  row.dend_save <- row.dend
    
  if(!is.null(lab_color_vect)) {
    lab_cols <- unlist(sapply(as.character(labels(row.dend)), function(x) {
      if(x %in% names(lab_color_vect)) return(as.character(lab_color_vect[x]))
      return("black")
    }))
  } else {
    lab_cols <- rep("black", length(labels(row.dend)))
  }
  labels(row.dend) <- paste0(paste0(rep(" ", 3), collapse = ''), labels(row.dend))
  if(!is.null(ylab_symbols)) {
    if(length(ylab_symbols) == length(labels(row.dend))) {
      lab_labels <- ylab_symbols
      labels(row.dend) <- ylab_symbols
    } else{
      lab_labels <- rep(ylab_symbols[1], length(labels(row.dend)))
      labels(row.dend) <- rep(ylab_symbols[1], length(labels(row.dend)))
    }
  }
  
  caller_labels <- gsub(" +", "", labels(row.dend))
  
  row.dend_save_lab <- row.dend
  
  # PLOT EMPTY LABELS
  if(addClusterDot) {
    labels(row.dend) <- rep("  \u25cf", length(labels(row.dend)))
  } else {
    labels(row.dend) <- rep("", length(labels(row.dend)))
  }
  
  lab_cex <- 1.2
  # gg_row_dendro <- ggplot(row.dend %>%
  #                  set('branches_col', 'black') %>%
  #                  set('branches_lwd', 0.6) %>%
  #                  set('labels_colors', rev(lab_cols)) %>%
  #                  set('labels_cex', lab_cex) %>%
  #                  rotate(rev(1:length(labels(row.dend)))),
  #                theme = theme_minimal(),
  #                horiz = TRUE)
  gg_row_dendro <- ggplot(row.dend %>%
                            set('branches_col', 'black') %>%
                            set('branches_lwd', 0.6) %>%
                            #set('labels_colors', rev(lab_cols)) %>%
                            set('labels_colors',lab_cols) %>%
                            set('labels_cex', lab_cex),# %>%
                            # rotate(rev(1:length(labels(row.dend)))),
                          theme = theme_minimal(),
                          horiz = TRUE)
  gg_row_dendro <- gg_row_dendro + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.text = element_blank(),
                       axis.text.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.title = element_blank())

  gg_row_dendro <- gg_row_dendro + theme(plot.margin = unit(rep(0, 4), "lines"))

  
  gg_row_dendro <- gg_row_dendro + scale_x_reverse(expand=c(0.5/length(caller_labels),0))
    

  gg_row_dendro$labels$hjust <- 1
  gg_row_dendro$labels$vjust <- 1
    
  labelsDT <- data.frame(x=1, y = 1:length(caller_labels), caller = rev(caller_labels))
  # labelsDT <- data.frame(x=1, y = 1:length(caller_labels), caller = caller_labels)
  
  if(addClusterDot) {
    labels_col <- "black"
  } else {
    labels_col <- rev(lab_cols)
  }
  # labels_col <- rev(lab_cols)
  
  gg_y_dendrolab <- ggplot(labelsDT, aes(x=x, y=y)) +
    theme_transparent()+
    scale_y_continuous(expand=c(0.5/length(caller_labels),0))+
    geom_text(aes(x-0.5, y, label =caller), 
              #size=dendroLabSize,
                size = axisLabelSize,
                hjust=1, vjust=0.5,colour=labels_col, fontface=labFontFace)+
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
  
  # ggd1And2 <- ggdraw() + # lower left corner coordinates
  #   draw_plot(ggd1, x=0, y=0, height=1, width=cladoWidth) +
  #   draw_plot(ggd2, x=ylabposition, y=0, height=1, width=labelWidth)
  # 
  # return(ggd1And2)
  
  
  ##########################################################################################################
  ##########################################################################################################  MATRIX HEATMAP CENTRAL PLOT
  ##########################################################################################################
  rowcol_order <- labels(row.dend_save)
  # my_ylab <- rowcol_order

  ## order of the dendros
  col.ord <- match(rowcol_order, colnames(x))
  row.ord <- match(rowcol_order, rownames(x))
  if(plotMap == "upper_right" | plotMap == "lower_left" | plotMap == "reverseSquare" | plotMap == "square") {
    xx <- x[rev(row.ord),col.ord]
  } else {
    xx <- x[row.ord,col.ord]  
  }

  if(plotMap != "square") {
    if(plotMap == "lower_left")  xx[lower_right_tri(xx)] <- NA
    if(plotMap == "lower_right")  xx[lower.tri(xx)] <- NA
    if(plotMap == "upper_left")  xx[upper.tri(xx)] <- NA
    if(plotMap == "upper_right")  xx[upper_left_tri(xx)] <- NA
  }
  xx_tmp <- xx
  xx_tmp <- melt(xx_tmp)
  colnames(xx_tmp)[1:2] <- c("X1", "X2")
    
  dimnames(xx) <- NULL
  xx <- melt(xx)
  stopifnot(ncol(xx) == 3)
  colnames(xx)[1:2] <- c("X1", "X2")
  xx$tmp_cond1 <- xx_tmp$X1 # !!! was xx_tmp$X1 and xx_tmp$X2 in the manuscript version, without renaming the columns ??? !!!
  xx$tmp_cond2 <- xx_tmp$X2
  

  if(is.null(low_limit_col)) low_limit_col <- min(xx$value, na.rm=T)
  if(is.null(high_limit_col)) high_limit_col <- max(xx$value, na.rm=T)
  
  stopifnot("value" %in% colnames(xx))
  xx$value <- round(xx$value,2)
  xx$value_format <- sprintf("%.2f", round(xx$value,2))
  
  botMarg <- 0
  rightMarg <- 2
  
  #### PREPARE THE Z-SCORE DATA AND THE LABLES FOR THE ENTRY OF THE MATRIX
  if(!is.null(zscoreDT)) {
    xx$comparison <- unlist(sapply(c(1:nrow(xx)), function(x) {
      paste0(sort(c(as.character(xx$tmp_cond1[x]), as.character(xx$tmp_cond2[x])))[1], "_", 
             sort(c(as.character(xx$tmp_cond1[x]), as.character(xx$tmp_cond2[x])))[2])
    }))
    xx$comparison <- as.character(xx$comparison)
    zscoreDT$comparison <- as.character(zscoreDT$comparison)
    
    xx <- left_join(xx, zscoreDT, by="comparison")
    xx$MoC <- round(xx$MoC, 2)
    xx$zscores <- round(xx$zscores, 2)   
    xx$zscores_format <- sprintf("%.2f", round(xx$zscores,2))
    stopifnot(all(na.omit(as.character(xx$MoC) == as.character(xx$value))))
    xx$mylabel <- unlist(sapply(c(1:length(xx$value)), function(x) ifelse(is.na(xx$value[x]), NA, paste0(xx$value_format[x], "\n(", xx$zscores_format[x], ")"))))
    xx$mylabel[xx$tmp_cond1 == xx$tmp_cond2] <- "1"
  } else {
    xx$mylabel <- xx$value_format
  }
  
  #### PREPARE MAIN PLOT CENTRAL HEATMAP
  if(annotateMean) {
    if(plotMap == "square") {
      xExpand <- 5.5
    } else {
      xExpand <- 2
    }
  }
  
  centre.plot <- ggplot(xx, aes(x = X2, y = X1)) + geom_tile(aes(fill=value), colour="white") +
    scale_x_continuous(breaks = 1:ncol(x),labels=rowcol_order, expand=c(0,0)) +
    scale_y_continuous(breaks = 1:ncol(x),labels=rev(rowcol_order), expand=c(0,0)) +
    scale_fill_gradient(low = low_values_col, high =high_values_col , na.value = na_values_col,
                        limit = c(low_limit_col,high_limit_col), 
                        space = "Lab",
                        name=fill_legName)+
    guides(fill = guide_colorbar(barwidth = 10, barheight = 1,title.position = "top", title.hjust = 0.5))+
				# top, right, bottom, and left
				# 7.11 changed left from 0 to -1
    theme(plot.margin = unit(c(0,rightMarg,botMarg,-1.5), "lines"),
          legend.justification = c(1, 0),
          legend.position = c(1, 0.9),
          # legend.position = c(1, 0.9),
          legend.title = element_text(size=legTitSize), 
          legend.text = element_text(size = legTextSize),
          legend.direction = "horizontal")
  
  if(annotateMat) 
    centre.plot <- centre.plot + geom_text(aes(X2, X1, label = mylabel), color = "black", size = cellSize)
  
  if(annotateMean){
    nCond <- length(unique(xx$tmp_cond1))
    stopifnot(is.numeric(xx$value[1]))
    # do not take the diago count.
    tmp_xx <- xx[xx$tmp_cond1 != xx$tmp_cond2,]
    mean_MoC <- unlist(sapply(rowcol_order, function(x) mean(tmp_xx$value[tmp_xx$tmp_cond1 == x | tmp_xx$tmp_cond2 == x], na.rm=T)))
	mean_MoC_format <- sprintf("%.2f", round(mean_MoC,2))
    mean_mean_MoC <- mean(mean_MoC, na.rm=T)
	mean_mean_MoC_format <- sprintf("%.2f", round(mean_mean_MoC,2))
	
	save(rowcol_order, file="rowcol_order.Rdata")

    if(!is.null(zscoreDT)) {
      stopifnot(is.numeric(xx$zscores[1]))
      mean_zscores_format <- sprintf("%.2f", round(mean_zscores,2))
      mean_zscores <- unlist(sapply(rowcol_order, function(x) mean(xx$zscores[xx$tmp_cond1 == x | xx$tmp_cond2 == x], na.rm=T)))
      meanRightLab <- paste0(as.character(mean_MoC_format), "\n(",as.character(mean_zscores_format), ")")
      meanRightLeg <- paste0(comparisonName, " mean ", fill_legName, " without diag (mean: ", mean_mean_MoC_format, ").\n(z-scores)")
    } else {
      meanRightLab <- paste0(as.character(mean_MoC_format))
      # meanRightLeg <- paste0(comparisonName, " mean MoC without diag (mean: ", round(mean_mean_MoC,2), ").")
      meanRightLeg <- paste0(comparisonName, " mean ", fill_legName, " without diag\n(mean: ", mean_mean_MoC_format, ")")
    }
	
	  names(meanRightLab) <- rowcol_order
	
    # add some space for the last label
    centre.plot <- centre.plot + expand_limits(x = nCond + xExpand) 
    
    if(plotMap == "square") {
      rowAnnotPos_x <- rep(nCond + 1 , nCond)
    }else if(plotMap == "lower_left") {
      rowAnnotPos_x <- c(1:nCond) + 1 
    } else{
      stop("not implemented\n")
    }
    
    save(meanRightLab, file="meanRightLab.Rdata")
    if(!is.null(perso_rightLab)) {
      save(meanRightLab, file ="meanRightLab.Rdata")
      save(perso_rightLab, file ="perso_rightLab.Rdata")
      stopifnot(setequal(names(meanRightLab), names(perso_rightLab)))
      perso_rightLab <- perso_rightLab[names(meanRightLab)]
      meanRightLab <- perso_rightLab 
    } 
    save(meanRightLeg, file="meanRightLeg.Rdata")
    if(!is.null(perso_rightLeg)) meanRightLeg <- perso_rightLeg
    
    centre.plot <- centre.plot + annotate("text", 
                                            size = cellSize,
                                           label = meanRightLab,
                                           fontface = 'italic',
                                           x = rowAnnotPos_x, y = rev(c(1:nCond)))
    # add the "legend" for the annotation
    rowLegPosX <- ifelse(plotMap == "square", nCond+3.7, nCond )
    centre.plot <- centre.plot + annotate("text", 
                            size = cellSize,
                            label = meanRightLeg,
                            fontface = 'italic',
                            # hjust = 1,
                            # x = nCond+2, 
                            hjust = 0.5,
                            x = rowLegPosX,
                            y = nCond-2)
  }
  
  ### IF PROVIDED CHANGE THE COLOR OF THE LABELS
  curr_xlab <- rowcol_order
  if(!is.null(lab_color_vect)) {
    curr_xlab_col <- unlist(sapply(as.character(curr_xlab), function(x) {
      if(x %in% names(lab_color_vect)) return(as.character(lab_color_vect[x]))
      return("black")
    }))
  } else {
    curr_xlab_col <- curr_ylab_col <- "black"
  }

  
  if(!is.null(legCategoryCols)){
    labelDot <- rep(" \u25cf", length(legCategoryCols))
    centre.plot <- centre.plot + annotate("text",
                                          size=12,
                                          label = labelDot,
                                          color = legCategoryCols,
                                          hjust = rep(0, length(legCategoryCols)),
                                          vjust = rep(0.5, length(legCategoryCols)),
                                          x = rowLegPosX-1,
                                          y=seq(nCond-5 , by=-1, length.out = length(legCategoryCols)))
                                          # y = c(nCond-5, nCond-6, nCond-7, nCond-8))
    labelDotTxt <- names(legCategoryCols)
    labelDotTxt[labelDotTxt == "statModels"] <- "stat. model"
    labelDotTxt[labelDotTxt == "oneDimScores"] <- "1D score"
    centre.plot <- centre.plot + annotate("text",
                                          size = cellSize,
                                          label = labelDotTxt,
                                          color = legCategoryCols,
                                          hjust = rep(0, length(legCategoryCols)),
                                          vjust = rep(1, length(legCategoryCols)),
                                          x = rowLegPosX-0.2,
                                          seq(nCond-5 , by=-1, length.out = length(legCategoryCols)))
                                          # y = c(nCond-5, nCond-6, nCond-7, nCond-8))
  }

  
  #### ADD XAXIS LABELS AS A SINGLE PLOT
  xaxis_DT <- data.frame(x = 1:max(xx$X1), y = rep(0, max(xx$X1)), value = rowcol_order)
  xaxis.plot <- ggplot(xaxis_DT, aes(x=x, y=y))+ geom_tile(fill="white", colour="white")+
    scale_x_continuous(breaks = 1:ncol(x),labels=rowcol_order, expand=c(0,0)) +
    scale_y_continuous(breaks = 1:ncol(x),labels=rev(rowcol_order), expand=c(0,0)) +
    geom_text(aes(x, y+0.4, label=value),size=axisLabelSize, hjust=1, angle=90, colour=curr_xlab_col, fontface=labFontFace)+
    # geom_text(aes(x, y+0.4, label=value),size=4, hjust=1, angle=45, colour=curr_xlab_col)+  
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0,rightMarg,botMarg,-1.5), "lines")
          # plot.margin = unit(c(0,rightMarg,0,0), "lines")
          )
  if(annotateMean){
    xaxis.plot <- xaxis.plot + expand_limits(x = nCond + xExpand) 
  }

  centre.plot <- centre.plot + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  
  gg_y_dendrolab <- gg_y_dendrolab + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  
  
  
  # if(plotMap == "square"){
  #   centre_legend <- g_legend(centre.plot)
  #   centre.plot <-  centre.plot + theme(legend.position = "none")
  # }
  ##########################################################################################################
  ##########################################################################################################  ASSEMBLE ALL THE PLOTS
  ##########################################################################################################
  # offsetXaxisPlot <- 0.1
  # widthDendro <- 0.1
  # widthLabels <- 0.5
  
  ### ASSEMBLE THE LEFT PART
  p_left <- ggdraw() +
         # draw the row dendrogram
         draw_plot(gg_row_dendro, x=0, y=0, width=0.18,height=1) +  # was 0.02 initially
         # draw the yaxis label
         draw_plot(gg_y_dendrolab, x=-0.01, y=0, width=0.5,height=1)  # -0.05 if left aligned
  
  #### ASSEMBLE THE CENTRAL PART with the left part
  p_centre <- ggdraw()+
         # draw the dendro
         draw_plot(p_left, x=0, y=0.0, width=0.7,height=1) +
         # draw the heatmap
         # draw_plot(centre.plot, x=0.22, y=0, width=0.78,height=1) 
        draw_plot(centre.plot, x=0.21, y=0, width=0.8,height=1) 
    
  # p_bottom <- ggdraw()+draw_plot(xaxis.plot, x=0.22, y=0, width=0.78,height=1)
  p_bottom <- ggdraw()+draw_plot(xaxis.plot, x=0.21, y=0, width=0.8,height=1)
  
  p_final <- ggdraw() + 
    # draw dendro + labels
    draw_plot(p_centre, x=0, y=0.12, width=1,height=0.88) +
    # draw heatmap + x axis
    draw_plot(p_bottom, x=0, y=0, width=1,height=0.12)
  
  return(p_final)
  
}
# 
