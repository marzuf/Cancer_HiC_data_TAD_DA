
printVar <- function(x){
  cat(paste0(x, " = ", eval(parse(text=x)), "\n"))
}


plot_cumMatch <- function(dt, tomatch){
  curr_match <- na.omit(dt[, tomatch])
  xvect <- seq_len(max(curr_match))
  yvect <- sapply(xvect, function(x){
    sum(curr_match >= x)
  })
  plot(x = xvect,
       y = yvect,
       xlab = paste0("# datasets in which matching signif. TAD"), 
       ylab = paste0("# query TAD"),
       type="l")
  mtext(side=3, text=paste0(tomatch))
}


ggplot_barplot_hicdsexprds <- function(barDT, xvar, yvar,
                                       xcolvar = NULL,
                                       myxlab="", myylab="", myTit="", mySub="", 
                                       barCol="dodgerblue3") {
  if(is.null(xcolvar)) {
    plotcols <- "black"
  }else {
    plotcols <- barDT[,xcolvar]
  }
  
  p_ref <- ggplot(barDT, aes_string(x = xvar, y = yvar)) +
    ggtitle(myTit, subtitle = mySub)+
    geom_bar(stat="identity", position = "dodge", fill = barCol)+
    scale_x_discrete(name=myxlab)+
    scale_y_continuous(name=myylab,
                       breaks = scales::pretty_breaks(n = 10))+
    labs(fill = "")+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=8, angle = 90, color = plotcols),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      # axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
    )
  return(p_ref)
}  

ggplot_boxplot_hicdsexprds <- function(barDT, xvar, yvar, colvar,
                                       myxlab="", myylab="", myTit="", mySub="", 
                                       barCol="dodgerblue3") {
  avg_barDT <- aggregate(as.formula(paste0(yvar, "~", xvar)), data = barDT, FUN=mean, na.rm=TRUE)
  xvar_order <- as.character(avg_barDT[,xvar][order(avg_barDT[,yvar], decreasing=TRUE)])
  stopifnot(!is.na(xvar_order))
  barDT[, xvar] <- factor(as.character(barDT[,xvar]), levels = xvar_order)
  plotDT <- na.omit(barDT)  
  
  if(is.null(colvar)){
    mycols <- "black"
  }else {
    mycols <- plotDT[,colvar]
  }
  
  p_ref <- ggplot(plotDT, aes_string(x = xvar, y = yvar)) +
    ggtitle(myTit, subtitle = mySub)+
    geom_boxplot(fill = barCol)+
    scale_x_discrete(name=myxlab)+
    scale_y_continuous(name=myylab,
                       breaks = scales::pretty_breaks(n = 10))+
    labs(fill = "")+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=8, angle = 90, color = mycols),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      # axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
    )
  return(p_ref)
}  
