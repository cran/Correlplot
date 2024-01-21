ggbplot <- function(A,B,main="",circle=TRUE,
                    xlab="",ylab="", main.size=8,
                    xlim=c(-1,1),ylim=c(-1,1),
                    rowcolor = "red", rowch = 1, 
                    colcolor = "blue", colch = 1, 
                    rowarrow=FALSE, colarrow = TRUE) {
  PA1 <- PA2 <- strings <- NULL		    
  a <- ggplot(A) + 
            coord_equal(xlim = xlim, ylim = ylim) +
            ggtitle(main) +
            xlab(xlab) + ylab(ylab) +
            theme(plot.title = element_text(hjust = 0.5, 
                                            size = main.size),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank())
  a <- a + geom_point(data = A, aes(x = PA1, y = PA2), 
                      colour = rowcolor, 
                      shape = rowch)
  a <- a + geom_text(aes(x = PA1, y = PA2, label = strings), 
                       size = 3, alpha = 1, vjust = 2)
  a <- a + geom_point(data = B, aes(x = PA1, y = PA2), 
                           colour = colcolor, 
                           shape = colch)
  a <- a + geom_text(data = B, aes(x = PA1, y = PA2, 
                                   label = strings), 
                        size = 3, alpha = 1, vjust = 2) 
  if(rowarrow) {
    a <- a + geom_segment(data = A, aes(x=0, y=0, xend=PA1, yend=PA2),
                         arrow=arrow(length = unit(0.2,"cm"), 
                                     angle = 10), 
                         alpha=0.75,
                         linewidth = 0.25,
                         color= rowcolor) 
  } 
  if(colarrow) {
        a <- a + geom_segment(data = B, 
                           aes(x = 0, y = 0, xend = PA1, 
                               yend = PA2),
                         arrow=arrow(length = unit(0.2,"cm"), 
                                     angle = 10), 
                         alpha=0.75,
                         color=colcolor)   
  } 
  if(circle) {
     xc <- 0
     yc <- 0
     r <- 1
     a <- a + annotate("path",
     x=xc+r*cos(seq(0,2*pi,length.out=100)),
     y=yc+r*sin(seq(0,2*pi,length.out=100)),linewidth=0.25) 
   }
   return(a)
}
