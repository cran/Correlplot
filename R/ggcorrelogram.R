ggcorrelogram <- function (R, labs = colnames(R), ifun = "cos", 
                           cex = 1, main = "", ntrials = 50, 
                           xlim = c(-1.2, 1.2), 
                           ylim = c(-1.2, 1.2), hjust = 1,
                           vjust = 2, size = 2, main.size = 8) 
{
    theta <- fit_angles(R, ifun = "cos", ntrials = ntrials)
#    X <- matrix(c(cos(theta), sin(theta)), ncol = 2)
    PA1 <- cos(theta)
    PA2 <- sin(theta)
    if (is.null(labs)) 
        labs <- paste("X",1:ncol(R),sep="")
    X <- data.frame(PA1,PA2,labs)
    Graph <- ggplot(X, aes(x = PA1, y = PA2)) + 
            geom_text(aes(label = labs),
                      size = size, alpha = 1,
                      vjust = vjust, hjust = hjust) +
            coord_equal(xlim = xlim, ylim = ylim) +
            xlab("") + ylab("") +  ggtitle(main) +
            theme(plot.title = element_text(hjust = 0.5, 
                                            size = main.size),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank())
    xc <- 0; yc <- 0; r <- 1
    Graph <- Graph + annotate("path",
     x=xc+r*cos(seq(0,2*pi,length.out=100)),
     y=yc+r*sin(seq(0,2*pi,length.out=100)),linewidth=0.25) 
    
    Graph <- Graph + geom_segment(aes(x=0, y=0, xend=PA1, yend=PA2),
                         arrow=arrow(length = unit(0.2,"cm"), 
                                     angle = 10), alpha=0.75,
                         linewidth = 0.25,
                         color="blue") + ggtitle(main)
    Graph$theta <- theta
    return(Graph)
}
