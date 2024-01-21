ggtally <- function (G, p1, adj = 0, values = seq(-1, 1, by = 0.2), dotsize = 0.1, 
    dotcolour = "black") 
{
    PA1 <- PA2 <- NULL
    if (is.vector(G)) {
        G <- matrix(G, ncol = 2)
        Dg <- matrix(G %*% t(G), 1, 1)
    }
    else {
        G <- G[, 1:2]
        Dg <- diag(diag(G %*% t(G)))
    }
    nr <- nrow(G)
    for (i in 1:nr) {
        Z <- NULL
        for (j in 1:length(values)) {
            DP <- (values[j] + adj) * G[i, 1:2]/Dg[i, i]
            Z <- rbind(Z, DP)
        }
        Zdf <- as.data.frame(Z)
        colnames(Zdf) <- c("PA1", "PA2")
        p1 <- p1 + geom_point(data = Zdf, aes(x = PA1, y = PA2), 
            size = dotsize, colour = dotcolour)
        Zt <- Z[values >= 0, ]
        Zt <- data.frame(Zt)
        colnames(Zt) <- c("PA1", "PA2")
        p1 <- p1 + geom_path(data = Zt, aes(x = PA1, y = PA2), colour = "blue", linewidth = 0.25)
        Zt <- Z[values <= 0, ]
        Zt <- data.frame(Zt)
        colnames(Zt) <- c("PA1", "PA2")
        p1 <- p1 + geom_path(data = Zt, aes(x = PA1, y = PA2), colour = "red", linewidth = 0.25)
    }
    return(p1)
}
