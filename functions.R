## Helper functions for computing and plotting posterior predictive p-values
kappalist <- function(v) {
    ## Creates a list of kappa values from a vector of BUGS output
    p <- max(as.numeric(gsub("(kappa\\[|,[[:digit:]]+\\])", "", grep("kappa", names(v), value=TRUE))))
    lapply(1:p, function(j) v[grep(paste("kappa\\[", j, ",[[:digit:]]+\\]", sep=""), names(v))])
}
newdata <- function(theta, alpha, kappa){
    ## Generates a new data set from a graded response model
    ## theta = vector (length n) of "ability" values
    ## alpha = vector (length p) of discrimination parameters
    ##         or a single value (vector of length 1)
    ## kappa = a list of kappa values
    n <- length(theta)
    p <- length(kappa)
    Ystar <- t(sapply(theta, function(th) rlogis(p, alpha*th, 1.0)))
    Y <- t(apply(Ystar, 1, function(y) sapply(seq(p), function(j) cut(y[j], breaks=c(-Inf, kappa[[j]], Inf), labels=FALSE))))
    dimnames(Y) <- NULL
    return(Y)
}
ppktau <- function(sims){
    require(plyr)
    parnames <- colnames(sims)
    dumfun <- function(v){
        theta <- v[grep("theta", parnames)]
        alpha <- v[grep("alpha", parnames)]
        kappa <- kappalist(v)
        p <- length(kappa)
        out <- cor(newdata(theta, alpha, kappa), method="pearson")[lower.tri(diag(p))]
        return(out)
    }
    out <- adply(sims, 1, dumfun, .progress="text")
    return(out[,-1])
}
densityShade <- function(x, from, to, col=rgb(0.5, 0.5, 0.5, 0.5), border=NA, ...){
    den <- density(x)
    from <- max(from, min(den$x))
    to <- min(to, max(den$x))
    den <- density(x, from=from, to=to)
    x <- den$x
    y <- den$y
    xx <- c(x, rev(x))
    yy <- c(y, rep(0, times=length(x)))
    polygon(xx, yy, col=col, border=border, ...)
}
makeNames <- function(name, ..., symmetric.matrix=FALSE, diag=TRUE){
    if (symmetric.matrix){
        idx <- which(upper.tri(diag(max(...)), diag=diag), arr.ind=TRUE)
    } else {
        idx <- expand.grid(rev(list(...)))
    }
    name <- paste(name, "[", sep="")
    if (ncol(idx)>1){
        for (j in ncol(idx):2){
            name <- paste(name, idx[,j], ",", sep="")
        }
    }
    name <- paste(name, idx[,1], "]", sep="")
    return(name)
}
plotpppval <- function(sims, true, xlim=NULL){
    if (!is.matrix(sims))
        sims <- as.matrix(sims)
    npar <- ncol(sims)
    ##mult.fig(npar, mar=c(3, 1.0, 1.0, 1.0), oma=rep(0, 4))
    mult.fig(npar, mar=c(2.5, 0.25, 0.25, 0.25), oma=rep(0, 4))
    den <- apply(sims, 2, density)
    if (is.null(xlim)){
        xl <- range(as.vector(sapply(den, function(d) range(d$x))))
        xl <- xl + c(-1, 1)*0.04*diff(xl)
    }
    else xl <- xlim
    ylmax <- max(sapply(den, function(d) max(d$y)))
    yl <- c(0, 1.04*ylmax)
    for(i in 1:npar){
        plot(density(sims[,i]), xlim=xl, ylim=yl, tcl=-0.2, mgp=c(1, 0.1, 0.0), xlab=colnames(sims)[i], ylab="", bty="n", yaxt="n", main="")
        if (mean(true[i]<sims[,i])>0.5){
            from <- -Inf
            to <- true[i]
            pval <- formatC(mean(true[i]>sims[,i]), digits=2, format="f")
        }
        else{
            from <- true[i]
            to <- Inf
            pval <- formatC(mean(true[i]<sims[,i]), digits=2, format="f")
        }
        densityShade(sims[,i],from=from, to=to)
        text(xl[1], 0.75*ylmax, labels=pval, pos=4)
        ##abline(v=true[i], col="red")
    }
}

