# PPMC Checks

ppmc_checks <- function(sim.data){


	x11(width=8,height=8)

	#Inter Item Correlations (Density)
	require("sfsmisc")
	p<-length(sim.data[[1]][1,])
	dumfun <- function(x){
		p <- length(x[1,])
		return(cor(x, method="pearson")[lower.tri(diag(p))])
	}

	pp.eq <- ldply(sim.data, dumfun, .progress="text")
	colnames(pp.eq) <- makeNames("r", 1:p, 1:p, symmetric.matrix=TRUE, diag=FALSE)
	ktau <- cor(Y, method="pearson")[lower.tri(diag(p))]
	plotpppval(pp.eq, ktau, xlim=c(0.0, 1.0))

	#Inter Item Correlations PPP	
	
	rs <- length(pp.eq[1,])

	ppp <- laply(1:rs, function(y) sum(pp.eq[,y]>ktau[y])/length(sim.data))
	ppp <- cut(ppp, breaks=c(-Inf,.05,.95,Inf), labels=c("<0.05","0.05<ppp<0.95",">0.95"))

	q1 <- as.integer(gsub("(r\\[|,[[:digit:]]+\\])", "",names(pp.eq)))
	q2 <- as.integer(gsub("r\\[[[:digit:]]+,|\\]","",names(pp.eq)))

	ppp.values <- data.frame(q1=q1,q2=q2,ppp=ppp)

	
	
	x11(width=4,height=4)
	p <- ggplot(ppp.values, aes(x=q1, y=q2,shape=ppp,colour=ppp))
	p <- p + layer(geom = "point")
	#p <- p + facet_grid(. ~ section) 
	p <- p + xlab("item")
	p <- p + ylab("item")
	p <- p + geom_point(size=4)
	p <- p + scale_colour_manual(values = c("darkblue", "grey","orange"))
	print(p)

	#Item Test Correlations
	x11()
	item.test.r <- laply(sim.data,  function(s1) aaply(s1,2, function(s) cor(s,rowSums(s1))))

	obs.r <- aaply(Y,2, function(s) cor(s,rowSums(Y)))
	obs.r.df <- cbind(1:length(Y[1,]),as.numeric(obs.r))
	obs <- data.frame(Sim=length(sim.data)+1,Item=obs.r.df[,1],R=obs.r.df[,2])
	
	df <- data.frame(melt(item.test.r))

	names(df) <- c("Sim","Item","R")
	p <- ggplot(df, aes(x=Item,y=R,colour="simulated"))
	p <- p + geom_point(position = "jitter")
	p <- p + geom_point(data=obs,aes(x=Item,y=R,colour="observed"))
	p <- p + scale_colour_manual(name="type",values = c("grey","red"))
	p <- p + scale_y_continuous(name="Item Test Correlation")
	print(p)
	
}

