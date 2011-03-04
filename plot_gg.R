plotCondSumScoreDist<-function(itm.scres,cssd,nm){
	
	#Compare models and observed score distributions
	require(ggplot2)

	x11()

	if(typeof(cssd)=="double"){
		cssd <- list(cssd)
		nm <- list(nm)
	}
	
	#total score
	ttl <- length(cssd[[1]][,1])-1
	#observed candidates
	cs <- length(itm.scres[,1])
	#sim
	n.sims <- length(cssd[[1]][1,])
	mdls <- matrix(ncol=4)
	mdls<-mdls[-1,]
	#observed
	scres <- rowSums(Y)
	obs <- matrix(cbind(0:ttl,table(factor(scres, levels = 0:ttl))),ncol=2)
	
	#models
	for (i in 1:length(cssd)){
		mdn <- apply(cssd[[i]],1,median)
		m05 <- apply(cssd[[i]],1,quantile,c(.05))
		m95 <- apply(cssd[[i]],1,quantile,c(.95))
		mdls<-rbind(mdls,cbind(i,1,0:ttl,mdn))
		mdls<-rbind(mdls,cbind(i,2,0:ttl,m05))
		mdls<-rbind(mdls,cbind(i,3,0:ttl,m95))
		mdls<-rbind(mdls,cbind(i,4,obs))
	}
	
	mdls <- data.frame(mdls)
	names(mdls) <- c("model","type","score","frequency")
	mdls$model <- factor(mdls$model, labels = nm)
	mdls$type <- factor(mdls$type, labels = c("median","05","95","observed"))
	
	p <- ggplot(mdls,aes(x=score,y=frequency,colour=type,linetype=type))
	p <- p + geom_line()
	p <- p + facet_wrap(~ model)
	p <- p + scale_colour_manual(values = c(rep("black",3),"blue")) 
	p <- p + scale_linetype_manual(values=c(rep(2,3),1)) 
	
	#ggsave(file=fn)
	print(p)

}

plotQQ <- function(itm.scres,cssd,nm){

	#Q-Q Plot to compare models and observed score distributions
	
	require(plyr)	
	require(ggplot2)

	x11() 

	if(typeof(cssd)=="double"){
		cssd <- list(cssd)
		nm <- list(nm)
	}
	ttl <- length(cssd[[1]][,1])-1
	#observed candidates
	cs <- length(itm.scres[,1])
	#observed scores
	scres <- sort(rowSums(itm.scres))
	mdls <- matrix(ncol=3)
	mdls<-mdls[-1,]
	for (i in 1:length(cssd)){
		pb <- runif(cs)
		mdn <- apply(cssd[[i]],1,median)/cs
		thresh <- c(cumsum(mdn),1)
		mdl.scres <- sort(laply(pb, function (j) min(which(thresh>=j))-1))
		mdls<-rbind(mdls,(cbind(i,scres,mdl.scres)))
	}
	
	mdls <- data.frame(mdls)
	names(mdls) <- c("model","observed","modelled")
	mdls$model <- factor(mdls$model,labels=nm)
	
	p <- ggplot(mdls,aes(x=observed,y=modelled))
	al <-  100/cs
	if(al>1){al<-1}
	p <- p + geom_jitter(alpha = al) + geom_abline(intercept = 0, slope = 1)
	p <- p + facet_wrap(~ model)
	
	#ggsave(file=fn)
	print(p)
}


