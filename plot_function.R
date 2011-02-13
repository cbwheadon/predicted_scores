plotCondSumScoreDist<-function(itm.scres,cssd,nm){

	#Compare up to 3 models and observed score distributions
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
	
	d <- density(scres)
	ylim <- max(d$y)*1.05
	
	plot(x=c(0,ttl+1),y=c(0,ylim), type="n",ylab="density",xlab="score")
		
	cols=c(rgb(0,0,255,20,maxColorValue=255),rgb(136,255,0,20,maxColorValue=255),rgb(0,170,0,20,maxColorValue=255))
	lcols=c(rgb(0,0,255,100,maxColorValue=255),rgb(136,255,0,100,maxColorValue=255),rgb(0,170,0,100,maxColorValue=255))
	
	#models
	for (i in 1:length(cssd)){
		for(j in n.sims:1){
			mdl <- round(cssd[[i]][,j],0)
			scres <- unlist(lapply(1:ttl+1, function(p) rep(p-1,mdl[p])))
			lines(density(scres), col = cols[i], lwd = 1,lty=4)	
		}
	}
	#observed
	scres <- rowSums(itm.scres)
	lines(density(scres), col = 1, lwd = 2)
	
	legend(x=0,y=ylim, legend=nm, lty=1,col=lcols[1:i])
}

plotQQ <- function(itm.scres,cssd,nm){
	#q-q plot according to median values
	ttl <- length(cssd[[1]][,1])-1
	#observed candidates
	cs <- length(itm.scres[,1])
	#observed scores
	scres <- rowSums(itm.scres)
	
	plot(x=c(0,ttl+1),y=c(0,ttl+1),xlim=c(0,ttl+1),ylim=c(0,ttl+1),xlab="Raw scores",ylab="Model scores",type="n")
	cols=c(rgb(0,0,255,50,maxColorValue=255),rgb(68,153,0,50,maxColorValue=255),rgb(0,170,0,50,maxColorValue=255))
	lcols=c(rgb(0,0,255,100,maxColorValue=255),rgb(68,153,0,100,maxColorValue=255),rgb(0,170,0,100,maxColorValue=255))
		
	for (i in 1:length(cssd)){
		pb <- runif(cs)
		mdn <- apply(cssd[[i]],1,median)/cs
		thresh <- c(cumsum(mdn),1)
		mdl.scres <- sort(laply(pb, function (j) min(which(thresh>=j))-1))
		points(sort(scres),mdl.scres,col=cols[i],pch=i)
	}
	abline(a=0,b=1,lty=4)
	legend(x=0,y=ttl, nm, pch=c(1:i),col=lcols[1:i])
}
