#Generalized Partial Credit Model

my.gpcm <- function(pars,abi,n){
	alpha <- as.numeric(pars[1])
	beta <- as.numeric(pars[-1])
	denor <- sum(exp(cumsum(alpha*(abi-beta))))
	ps <- exp(cumsum(alpha*(abi-beta)))/(1+denor)
	z.p <- 1 - sum(ps)	
	ps <- c(z.p,ps,rep(0,n-length(ps)-1))
	return(ps)
}

parlist <- function(v) {
	p <- max(as.numeric(gsub("(beta\\[|,[[:digit:]]+\\])","",grep("beta", names(v), value=TRUE))))
	alpha <-  v[grep("alpha", names(v))]
	if(length(alpha)==0){
		#Partial Credit Model
		alpha<- rep(1,p)
	}
	lapply(1:p, function(j) c(alpha[j],v[grep(paste("beta\\[", j, ",[[:digit:]]+\\]", sep=""), names(v))]))
}

cond.probs <- function(v){
	theta <- as.numeric(v[grep("theta", names(v))])
	pars <- parlist(v)
	K <- as.matrix(table(as.numeric(gsub("(beta\\[|,[[:digit:]]+\\])","",grep("beta", names(v), value=TRUE)))))
	K <- K + 1
	cts <- max(K)
	ps <- llply(theta, function(th) laply(pars,my.gpcm,th,cts),.progress="none")
	f.ps <- laply(ps, function (pss) wLord(as.matrix(pss),K),.progress="none")
	sum.cond.sum.sc.d <- colSums(f.ps)
	return(sum.cond.sum.sc.d)
}

f.probs <- function(v,cts){
	theta <- as.numeric(v[grep("theta", names(v))])
	pars <- parlist(v)
	ps <- mlply(theta, function(th) laply(pars,my.gpcm,th,cts),.progress="none")
}

sum.prob <- function(sims,cts,cats,max.scr){ 
	
	n.sims <- length(sims[,1])
	scssd <- matrix(0,nrow=max.scr+1,ncol=n.sims)
	pb <- txtProgressBar(min = 1, max = n.sims, style = 3)
	for(i in 1:n.sims){
		v <- sims[i,]
		theta <- as.numeric(v[grep("theta", names(v))])
		pars <- parlist(v)
		cssd <- matrix(0,nrow=max.scr+1,ncol=n)
	
		for(j in 1:n){
	
			#For each theta calculate probability distribution
			ps <- numeric(cts*p)
			resp <- numeric(p)

			for(k in 1:p){
				#For each question calculate probability
				strt <- cts*(k-1)+1
				ed <- strt+cts-1 
				ps[strt:ed]<-my.gpcm(as.numeric(pars[[k]]),theta[j],cts)
			}

			#Conditional Probability Distribution
			ps.t <- matrix(ps,nrow=p,ncol=cts,byrow=TRUE)
			cssd[,j] <- wLord(ps.t,cats)
		}

		scssd[,i] <- rowSums(cssd)
		setTxtProgressBar(pb, i)
	}
	return(scssd)
}