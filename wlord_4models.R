## Comparison of 4 dichotomous models
rm(list=ls())

#Make sure this is correct so the bugs syntax file is found
setwd(file.path(getwd(), "/predicted_scores"))

library("R2WinBUGS")
library("mcmcplots")
library("inline")
source("wl_cpp.R")
source("plot_gg.R")
source("pfunctions.R")

#Change this to your bugs directory
bugs.directory = "C:/Users/User/Winbugs/WinBUGS14"

Ypath <- file.path(getwd(), "data/chem.csv")
Y <- as.matrix(read.csv(file=Ypath,header=FALSE,nrows=1000))

n <- nrow(Y)
p <- ncol(Y)

m.delta <- 0.0
s.delta <- 1.0

iter <- 5000
burnin <- 4000
thin <- 10

data <- list("Y", "n", "p",
		"m.delta", "s.delta")
monitor <- c("delta", "theta")
bugs.file <- file.path(getwd(), "bugs/rasch.bug")

system.time(winbugsout <- bugs(data=data, inits=NULL, parameters.to.save=monitor,
				model.file=bugs.file,
				n.iter=iter, n.thin=thin, n.burnin=burnin,bugs.directory=bugs.directory))

sims <- winbugsout$sims.matrix
n.sims <- length(sims[,1])
parnames <- colnames(sims)

cats <- matrix(rep(2,p),ncol=1)

scssd <- matrix(0,nrow=p+1,ncol=n.sims)

pb <- txtProgressBar(min = 1, max = n.sims, style = 3)

for (i in 1:n.sims){
	v <- sims[i,]
	theta <- v[grep("theta", parnames)]
	delta <- v[grep("delta", parnames)]
	
	probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) exp(th-delta[j])/(1+exp(th-delta[j]))))) 
	
	for (s in 1:n){
		ps <- matrix(cbind(1-probs[s,],probs[s,]),ncol=2)
		cond.sum.sc.d <- wLord(ps,cats)
		if(s==1){
			scssd[,i] <- cond.sum.sc.d
		} else {
			scssd[,i] <- scssd[,i] + cond.sum.sc.d
		}
	}
	setTxtProgressBar(pb, i)
}

close(pb)

rasch.scssd <- scssd

m.alpha <- 1.0
s.alpha <- 1.0

iter <- 5000
burnin <- 4000
thin <- 10

data <- list("Y", "n", "p",
		"m.alpha", "s.alpha",
		"m.delta", "s.delta")
monitor <- c("delta", "theta", "alpha")
bugs.file <- file.path(getwd(), "bugs/tpl.bug")

system.time(winbugsout <- bugs(data=data, inits=NULL, parameters.to.save=monitor,
				model.file=bugs.file,
				n.iter=iter, n.thin=thin, n.burnin=burnin,bugs.directory=bugs.directory))

sims <- winbugsout$sims.matrix
n.sims <- length(sims[,1])
parnames <- colnames(sims)

scssd <- matrix(0,nrow=p+1,ncol=n.sims)

pb <- txtProgressBar(min = 1, max = n.sims, style = 3)

for (i in 1:n.sims){
	v <- sims[i,]
	theta <- v[grep("theta", parnames)]
	delta <- v[grep("delta", parnames)]
	alpha <- v[grep("alpha", parnames)]
	
	probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) exp(alpha[j]*(th-delta[j]))/(1+exp(alpha[j]*(th-delta[j])))))) 
	
	for (s in 1:n){
		ps <- matrix(cbind(1-probs[s,],probs[s,]),ncol=2)
		cond.sum.sc.d <- wLord(ps,cats)
		if(s==1){
			scssd[,i] <- cond.sum.sc.d
		} else {
			scssd[,i] <- scssd[,i] + cond.sum.sc.d
		}
	}
	setTxtProgressBar(pb, i)
}

close(pb)
tpl.scssd <- scssd

a.eta <- 1.0
b.eta <- 1.0

guess.ind <- rep(1, p)

iter <- 10000
burnin <- 9000
thin <- 10

data <- list("Y", "n", "p", "guess.ind",
		"m.alpha", "s.alpha",
		"m.delta", "s.delta",
		"a.eta", "b.eta")

monitor <- c("alpha", "delta", "theta", "eta")
bugs.file <- file.path(getwd(), "bugs/thpl.bug")

system.time(winbugsout <- bugs(data=data, inits=NULL, parameters.to.save=monitor,
				model.file=bugs.file,
				n.iter=iter, n.thin=thin, n.burnin=burnin,bugs.directory=bugs.directory))

sims <- winbugsout$sims.matrix
n.sims <- length(sims[,1])
parnames <- colnames(sims)

cats <- matrix(rep(2,p),ncol=1)

scssd <- matrix(0,nrow=p+1,ncol=n.sims)

pb <- txtProgressBar(min = 1, max = n.sims, style = 3)

for (i in 1:n.sims){
	v <- sims[i,]
	theta <- v[grep("theta", parnames)]
	delta <- v[grep("delta", parnames)]
	alpha <- v[grep("alpha", parnames)]
	eta <- v[grep("^eta", parnames)]
	
	probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) eta[j]+((1-eta[j]) * (exp(alpha[j]*(th-delta[j]))/(1+exp(alpha[j]*(th-delta[j])))))))) 
	
	for (s in 1:n){
		ps <- matrix(cbind(1-probs[s,],probs[s,]),ncol=2)
		cond.sum.sc.d <- wLord(ps,cats)
		if(s==1){
			scssd[,i] <- cond.sum.sc.d
		} else {
			scssd[,i] <- scssd[,i] + cond.sum.sc.d
		}
	}
	setTxtProgressBar(pb, i)
}

close(pb)

thpl.scssd <- scssd

a.sigsq.gamma <- 1.0
b.sigsq.gamma <- 1.0

#Change this to represent testlet structure
d <- c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4))

#Number of testlets
n.t <- max(d) 

iter <- 5000
burnin <- 4000
thin <- 10

data <- list("Y", "n", "p", "d", "n.t",
		"m.alpha", "s.alpha",
		"m.delta", "s.delta",
		"a.sigsq.gamma", "b.sigsq.gamma")
monitor <- c("alpha", "delta", "theta", "gamma", "sigsq.gamma")
bugs.file <- file.path(getwd(), "bugs/tlet.bug")

system.time(winbugsout <- bugs(data=data, inits=NULL, parameters.to.save=monitor,
				model.file=bugs.file,
				n.iter=iter, n.thin=thin, n.burnin=burnin,bugs.directory=bugs.directory))

sims <- winbugsout$sims.matrix
n.sims <- length(sims[,1])
parnames <- colnames(sims)

cats <- matrix(rep(2,p),ncol=1)

scssd <- matrix(0,nrow=p+1,ncol=n.sims)

pb <- txtProgressBar(min = 1, max = n.sims, style = 3)

for (i in 1:n.sims){
	v <- sims[i,]
	theta <- v[grep("theta", parnames)]
	delta <- v[grep("delta", parnames)]
	alpha <- v[grep("alpha", parnames)]
	gamma <- v[grep("^gamma", parnames)]
	
	probs <- matrix(0,ncol=p,nrow=n)
	for (j in 1:n){
		for (k in 1:p){
			tstlt <- n.t * (j - 1) + d[k]
			probs[j,k] <- exp(alpha[k]*(theta[j]-delta[k]-gamma[tstlt]))/(1+exp(alpha[k]*(theta[j]-delta[k]-gamma[tstlt])))
		}
	} 
	
	for (s in 1:n){
		ps <- matrix(cbind(1-probs[s,],probs[s,]),ncol=2)
		cond.sum.sc.d <- wLord(ps,cats)
		if(s==1){
			scssd[,i] <- cond.sum.sc.d
		} else {
			scssd[,i] <- scssd[,i] + cond.sum.sc.d
		}
	}
	setTxtProgressBar(pb, i)
}

close(pb)

tplt.scssd <- scssd

plotCondSumScoreDist(Y,list(rasch.scssd,tpl.scssd,tplt.scssd,thpl.scssd),list("rasch","2pl","2plt","3pl"))
plotQQ(Y,list(rasch.scssd,tpl.scssd,tplt.scssd,thpl.scssd),list("rasch","2pl","2plt","3pl"))
chi_sq(Y,list(rasch.scssd,tpl.scssd,tplt.scssd,thpl.scssd),list("rasch","2pl","2plt","3pl"))

#save.image(file = "mdls.RData")