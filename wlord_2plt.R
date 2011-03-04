## 2 parameter testlet model
rm(list=ls())

#Make sure this is correct so the bugs syntax file is found
setwd(file.path(getwd(), "/predicted_scores"))

library("R2WinBUGS")
library("mcmcplots")
library("inline")
source("wl_cpp.R")
source("plot_gg.R")

#Change this to your bugs directory
bugs.directory = "C:/Users/User/Winbugs/WinBUGS14"

Ypath <- file.path(getwd(), "data/chem.csv")
Y <- as.matrix(read.csv(file=Ypath,header=FALSE,nrows=500))

n <- nrow(Y)
p <- ncol(Y)

m.alpha <- 1.0
s.alpha <- 1.0
m.delta <- 0.0
s.delta <- 1.0
a.sigsq.gamma <- 1.0
b.sigsq.gamma <- 1.0

#Change this to represent testlet structure
d <- c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4),rep(7,4))

#Number of testlets
n.t <- max(d) 

iter <- 2000
burnin <- 1000
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

#plot(winbugsout)
mcmcplot(winbugsout, random=20)

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
plotCondSumScoreDist(Y,scssd,"2plt")
plotQQ(Y,scssd,"2plt")