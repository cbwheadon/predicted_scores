## Partial Credit Model
rm(list=ls())

#Make sure this is correct so the bugs syntax file is found
setwd(file.path(getwd(), "/predicted_scores"))

library("R2WinBUGS")
library("mcmcplots")
library("inline")
source("wl_cpp.R")
source("plot_function.R")
source("pfunctions.R")
require("plyr")
source("ppmc.R")
source("functions.R")

#Change this to your bugs directory
bugs.directory = "C:/Program Files/WinBUGS14"

Ypath <- file.path(getwd(), "data/geog.dat")
Y <- as.matrix(read.fwf(file=Ypath,header=FALSE,w=rep(1,10)))

#Convert scores to categories
Y <-  Y + 1

n <- nrow(Y)
p <- ncol(Y)

m.beta <- 0.0
s.beta <- 1.0

#Specify response categories
K <- as.numeric(apply(Y,2,max))

data <- list("Y", "n", "p", "K", "m.beta", "s.beta")
monitor <- c("beta", "theta")
bugs.file <- file.path(getwd(), "bugs/pcm.bug")

iter <- 2000
burnin <- 1000
thin <- 10
n.sims <- (iter - burnin)/thin

system.time(winbugsout <- bugs(data=data, inits=NULL, parameters.to.save=monitor,
				model.file=bugs.file,
				n.iter=iter, n.thin=thin, n.burnin=burnin,bugs.directory=bugs.directory))

#Convert categories back to scores
Y <-  Y - 1

plot(winbugsout)
mcmcplot(winbugsout, random=20)

sims <- winbugsout$sims.matrix

#Clean up after winbugs run
rm(winbugsout)

#Loop through sims
n.sims <- length(sims[,1])
cts <- max(K)
cats <- as.matrix(K)
sim.scores <- vector("list", n.sims)
max.scr <- sum(K-1)
scssd <- matrix(0,nrow=max.scr+1,ncol=n.sims)

for(i in 1:n.sims){
	v <- sims[i,]
	theta <- as.numeric(v[grep("theta", names(v))])
	pars <- parlist(v)
	cssd <- matrix(0,nrow=max.scr+1,ncol=n)
	sim.score <- numeric(p*n)

	for(j in 1:n){

		#For each theta calculate probability distribution and simulated score
		ps <- numeric(cts*p)
		resp <- numeric(p)

		for(k in 1:p){
			#For each question calculate probability
			strt <- cts*(k-1)+1
			ed <- strt+cts-1 
			ps[strt:ed]<-my.gpcm(as.numeric(pars[[k]]),theta[j],cts)

			#Simulate score
			pb <- runif(1)
			thresh <- c(cumsum(ps[strt:ed]),1)
			resp[k] <- min(which(thresh>=pb))-1

		}

		#Conditional Probability Distribution
		ps.t <- matrix(ps,nrow=p,ncol=cts,byrow=TRUE)
		cssd[,j] <- wLord(ps.t,cats)
		strt <- (j-1) * p + 1
		ed <- (j-1) * p + p
		sim.score[strt:ed] <- resp
	}

	scssd[,i] <- rowSums(cssd)
	sim.scores[[i]] <- matrix(sim.score,nrow=n,ncol=p,byrow=TRUE)
	cat(round((i/n.sims)*100,0),"%","\n",sep="")
}

#Score distributions
plotCondSumScoreDist(Y,scssd)
plotQQ(Y,scssd)

#Plot PPMC checks
ppmc_checks(sim.scores)

