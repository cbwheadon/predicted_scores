## Rasch Model
rm(list=ls())

#Make sure this is correct so the bugs syntax file is found
setwd(file.path(getwd(), "/predicted_scores"))

library("R2WinBUGS")
library("mcmcplots")
library("inline")
source("wl_cpp.R")
source("plot_function.R")

#Change this to your bugs directory
bugs.directory = "C:/Program Files/WinBUGS14"

Ypath <- file.path(getwd(), "data/chem.csv")
Y <- as.matrix(read.csv(file=Ypath,header=FALSE,nrows=200))

n <- nrow(Y)
p <- ncol(Y)

m.delta <- 0.0
s.delta <- 1.0
m.alpha <- 1.0
s.alpha <- 1.0

iter <- 2000
burnin <- 1000
thin <- 10

data <- list("Y", "n", "p",
		"m.alpha", "s.alpha",
		"m.delta", "s.delta")
monitor <- c("delta", "theta", "alpha")
bugs.file <- file.path(getwd(), "bugs/tpl.bug")

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
	
	probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) exp(alpha[j]*(th-delta[j]))/(1+exp(alpha[j]*(th-delta[j])))))) 
	
	for (s in 1:n){
		ps <- matrix(cbind(1-probs[s,],probs[s,]),ncol=2)
		cond.sum.sc.d <- wLord(ps,cats)
		if(s==1){
			scssd[,1] <- cond.sum.sc.d
		} else {
			scssd[,i] <- scssd[,i] + cond.sum.sc.d
		}
	}
	setTxtProgressBar(pb, i)
}

close(pb)
plotCondSumScoreDist(Y,scssd,"2pl")
plotQQ(Y,scssd,"2pl")