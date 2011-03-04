## Partial Credit Model
rm(list=ls())

#Make sure this is correct so the bugs syntax file is found
setwd(file.path(getwd(), "/predicted_scores"))

library("R2WinBUGS")
library("mcmcplots")
library("inline")
source("wl_cpp.R")
source("plot_gg.R")
source("pfunctions.R")
require("plyr")

#Change this to your bugs directory
bugs.directory = "C:/Users/User/Winbugs/WinBUGS14"

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

data <- list("Y", "n", "p", "K",
		"m.beta", "s.beta")
monitor <- c("beta", "theta")
bugs.file <- file.path(getwd(), "bugs/pcm.bug")

iter <- 2000
burnin <- 1000
thin <- 10

system.time(winbugsout <- bugs(data=data, inits=NULL, parameters.to.save=monitor,
				model.file=bugs.file,
				n.iter=iter, n.thin=thin, n.burnin=burnin,bugs.directory=bugs.directory))

plot(winbugsout)
mcmcplot(winbugsout, random=20)

sims <- winbugsout$sims.matrix

##Convert categories back to scores
Y <-  Y - 1

#Loop through sims
n.sims <- length(sims[,1])
cts <- max(K)
cats <- as.matrix(K)
max.scr <- sum(K-1)

#Score distributions
scssd <- sum.prob(sims,cts,cats,max.scr)

#Score distributions
plotCondSumScoreDist(Y,scssd,"pcm")
plotQQ(Y,scssd,"pcm")