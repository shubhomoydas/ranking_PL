
# Plots the ranks for actual anomalies to see if there are distinct
# patterns in ranking among the detectors.

rm(list=ls())

# install.packages("plotrix",lib=c("/nfs/guille/u2/d/dassh/R"))

#require("plotrix") # Use this on Windows local
library(lattice)
require("plotrix",lib=c("/nfs/guille/u2/d/dassh/R")) # Use this in Linux Env

inpath <- '/nfs/guille/u2/d/dassh/codebase/repository/ADAMS/code/RankingBTPL/benchmark'
outpath <- '/nfs/guille/u2/d/dassh/classes/ECE599/report'

benchmarkdata <- as.matrix(read.csv(file.path(inpath,'shuttle-anomalies.csv'),header=T))
ensembledata <- as.matrix(read.csv(file.path(inpath,'shuttle-gamma-ranks.csv'),header=F))

benchmarkdata <- benchmarkdata[order(benchmarkdata[,2]),]
ensembledata <- ensembledata[order(ensembledata[,1]),]

benchmarkdata[1:20,1:3]
ensembledata[1:20,]

anomalyids <- as.matrix(read.csv(file.path(inpath,'shuttle-anomaly-ids.csv'),header=T))
anomalyids <- anomalyids[order(anomalyids[,1]),]

D = ncol(benchmarkdata) - 3 # first three columns are not to be used here...
N <- nrow(benchmarkdata)

rankings <- matrix(0,nrow=N,ncol=1+D)
rankings[,1] <- benchmarkdata[,2]

T <- matrix(0,nrow=N,ncol=2)
for (d in 1:D) {
  T[,1] <- benchmarkdata[order(benchmarkdata[,d+3]),2]
  T[,2] <- 1:N
  T <- T[order(T[,1]),]
  rankings[,d+1] <- T[,2]
}

anomrankings <- rankings[anomalyids,]

# Following is a lazy 'hack' to get D distinct colors
allcolors <- colors()
allcolors <- allcolors[-grep("grey",allcolors)]
allcolors <- allcolors[-grep("gray",allcolors)]
allcolors <- allcolors[-grep("black",allcolors)]
allcolors <- allcolors[-grep("white",allcolors)]
set.seed(32767)
colidxs <- sample(1:length(allcolors),D)


ltys <- c("dashed","dotted","dotdash","twodash","longdash") #"solid",

maxrank <- max(anomalyids)
tmpranks <- rbind(rep(1,D),rep(maxrank,D),anomrankings[,2:(D+1)])
cols <- c("white","white",allcolors[colidxs])
detnames <- paste("",1:D,sep="")
anomnames <- c("","",paste("Item ",1:D,sep=""))

#png(file.path(outpath,"anomalies-parallel.png"),width=8,height=4,units="in",res=600)
pdf(file.path(outpath,"anomalies-parallel.pdf"),width=8,height=4)
  par(cex.lab=1)
  print(
    parallelplot(~tmpranks, horizontal.axis=FALSE, main="Anomaly Detector Ranks",
                 col=cols, lty=ltys, 
                 cex.lab=1, cex=1, cex.axis=1, cex.sub=1, cex.main=1,
                 xlab="Detector",ylab="Ranks",
                 scales=list(cex=1,cex.lab=1,cex.axis=1,cex.sub=1),
                 lwd=2, varnames=detnames,ylim=c(0,1)
                 #,auto.key=F
                 ,key=list(
                   corner=c(1,0),
                   cex=0.5,
                   corner=c(0.1,0.99),
                   border=TRUE,
                   lines=list(
                     lty=ltys,lwd=2,col=cols
                   )
                   ,text=list(anomnames)
                 )
    )
  )
dev.off()

