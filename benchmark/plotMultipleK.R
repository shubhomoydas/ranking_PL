rm(list=ls())

inpath <- 'report'
outpath <- 'report'
anomalypath <- 'benchmark'

filenames <- c(
  #'PL_Shuttle_Gs_K1_N4000-inferred-GMM.csv',
  #'PL_Shuttle_Gs_K1_N4000-inferred-noprior.csv',
  #'PL_Shuttle_Gs_K1_N4000-inferred-prior.csv',
  'PL_Shuttle_Gs_K2_N4000-inferred-noprior.csv'
  #'PL_Shuttle_Gs_K2_N4000-inferred-prior.csv'
)

anomalyids <- as.matrix(read.csv(file.path(anomalypath,'shuttle-anomaly-ids.csv'),header=T))
anomalyids <- anomalyids[order(anomalyids[,1]),]

D <- length(filenames)
#d <- 1
for (d in 1:D) {
  data <- read.csv(file.path(inpath,filenames[d]),header=F)
  maxepoch <- max(data[,1])
  K <- ncol(data)-1
  gammas <- as.matrix(data[data[,1]==maxepoch,2:(K+1)])
  N <- nrow(gammas)
  ss <- diag(K)
  diag(ss) <- 1/as.matrix(apply(gammas,2,sum))
  gammas <- gammas %*% ss

  scores <- matrix(0,nrow=N,ncol=2)
  scores[,1] <- 1:N
  scores[,2] <- apply(gammas,1,max)
  scores <- scores[order(-scores[,2]),]
  scores[,2] <- 1:N

  if (d == 1) {
    rankings <- matrix(0,nrow=N,ncol=D+1)
    rankings[,1] <- 1:N
  }

  if (K > 1) {
    subranks <- matrix(0,nrow=N,ncol=K+1)
    subranks[,1] <- 1:N
    for (k in 1:K) {
      tmp <- 1:N
      tmp <- tmp[order(-gammas[,k])]
      subranks[,(k+1)] <- subranks[order(tmp),1]
    }
    write.table(subranks,
      file=file.path(outpath,paste("SUBRANKS_",filenames[d],sep="")),
      col.names=F,
      row.names=F,sep=",",quote=F)
    if (K == 2) {
      pdf(file.path(outpath,paste("Scatter-",filenames[d],".pdf",sep="")))
      plot(subranks[,2:3],typ="p",xlab="Ranking 1",ylab="Ranking 2")
      points(subranks[anomalyids,2:3],col="red",pch="*",cex=2)
      dev.off()
    }
  }

  rankings[,1+d] <- scores[order(scores[,1]),2]
}

