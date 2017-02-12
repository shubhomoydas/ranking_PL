
# Load the individual scores from the egmm output for Shuttle data and convert them
# to ranks so that each component in the ensemble can act as a separate detector.

inpath <- 'report'
outpath <- 'report'

filenames <- c(
  'PL_Shuttle_Gs_K1_N4000-inferred-GMM.csv',
  'PL_Shuttle_Gs_K1_N4000-inferred-noprior.csv',
  'PL_Shuttle_Gs_K1_N4000-inferred-prior.csv',
  'PL_Shuttle_Gs_K2_N4000-inferred-noprior.csv',
  'PL_Shuttle_Gs_K2_N4000-inferred-prior.csv'
)

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

  rankings[,1+d] <- scores[order(scores[,1]),2]
}

write.table(rankings,
  file=file.path(outpath,"shuttle-gamma-ranks.csv"),
  col.names=F,
  row.names=F,sep=",",quote=F)

