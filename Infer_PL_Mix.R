#rm(list=ls())

# drchrnd has been borrowed from from:
#  https://stat.ethz.ch/pipermail/r-help/2006-October/115148.html
#
# One way to simulate X ~ Dirichlet(A1, A2, ..., Ad) 
# is to generate d independent gamma variables having equal rate parameter
# (doesn't matter, so why not 1) and shape parameters  A1, A2, ..., Ad
# Then the vector of components divided by their sum is the desired Dirichlet
#a <- c(1,1,1,1)
#n <- 5
# drchrnd(c(1,1,1,1),5)
drchrnd <- function(a,n) {
  p = length(a)
  r = matrix(0, n, p)
  for (k in 1:p)  r[,k] <- rgamma(n, shape=a[k], rate=1)
  y = apply(r,1,sum)
  y[y == 0] = 1
  r = r / y
  return(r)
}

data <- read.csv("PL-mix.csv",header=F)
#data <- data[1:200,]

K <- 3
S <- as.matrix(data[,3:ncol(data)])
D <- nrow(S)
N <- ncol(S)
Z <- matrix(0,ncol=K,nrow=D)

# Initialize Z with Dirichlet samples
set.seed(32567)
Z <- drchrnd(rep(1,K),D)

# Initialize mixture proportions
Pi <- matrix(0,nrow=K,ncol=1)
Pi[,1] <- (apply(Z,2,sum)/D)

# Initialize V's to random
V_saved <- c()
V <- t(drchrnd(rep(1,N),K))
V_tmp <- matrix(0,nrow=N,ncol=K)

Items <- 1:N
Ordered_ranks <- matrix(0,nrow=D,ncol=N)
for (d in 1:D) {
  Ordered_ranks[d,] <- Items[order(S[d,])]
}

maxepochs <- 600
for (epoch in 1:maxepochs) {
  
  # M-Step
  # ------
  
  # MLE of Pi
  Pi[,1] <- (apply(Z,2,sum)/D)

  # MLE of v_kn
  for (k in 1:K) {
    for (n in 1:N) {
      # we assume that we get complete rankings on N items from
      # each detector. Therefore, we can simplify the indicator
      # function in the numerator of the derivation.
      nv <- sum(Z[,k])
      #for (d in 1:D) {
      #  tv <- tv + Z[d,k]*sum(<Indicator function>)
      #}
      tv <- 0
      for (d in 1:D) {
        s_l <- 0
        for (i in 1:N) {
          #d_iN <- S[d,i:N]
          #a <- length(which(d_iN==n))
          #if (a > 0) {
          if (Ordered_ranks[d,n] >= i) { # this check is equivalent to the 'which'
            d_iN <- S[d,i:N]
            s_l <- s_l + (1/sum(V[S[d,i:N],k]))
          }
        }
        tv <- tv + Z[d,k]*s_l
      }
      V_tmp[n,k] <- nv / tv
    }
  }
  V[,] <- V_tmp[,]
  
  # E-Step
  # ------
  for (d in 1:D) {
    for (k in 1:K) {
      tz = Pi[k]
      for (i in 1:N) {
        tz <- tz * V[S[d,i],k] / sum(V[S[d,i:N],k])
      }
      Z[d,k] <- tz
    }
  }
  Z <- Z / (apply(Z,1,sum)) # normalize to 1
  
  if (epoch %% 10 == 0) {
    V_saved <- rbind(V_saved,cbind(rep(epoch,N),V))
    cat("Epoch", epoch, "\n")
  }
}

outpath <- "results"
write.table(V_saved,
            file=file.path(outpath,paste("PL_V_K",as.character(K),".csv",sep="")), 
            col.names=c("Iter",paste("K",1:(ncol(V_saved)-1),sep="")),
            row.names=F,sep=",",quote=F)

# Original V
V_orig <- matrix(c(0.10,0.56,0.30,0.89,0.02,
                    0.50,0.10,0.15,0.05,0.20,
                    0.25,0.30,0.10,0.10,0.25),byrow=T,nrow=3)
V_orig <- t(V_orig) %*% diag((1/apply(V_orig,1,sum)))
V_orig
V %*% diag((1/apply(V,2,sum)))
