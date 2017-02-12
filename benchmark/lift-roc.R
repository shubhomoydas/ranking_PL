#########################################################
##  Computes lift and AUCs
## rankings:
##   Col_1 - Id
##   Col_2 - Ranks under algo 1
##   Col_3 - Ranks under algo 2
##   Col_4 - ...
#########################################################

fnCalcLiftAndDispROC <- function (rankings, answerfile, plotROC, rocfile) {

  cat("Computing Lift/AUCs...\n")
  anomalydata <- read.csv(answerfile, header=T, sep=",")
  anomalyids <- anomalydata[,1]
  #print(anomalyids)
  
  D <- ncol(rankings)-1
  N <- nrow(rankings)
  # Reorder the data on userid after ranking
  rankingsOrderedByIds <- rankings[order(rankings[,1]),]
  
  # Creating a separate variable tdata so that if we want, we
  # can append results of other algorithms as columns on the right.
  # For now, we have only one input datafile and hence only one algorithm data
  # Col 1 - Id
  # Col 2 - Ground Truth Anomaly/Non-Anomaly
  # Col 3 - Detector 1 ranks
  # Col 4 - Detector 2 ranks...
  # ...
  tdata <- matrix(0,nrow=N,ncol=D+2)
  tdata[,1] <- as.numeric(rankingsOrderedByIds[,1]) # ids
  for (d in 1:D) {
    tdata[,d+2] <- rankingsOrderedByIds[,d+1]
  }
  
  # At this point, mark all the anomalous instances
  for (i in 1:length(anomalyids)) {
    tdata[rankingsOrderedByIds[,1] == anomalyids[i],2] <- 1 # mark the anomalies
  }
  
  # Show the marked anomalies to make sure we did it right
  cat("Following are the actual anomalies:\n")
  print(tdata[tdata[,2]==1,])
  
  # average lifts will be stored here
  avglifts = rep(0,D)
  n <- sum(as.numeric(tdata[,2])) # total number of anomalies
  
  for (d in  1:D) {
    ranked <- tdata[order(tdata[,2+d]),] # order the data on the basis of ranking
    found  <- ranked[ranked[,2]==1,] # grab the rows for anomalous users
    found  <- found[order(found[,2+d]),]
    r <- cumsum(ranked[,2]) / n # lift for uniform random algorithm
    
    lifts = cbind(r[found[,2+d]],found[,2+d],n*r[found[,2+d]] / (n*found[,2+d]/length(r)))
    avglifts[d] <- mean(lifts[,3])
  }
  
  # dispay the avg lifts
  cat("Lifts:\n")
  print(avglifts)
  
  if (plotROC) {
    # required for ROC curves
    require("pROC")
    
    ltys <- c("solid","solid",
              "solid", #"dashed",
              "solid"  #"dashed"
              )
    cols <- c("black","red",
              "blue", #"blue",
              "green" #,"green"
              )
    lwds <- rep(2,D)
    
    algos <- c(
      "Benchmark"
      , "GMM"
      , "EM-K1-NoPrior" #, "EM-K1-Prior"
      , "EM-K2-NoPrior" #, "EM-K2-Prior"
    )
    
    aucs <- rep(0,D)
    # Plot ROC
    png(rocfile)
    # more options, CI and plotting
    roc1 <- roc(tdata[,1+1],
                tdata[,2+1], percent=TRUE,
                smooth=FALSE,
                # arguments for auc
                auc=TRUE,
                #partial.auc=c(100, 90), partial.auc.correct=TRUE,
                #partial.auc.focus="sens",
                # arguments for ci
                #ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=TRUE,
                print.auc=FALSE,lty=ltys[1],col=cols[1],lwd=lwds[1],cex.lab=2)
    aucs[1] <- auc(roc1)
    
    if (D > 1) {
      allROCs <- c(roc1)
      for (d in 2:D) {
        roc2 <- roc(tdata[,1+1], tdata[,2+d],
                    plot=TRUE, add=TRUE, 
                    percent=roc1$percent, print.auc=FALSE,
                    lty=ltys[d],col=cols[d],lwd=lwds[d],cex.lab=2)
        aucs[d] <- auc(roc2)
        allROCs <- c(allROCs, roc2)
      }
    }
    
    legends <- algos
    for (d in 1:D) {
      legends[d] <- paste(legends[d]," (",round(aucs[d],2),")",sep="")
    }
    legend(
      "bottomright",cex=1.5,
      legend=legends,
      lty=ltys,
      col=cols,
      lwd=2
    )
    dev.off()
    
    # display AUC
    cat("AUCs:\n")
    print(aucs)
  }
  
}
