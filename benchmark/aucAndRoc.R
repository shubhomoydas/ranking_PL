rm(list=ls())

inpath <- 'C:/smd/git/repository/ADAMS/code/RankingBTPL/benchmark'
outpath <- 'C:/smd/git/repository/ADAMS/code/RankingBTPL/benchmark'

#inpath <- '/nfs/guille/u2/d/dassh/codebase/repository/ADAMS/code/RankingBTPL/benchmark'
#outpath <- '/nfs/guille/u2/d/dassh/classes/ECE599/report'

# For benchmark csv file, following is the format of rows:
# First row is header and the order of columns is:
#   aggRank,id,meanpdf,pdfGMM_1,pdfGMM_2,...,pdfGMM_D
# here, 'GMM' - Density Estimator using Gaussian Mixture Model
benchmarkdata <- read.csv(file.path(inpath,'shuttle-anomalies.csv'),header=T)

# Ensemble data has NO header row.
# The order of columns is:
#   id,rankGMM_1,rankGMM_2,...,rankGMM_D
ensembledata <- read.csv(file.path(inpath,'shuttle-gamma-ranks.csv'),header=F)

answerfile <- file.path(inpath,'shuttle-anomaly-ids.csv')
rocfile <- file.path(outpath,'shuttle-anomaly-roc.png')

D <- (ncol(ensembledata) - 1) + 1
N <- nrow(benchmarkdata)

plotROC <- T

# Order by id
benchmarkdata = benchmarkdata[order(benchmarkdata[,2]),] # second column is id in benchmark
ensembledata = ensembledata[order(ensembledata[,1]),] # first column is id in ensemble

rankings <- matrix(0,nrow=N,ncol=1+D)
rankings[,1] <- benchmarkdata[,2]
rankings[,2] <- benchmarkdata[,1]

for (d in 1:(D-1)) {
  rankings[,2+d] <- ensembledata[,1+d]
}
rankings <- rankings[,-c(5,7)] # exclude the EM-Priors

source(file.path(inpath,'lift-roc.R'))
fnCalcLiftAndDispROC(rankings,answerfile,plotROC,rocfile)

