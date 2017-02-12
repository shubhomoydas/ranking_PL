
# Load the individual scores from the egmm output for Shuttle data and convert them
# to ranks so that each component in the ensemble can act as a separate detector.

inpath <- 'benchmark'
outpath <- 'benchmark'
data <- read.csv(file.path(inpath,'shuttle-anomalies.csv'),header=T)

# First three columns should not be converted as scores.
# Second column is the id.
D = ncol(data) - 3
N = nrow(data)

ranked <- matrix(nrow=D,ncol=N)

for (d in 1:D) {
  ranked[d,] <- data[order(data[,d+3]),2]
}

write.table(cbind(rep(1,D),1:D,ranked),
            file=file.path(outpath,"shuttle-ranks.csv"),
            col.names=F,
            row.names=F,sep=",",quote=F)

# Convert data into a format that can be used by ADAMS ensembling
sorteddata <- -data[order(data[,2]),-c(1,3)]
sorteddata[,1] <- -sorteddata[,1]
#head(sorteddata)
write.table(sorteddata,
            file=file.path(outpath,"shuttle-scores.csv"),
            col.names=T,
            row.names=F,sep=",",quote=F)
