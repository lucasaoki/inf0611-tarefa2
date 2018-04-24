####################################
## INF-0611 - Trabalho 2 - Questao 1
## Aluno: Felipe Wolff Ramos
####################################

library(ggplot2)

A <- c(21.7, 21.7, 21.6, 21.6, 21.7, 21.7, 21.7, 21.6, 21.5, 21.5, 21.4, 21.2,
       21.2, 21.1, 21.0, 20.9, 20.9, 21.0, 20.9, 20.9, 20.8, 20.7, 20.6, 20.6,
       20.5, 20.5, 20.5, 20.5, 20.5, 20.4, 20.3, 20.2, 20.1, 20.0, 20.0, 20.0,
       20.0, 19.9, 19.8, 19.8, 19.8, 20.0, 20.3, 20.8, 21.1, 21.7, 22.3, 22.6,
       23.0, 23.8, 24.4, 24.8, 24.7, 25.1, 25.8, 26.3, 26.6, 26.5, 27.0, 27.2,
       27.6, 27.6, 27.9, 28.1, 28.2, 28.2, 28.6, 29.0, 29.0, 29.1, 29.4, 29.4,
       29.5, 29.5, 29.6, 30.1, 30.1, 30.4, 30.2, 30.5, 30.6, 30.4, 30.6, 30.2,
       30.4, 30.6, 30.1, 30.2, 30.3, 30.2, 30.3, 30.5, 30.1, 30.0, 30.3, 31.1,
       31.2, 31.1, 31.2, 31.3, 31.6, 31.3, 30.8, 30.0, 30.5, 29.9, 29.7, 29.9,
       29.2, 28.7, 28.4, 28.2, 26.4, 25.0, 24.4, 23.9, 23.7, 23.7, 23.8, 23.9,
       23.9, 23.8, 24.0, 24.1, 24.2, 24.2, 24.1, 24.1, 24.0, 24.0, 24.0, 24.0,
       23.9, 23.6, 23.4, 23.4, 23.4, 23.3, 23.2, 23.1, 23.0, 22.9, 22.9, 22.8)
B <- c(21.4, 21.3, 21.3, 20.9, 20.4, 20.0, 19.8, 19.9, 19.9, 19.7, 20.0, 19.8, 
       19.7, 20.1, 20.1, 19.9, 19.7, 18.8, 19.0, 18.3, 18.0, 17.5, 17.4, 17.5, 
       17.7, 18.0, 18.0, 17.5, 17.5, 17.7, 18.1, 18.0, 17.9, 17.6, 17.2, 17.3, 
       17.5, 17.1, 17.2, 17.5, 17.4, 17.7, 18.0, 18.0, 17.8, 17.7, 17.6, 17.9, 
       19.3, 20.2, 20.6, 21.6, 22.3, 21.7, 21.5, 21.7, 22.2, 22.4, 22.6, 23.1, 
       23.4, 24.0, 24.1, 24.5, 24.8, 25.0, 25.7, 25.8, 25.8, 26.4, 26.6, 27.0, 
       26.8, 26.9, 27.0, 27.3, 27.1, 27.8, 28.0, 28.2, 28.2, 27.9, 27.4, 27.2, 
       27.2, 27.3, 27.2, 27.1, 27.4, 27.7, 27.4, 27.3, 27.2, 27.7, 27.8, 28.2, 
       28.0, 27.8, 27.7, 27.7, 27.7, 27.8, 27.5, 26.6, 25.7, 25.0, 24.2, 23.5, 
       23.2, 22.9, 22.5, 22.3, 22.0, 21.6, 21.3, 21.0, 20.8, 20.4, 20.3, 20.0, 
       19.7, 19.5, 19.3, 19.1, 19.0, 18.9, 18.7, 18.6, 18.5, 18.4, 18.4, 18.4, 
       18.4, 18.3, 18.3, 18.4, 18.4, 18.4, 18.4, 18.3, 18.3, 18.3, 18.4, 18.3)


gaussianEqProbPartitions <- function(numPart) {
  partitions <- list(c(), 
                     c(), 
                     c(-0.43, 0.43), 
                     c(-0.67, 0, 0.67), 
                     c(-0.84, -0.25, 0.25, 0.84),
                     c(-0.97, -0.43, 0, 0.43, 0.97),
                     c(-1.07, -0.57, -0.18, 0.18, 0.57, 1.07))
  return(partitions[[numPart]])
}

normalize <- function(data) {
  meanData <- mean(data)
  sdData <- sd(data)
  return((data - meanData) / sdData)
}

paa <- function(data, dimens) {
  combinationSize <- length(data) / dimens
  result <- c()
  for (i in c(0:(dimens-1))) {
    tmp <- data[(1+i*combinationSize):((1+i)*combinationSize)]
    result <- c(result,(sum(tmp)/combinationSize))
  }
  return(result)
}

sumIfInRange <- function(value, cmp, min = NULL, max = NULL, by=1) {
  if (is.null(min) && is.null(max)) return(value)
  if (!is.null(min) && !is.null(max)) {
    if (cmp >= min && cmp < max) return(value + by)
  } else if (is.null(max)) {
    if (cmp >= min) return(value + by)
  } else {
    if (cmp < max) return (value + by)
  }
  return(value)
}

sumIfInRange <- Vectorize(sumIfInRange)

convertToEqualProbSymbols <- function(data, partitions) {
  partLimits <- gaussianEqProbPartitions(partitions)
  result <- seq(from=1, by=0, length.out=length(data))
  for (i in c(1:length(partLimits))) {
    result <- sumIfInRange(result, data, min=partLimits[i])
  }
  result <- letters[result]
  return(result)
}

symbolsDist <- function(val1 , val2, partitions) {
  partLimits <- gaussianEqProbPartitions(partitions)
  val1 <- match(val1, letters)
  val2 <- match(val2, letters)
  if (abs(val1 - val2) <= 1) {
    return(0)
  } else {
    idxMax <- max(val1,val2) - 1
    idxMin <- min(val1, val2)
    return(partLimits[idxMax] - partLimits[idxMin])
  }  
}

minDist <- function(data1, data2, dimens, origLength, partitions) {
  combinationSize <- origLength / dimens
  result <- 0
  for (i in c(1:length(data1))) {
    result <- result + (symbolsDist(data1[i],data2[i],partitions))^2
  }
  result <- sqrt(result)*sqrt(combinationSize)
  return(result)
}

genLabelsPerPartition <- function(data, partitions, totalSize) {
  symbols <- convertToEqualProbSymbols(data, partitions)
  reps <- totalSize / length(data)
  result <- c()
  for (i in c(1:length(symbols))) {
    tmp <- c()
    for (j in c(1:reps)) {
      if (j == reps/2) {
        tmp <- c(tmp,symbols[i])
      } else {
        tmp <- c(tmp,"")
      }
    }
    result <- c(result,tmp)
  }
  return(result)
}

dimens <- 24

ANorm <- normalize(A)
BNorm <- normalize(B)

APaa <- paa(ANorm, dimens)
BPaa <- paa(BNorm, dimens)

for (i in c(4:7)) {
  ASym <- convertToEqualProbSymbols(APaa,i)
  BSym <- convertToEqualProbSymbols(BPaa,i)
  
  print(minDist(ASym,BSym,dimens,length(A),i))
}

ASym <- convertToEqualProbSymbols(APaa,7)
BSym <- convertToEqualProbSymbols(BPaa,7)

minDist(ASym,BSym,dimens,length(A),7)


#points <- rep(seq(from=0,by = length(A)/dimens,length=length(APaa)), times=1, each=2)
#points <- c(points[2:length(points)],length(A))

#APaaPlot <- rep(APaa, times=1, each=2)
#BPaaPlot <- rep(BPaa, times=1, each=2)

points <- c(1:length(A))
Asym <- convertToEqualProbSymbols(APaa,i)
APaaPlot <- rep(APaa, times=1, each=6)
BPaaPlot <- rep(BPaa, times=1, each=6)
text <- genLabelsPerPartition(APaa,4,length(A))

dFrame <- data.frame(APaaPlot, BPaaPlot ,points, text)
g <- ggplot(data = dFrame)  + geom_line(aes(x=points,y=APaaPlot)) + geom_text(nudge_y = 0.08,aes(x=points, y=APaaPlot, label=text))
