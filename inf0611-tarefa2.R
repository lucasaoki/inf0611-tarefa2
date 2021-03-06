#######################################
## INF-0611 - Trabalho 2 - Questao 1 ##
## Alunos: Felipe Wolff Ramos        ##
##         Lucas Aoki Heredia        ##
#######################################

library(ggplot2)

## Séries temporais fornecidas
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

## Tabela pre-definida para discretização das séries normalizadas
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

## Função para normalização das séries
normalize <- function(data) {
  meanData <- mean(data)
  sdData <- sd(data)
  return((data - meanData) / sdData)
}

## Função de algoritmo do Piecewise Aggregate Approximation
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

## Função para discretização em símbolos
convertToEqualProbSymbols <- function(data, partitions) {
  partLimits <- gaussianEqProbPartitions(partitions)
  result <- seq(from=1, by=0, length.out=length(data))
  for (i in c(1:length(partLimits))) {
    result <- sumIfInRange(result, data, min=partLimits[i])
  }
  result <- letters[result]
  return(result)
}

## Calculo da distancia entre símbolos
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

## Função para calculo da MINDIST entre series
minDist <- function(data1, data2, dimens, origLength, partitions) {
  combinationSize <- origLength / dimens
  result <- 0
  for (i in c(1:length(data1))) {
    result <- result + (symbolsDist(data1[i],data2[i],partitions))^2
  }
  result <- sqrt(result)*sqrt(combinationSize)
  return(result)
}

## Função para gerar as labels dos gráficos
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

## Função para gerar os gráficos de comparação das séries dado o número de 
## partições e dimensões.
plotDiscrete <- function(ANorm, BNorm, numPartitions, dimens) {
  size <- length(ANorm)
  points <- c(1:size)
  partitions <- gaussianEqProbPartitions(numPartitions)
  # aplica PAA
  APaa <- paa(ANorm, dimens)
  BPaa <- paa(BNorm, dimens)
  
  APaaPlot <- rep(APaa, times=1, each=6)
  BPaaPlot <- rep(BPaa, times=1, each=6)
  Alabels <- genLabelsPerPartition(APaa, numPartitions, size)
  Blabels <- genLabelsPerPartition(BPaa, numPartitions, size)
  
  Norm <- c(ANorm, BNorm)
  PaaPlot <- c(APaaPlot, BPaaPlot)
  labels <- c(Alabels, Blabels)
  keys <- rep(c("A","B"), times=1, each=length(A))
  points <- c(points,points)
  dFrame <- data.frame(Norm,PaaPlot,labels,points,keys)
  
  g <- ggplot(data = dFrame, aes(x=points)) + 
    labs(title=paste("Séries com", numPartitions, "símbolos"), y="Valores Normalizados", x="Número da medição", colour="Séries") +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    geom_hline(yintercept = partitions, alpha = 0.4, colour="#003366") +
    scale_y_continuous(breaks=partitions) + 
    scale_colour_manual(values=c("#990000", "#006633"))

  g <- g +
    geom_line(aes(y=PaaPlot, colour=keys)) +
    geom_line(aes(x=points, y=Norm, colour=keys), alpha=0.5) +
    geom_text(aes(x=points, y=PaaPlot, label=labels, colour=keys), nudge_y = 0.08, show.legend = F)
  
  return(g)
}

# tamanho final da dimensão
dimens <- 24
# séries normalizadas
ANorm <- normalize(A)
BNorm <- normalize(B)
# PAA aplicada
APaa <- paa(ANorm, dimens)
BPaa <- paa(BNorm, dimens)

# discretizar e calcular distancias
for (i in c(4:7)) {
  ASym <- convertToEqualProbSymbols(APaa,i)
  BSym <- convertToEqualProbSymbols(BPaa,i)
  
  print(minDist(ASym,BSym,dimens,length(A),i))
}

# gerar gráficos
plotDiscrete(ANorm, BNorm, 4, dimens)
plotDiscrete(ANorm, BNorm, 5, dimens)
plotDiscrete(ANorm, BNorm, 6, dimens)
plotDiscrete(ANorm, BNorm, 7, dimens)
