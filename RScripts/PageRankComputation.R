#
# Example implementations to compute the PageRank vector, and to display how the theoretical proposed 
# solutions and acceleration methods may be implemented in reality.
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#------------------------------------------------------------------------------------------------------------
# Algorithms to compute the PageRank Vector

# computes the PageRank vector to the given error > 0 with the absolute norm.
# ALPHA must be chosen from [0,1] and is the weight for the LinkMatrix.
computePageRankByError <- function(LinkMatrix, ALPHA, ERROR) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankVector <- matrix(1/numberOfPages, 1, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  currentError <- ERROR + 1
  numberOfIterations <- 0
  while(ERROR <= currentError) {
    tmp <- PageRankVector
    PageRankVector <- pageRankStep(PageRankVector, LinkMatrix, danglingPagesIndicator, numberOfPages,ALPHA)
    currentError <- sum(abs(tmp - PageRankVector))
    numberOfIterations <- numberOfIterations + 1
  }
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("POWER ITERATION COMPLETE")
  return(PageRankVector)
}

computePageRankAitken <- function(LinkMatrix, ALPHA, ERROR, ExtrapolateAtIterationStep) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  currentError <- ERROR + 1
  numberOfIterations <- 0
  
  if(ERROR <= currentError) {
    while(ERROR <= currentError) {
      if(numberOfIterations == ExtrapolateAtIterationStep) {
        PageRankSaver[3,] <- aitkenExtrapolation(PageRankSaver[1,], PageRankSaver[2,], PageRankSaver[3,])
      }
      PageRankSaver[1,] <- PageRankSaver[2,]
      PageRankSaver[2,] <- PageRankSaver[3,]
      PageRankSaver[3,] <- pageRankStep(PageRankSaver[3,], LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
      currentError <- sum(abs(PageRankSaver[3,] - PageRankSaver[2,]))
      numberOfIterations <- numberOfIterations + 1
    }
  }
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH ONE AITKEN EXTRAPOLATION COMPLETE")
  return(PageRankSaver[3,])
}

# computes the PageRank with t
computePageRankLeastSquare <- function(LinkMatrix, ALPHA, ERROR, vectorSaves, leastSquareFrequency, Extrapolsteps) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, vectorSaves, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  frequencyCounter <- 1
  currentError <- ERROR + 1
  numberOfIterations <- 0
  numberExtrapolSteps <- 0
  while(ERROR <= currentError){
    for(i in 2:vectorSaves) {
      PageRankSaver[i - 1, ] <- PageRankSaver[i,]
    }
    PageRankSaver[vectorSaves,] <- pageRankStep(PageRankSaver[vectorSaves,], LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
    currentError <- sum(abs(PageRankSaver[vectorSaves,] - PageRankSaver[vectorSaves - 1,]))
    
    if((numberOfIterations  + vectorSaves + 1) %%  leastSquareFrequency == 0 && numberOfIterations > 0 && numberExtrapolSteps < Extrapolsteps) {
      PageRankSaver[vectorSaves,] <- quadraticExtrapolation(PageRankSaver)
      numberExtrapolSteps <- numberExtrapolSteps + 1
    }
    numberOfIterations <- numberOfIterations + 1
  }

  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH LEAST SQUARE COMPLETE")
  return(PageRankSaver[vectorSaves,])
}

# 
computePageRankTSS <- function(LinkMatrix, ALPHA, BETA, ERROR, ErrorInner) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[1,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  currentError <- ERROR + 1
  
  while(ERROR <= currentError) {
    PageRankSaver[1,] <- PageRankSaver[2,]
    reverseSplitting <- 
      Residual <- ErrorInner + 1
    while(ErrorInner <= sum(abs(Residual - PageRankSaver[2,]))) {
      PageRankSaver[2,] <- (ALPHA * PageRankSaver[3,] 
                            + ALPHA * PageRankSaver[3,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
                            + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages))
      PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
      Residual <- (ALPHA - BETA) * PageRankSaver[3,] + (1 - ALPHA) * matrix(1 / numberOfPages, 1, numberOfPages)
      numberOfInnerIterations <- numberOfInnerIterations + 1
    }
    PageRankSaver[2,] <- Residual + BETA * PageRankSaver[3,]
    PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
    currentError <- sum(abs(PageRankSaver[2,] - PageRankSaver[1,]))
    
    numberOfIterations <- numberOfIterations + 1
  }
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH TWO STEP SPLITTING ITERATION COMPLETE")
  return(PageRankSaver[3,])
}

# 
computePageRankRTSS <- function(LinkMatrix, ALPHA, BETA, GAMMA, ERROR, ErrorInner) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[1,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  currentError <- ERROR + 1
  
  while(ERROR <= currentError) {
    PageRankSaver[1,] <- PageRankSaver[2,] 
    Residual <- ErrorInner + 1
    while(ErrorInner <= sum(abs(Residual - PageRankSaver[2,]))) {
      PageRankSaver[2,] <- (ALPHA * PageRankSaver[3,] 
                            + ALPHA * PageRankSaver[3,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
                            + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages))
      PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
      Residual <- ((GAMMA - 1) / GAMMA * PageRankSaver[2,] 
                   + (ALPHA - BETA) / GAMMA * PageRankSaver[3,] 
                   + (1 - ALPHA) / GAMMA * matrix(1/numberOfPages, 1, numberOfPages))
      numberOfInnerIterations <- numberOfInnerIterations + 1
    }
    PageRankSaver[2,] <- Residual + BETA / GAMMA * PageRankSaver[3,]
    PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
    currentError <- sum(abs(PageRankSaver[2,] - PageRankSaver[1,]))
    
    numberOfIterations <- numberOfIterations + 1
  }
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH RELAXED TWO STEP SPLITTING ITERATION COMPLETE")
  return(PageRankSaver[3,])
}

# computes the PageRank by a number of given iterations. Simple method for acquiring a close approximation
# of the PageRank Vector.
computePageRankByIterations <- function(LinkMatrix, ALPHA, nIterations) {
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 1, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  for(i in 1:nIterations) {
    PageRankSaver <- pageRankStep(PageRankSaver, LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
  }
  return(PageRankSaver)
}

