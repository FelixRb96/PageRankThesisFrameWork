#
# Example implementations to compute the PageRank vector, and to display how the theoretical proposed 
# solutions and acceleration methods may be implemented in reality.
#
# All methods are stripped of the tracking tools to provide optimal performance.
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#------------------------------------------------------------------------------------------------------------
# Algorithms to compute the PageRank Vector

# computes the PageRank vector to the given error > 0 with the absolute norm.
# ALPHA must be chosen from [0,1] and is the weight for the LinkMatrix.
computePageRankByError <- function(LinkMatrix, ALPHA, ERROR) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
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

# computes the PageRank vector by applying ONE Aitken delta-squared
# step at the given iteration step.
computePageRankAitken <- function(LinkMatrix, ALPHA, ERROR, ExtrapolateAtIterationStep) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
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

# computes the PageRank with the least-square extrapolation.
# can be done with a chosen number of extrapolation vectors, a maximum number of extrapolation
# steps is necessary to provide shure convergence.
#
# WARNING: choosing a high number for vectorsaves can cause ERRORS for the least-square step, since the 
# saved vectors are in a cauchy-sequence, thus approach each other leading to stability problems.
computePageRankLeastSquare <- function(LinkMatrix, ALPHA, ERROR, vectorSaves, leastSquareFrequency, Extrapolsteps) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
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


# computes the PageRank vector with the implementation of the TSS ALgorithm
computePageRankTSS <- function(LinkMatrix, ALPHA, BETA, ERROR, ERRORINNER) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * matrix(1/numberOfPages, 1, numberOfPages))
  
  Residual <- c()
  
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  currentError <- ERROR + 1
  
  while(ERROR <= currentError) {
    
    PageRankSaver[2,] <- ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages)
    PageRankSaver[2,] <- PageRankSaver[2,] / sum(abs(PageRankSaver[2,]))
    PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * 1/numberOfPages)
    PageRankSaver[3,] <- PageRankSaver[3,] / sum(abs(PageRankSaver[3,]))
    Residual <- (ALPHA - BETA) * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages)
    
    currentErrorInner <- ERRORINNER + 1
    
    while(ERRORINNER <= currentErrorInner) {
      PageRankSaver[1,] <- PageRankSaver[2,]
      PageRankSaver[2,] <- Residual + BETA* PageRankSaver[3,]
      PageRankSaver[2,] <- PageRankSaver[2,] / sum(abs(PageRankSaver[2,]))
      PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * 1/numberOfPages)
      PageRankSaver[3,] <- PageRankSaver[3,] / sum(abs(PageRankSaver[3,]))
      currentErrorInner <- sum(abs(Residual + BETA * PageRankSaver[3,] - PageRankSaver[2,]))
 
      numberOfInnerIterations <- numberOfInnerIterations + 1
    }
    
    currentError <- sum(abs(ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages) - PageRankSaver[2,]))

    numberOfIterations <- numberOfIterations + 1
  }
  
  PageRankSaver[2,] <- ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages)
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH TWO STEP SPLITTING ITERATION COMPLETE")
  return(PageRankSaver[2,] / sum(abs(PageRankSaver[2,])))
}

# computes the PageRank vector with the implementation of the RTSS ALgorithm
# choosing 0 < BETA < ALPHA < GAMMA <= 1 is necessary.
#
# choosing GAMMA = 1 produces the TSS algorithm with some additional computation.
computePageRankRTSS <- function(LinkMatrix, ALPHA, BETA, GAMMA, ERROR, ERRORINNER) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * matrix(1/numberOfPages, 1, numberOfPages))
  
  Residual <- c()
  
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  currentError <- ERROR + 1
  
  while(ERROR <= currentError) {
    
    PageRankSaver[2,] <- ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages)
    PageRankSaver[2,] <- PageRankSaver[2,] / sum(abs(PageRankSaver[2,]))
    PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * 1/numberOfPages)
    PageRankSaver[3,] <- PageRankSaver[3,] / sum(abs(PageRankSaver[3,]))
    Residual <- (GAMMA - 1) / GAMMA * PageRankSaver[2,] + (ALPHA - BETA) / GAMMA * PageRankSaver[3,] + (1 - ALPHA) / GAMMA * matrix(1/numberOfPages, 1, numberOfPages)
    
    currentErrorInner <- ERRORINNER + 1
    
    while(ERRORINNER <= currentErrorInner) {
      PageRankSaver[1,] <- PageRankSaver[2,]
      PageRankSaver[2,] <- Residual + BETA / GAMMA * PageRankSaver[3,]
      PageRankSaver[2,] <- PageRankSaver[2,] / sum(abs(PageRankSaver[2,]))
      PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * 1/numberOfPages)
      PageRankSaver[3,] <- PageRankSaver[3,] / sum(abs(PageRankSaver[3,]))
      currentErrorInner <- sum(abs(Residual + BETA / GAMMA * PageRankSaver[3,] - PageRankSaver[2,]))
      
      numberOfInnerIterations <- numberOfInnerIterations + 1
    }
    
    currentError <- sum(abs(ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages) - PageRankSaver[2,]))
    
    numberOfIterations <- numberOfIterations + 1
  }
  
  PageRankSaver[2,] <- ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages)
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH RELAXED TWO STEP SPLITTING ITERATION COMPLETE")
  return(PageRankSaver[2,] / sum(abs(PageRankSaver[2,])))
}


# computes the PageRank by a number of given iterations. Simple method for acquiring a close approximation
# of the PageRank Vector.
computePageRankByIterations <- function(LinkMatrix, ALPHA, nIterations) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 1, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  for(i in 1:nIterations) {
    PageRankSaver <- pageRankStep(PageRankSaver, LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
  }
  return(PageRankSaver)
}

