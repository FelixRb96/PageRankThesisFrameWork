#
# Same implementations of the PageRank algorithms as in @PageRankComputation, but returning an array 
# of the recorded errors and some further methods to visualize and compare the performance.
# 
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#---------------------------------------------------------------------------------

# returns the error for the standard PageRank algorithm.
recordErrorsForpageRankSteps <- function(LinkMatrix, ALPHA, ERROR) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankVector <- matrix(1/numberOfPages, 1, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  errorVector <- c()
  currentError <- ERROR + 1
  numberOfIterations <- 0
  while(ERROR <= currentError) {
    tmp <- PageRankVector
    PageRankVector <- pageRankStep(PageRankVector, LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
    currentError <- sum(abs(tmp - PageRankVector))
    errorVector <- cbind(errorVector, currentError)
    numberOfIterations <- numberOfIterations + 1
  }
  
  timeTracker.finish()
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("CRUDE POWER ITERATION COMPLETE")
  return(errorVector)
}


#--------------------ADAPTIVE PAGERANK------------------------------------------

# This is a sample implementation for the use of adaptive PageRank.
#
# NOT RECOMMMENDED TO USE, since it is very slow due to the implementation
# of the matrix data type in R.
# THEORETICALLY it can outperform the normal PageRank by a lot, if 
# efficiently implemented in C++.
recordErrorAdaptivePageRank <- function(LinkMatrix, ALPHA, ERROR) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankVector <- matrix(1/numberOfPages, 1, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  convergedIndicator <- matrix(0,1, numberOfPages)
  errorVector <- c()
  currentError <- ERROR + 1
  numberOfIterations <- 0
  while(ERROR <= currentError) {
    tmp <- PageRankVector
    PageRankVector <- adaptivePageRankStep(PageRankVector, LinkMatrix, danglingPagesIndicator, convergedIndicator, numberOfPages, ALPHA)
    convergedIndicator <- convergedPagesIndicator(PageRankVector, tmp, numberOfPages, ERROR / numberOfPages)
    currentError <- sum(abs(tmp - PageRankVector))
    errorVector <- cbind(errorVector, currentError)
    numberOfIterations <- numberOfIterations + 1
  }
  
  timeTracker.finish()
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ADAPTIVE POWER ITERATION COMPLETE")
  return(errorVector)
}


#-------------------------EXTRAPOLATION ALGORITHMS-----------------------------------

#returns the errors for the PageRank algorithm with an Aitken extrapolation at a given step.
recordErrorForAitken <- function(LinkMatrix, ALPHA, ERROR, ExtrapolateAtIterationStep) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  currentError <- ERROR + 1
  
  errorVector <- c()
  
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
      errorVector <- cbind(errorVector, currentError)
      numberOfIterations <- numberOfIterations + 1
    }
  }
  
  timeTracker.finish()
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH ONE AITKEN EXTRAPOLATION COMPLETE")
  return(errorVector)
}

# records the Errors for the vector Valued Epsilon ALogrithm for k = 2.
recordErrorForEpsilonVector2 <- function(LinkMatrix, ALPHA, ERROR, ExtrapolateAtIterationStep) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  currentError <- ERROR + 1
  
  errorVector <- c()
  
  numberOfIterations <- 0
  if(ERROR <= currentError) {
    while(ERROR <= currentError) {
      if(numberOfIterations == ExtrapolateAtIterationStep) {
        PageRankSaver[3,] <- vectorEpsilon2(PageRankSaver[1,], PageRankSaver[2,], PageRankSaver[3,])
      }
      PageRankSaver[1,] <- PageRankSaver[2,]
      PageRankSaver[2,] <- PageRankSaver[3,]
      PageRankSaver[3,] <- pageRankStep(PageRankSaver[3,], LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
      currentError <- sum(abs(PageRankSaver[3,] - PageRankSaver[2,]))
      errorVector <- cbind(errorVector, currentError)
      numberOfIterations <- numberOfIterations + 1
    }
  }
  
  timeTracker.finish()
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH ONE AITKEN EXTRAPOLATION COMPLETE")
  return(errorVector)
}

# returns the errors for the PageRank algorithm with an Aitken extrapolation at a given step.
# 
# WARNING TURNS UNSTABLE REALLY FAST
recordErrorForIterativeAitken <- function(LinkMatrix, ALPHA, ERROR, nIterations) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, nIterations, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  currentError <- ERROR + 1
  
  errorVector <- c()
  
  numberOfIterations <- 0
  if(ERROR <= currentError) {
    while(ERROR <= currentError) {
      if(numberOfIterations + 1 == nIterations) {
        PageRankSaver[nIterations,] <- iterativeAitken(PageRankSaver)
      }
      for(i in 1 : (length(PageRankSaver[,1]) - 1)) {
        PageRankSaver[i,] <- PageRankSaver[i+1,]
      }
      PageRankSaver[nIterations, ] <- pageRankStep(PageRankSaver[nIterations,], LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
      currentError <- sum(abs(PageRankSaver[nIterations,] - PageRankSaver[nIterations - 1,]))
      errorVector <- cbind(errorVector, currentError)
      numberOfIterations <- numberOfIterations + 1
    }
  }
  
  timeTracker.finish()
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH ONE AITKEN EXTRAPOLATION COMPLETE")
  return(errorVector)
}


# implementation for recording the errors in the least-square PageRank.
# also displays the elpased time and prints out the sum of the extrapolated 
# coefficients. sums of the coefficients close to zero relate to 
# a better approximation of the true PageRank.
recordErrorForLeastSquare <- function(LinkMatrix, ALPHA, ERROR, vectorSaves, freqLeastSquare, ExtrapolSteps) {
  
  print("===============================")

  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, vectorSaves, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  errorVector <- c()
  
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
    errorVector <- cbind(errorVector, currentError)
    
    if((numberOfIterations  + vectorSaves + 1) %%  freqLeastSquare == 0 && numberOfIterations > 0 && numberExtrapolSteps < ExtrapolSteps) {
      PageRankSaver[vectorSaves,] <- quadraticExtrapolation(PageRankSaver)
      numberExtrapolSteps <- numberExtrapolSteps + 1
    }
    numberOfIterations <- numberOfIterations + 1
  }
  
  timeTracker.finish()
  
  print("ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH LEAST SQUARE COMPLETE")
  return(errorVector)
}


#-----------------TWO STEP SPLITTING ALGORITHMS----------------------

# computes the PageRank vector with the implementation of the TSS ALgorithm
# proposed by:
# Chuanqing Gua,Fei Xie, Ke Zhang
recordErrorsTSS <- function(LinkMatrix, ALPHA, BETA, ERROR, ERRORINNER) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * matrix(1/numberOfPages, 1, numberOfPages))
  
  Residual <- c()
  
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  currentError <- ERROR + 1
  errorSaver <- c()
  
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
    errorSaver <- cbind(errorSaver, currentError)
    
    numberOfIterations <- numberOfIterations + 1
  }
  
  PageRankSaver[2,] <- ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages) - PageRankSaver[2,]
  
  
  timeTracker.finish()
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH TWO STEP SPLITTING ITERATION COMPLETE")
  return(errorSaver)
}

# computes the PageRank vector with the implementation of the RTSS ALgorithm
# stable implementation by refining the TSS algorithm by:
# Chuanqing Gua,Fei Xie, Ke Zhang
recordErrorsRTSS <- function(LinkMatrix, ALPHA, BETA, GAMMA, ERROR, ERRORINNER) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] * sum(danglingPagesIndicator * matrix(1/numberOfPages, 1, numberOfPages))
  
  Residual <- c()
  
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  currentError <- ERROR + 1
  errorSaver <- c()
  
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
    errorSaver <- cbind(errorSaver, currentError)
    
    numberOfIterations <- numberOfIterations + 1
  }
  
  PageRankSaver[2,] <- ALPHA * PageRankSaver[3,] + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages) - PageRankSaver[2,]
  
  timeTracker.finish()
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH RELAXED TWO STEP SPLITTING ITERATION COMPLETE")
  return(errorSaver)
}
