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

#returns the errors for the PageRank algorithm with an Aitken extrapolation at a given step.
recordErrorForAitken <- function(LinkMatrix, ALPHA, ERROR, ExtrapolateAtIterationStep) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
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

recordErrorForLeastSquare <- function(LinkMatrix, ALPHA, ERROR, vectorSaves, freqLeastSquare, ExtrapolSteps) {
  
  print("===============================")

  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
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

# 
recordErrorForTSS <- function(LinkMatrix, ALPHA, BETA, ERROR, ErrorInner) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[1,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  errorVector <- c()
  currentError <- ERROR + 1
  
 while(ERROR <= currentError) {
   PageRankSaver[1,] <- PageRankSaver[2,]
   reverseSplitting <- 
  Residual <- ErrorInner + 1
   while(ErrorInner <= sum(abs(Residual - PageRankSaver[2,]))) {
     PageRankSaver[2,] <- (ALPHA * PageRankSaver[3,] 
     + ALPHA * PageRankSaver[3,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
     + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages))
     PageRankSaver[2,] <- PageRankSaver[2,] / sum(abs(PageRankSaver[2,]))
     PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
     PageRankSaver[3,] <- PageRankSaver[3,] / sum(abs(PageRankSaver[3,]))
     Residual <- (ALPHA - BETA) * PageRankSaver[3,] + (1 - ALPHA) * matrix(1 / numberOfPages, 1, numberOfPages)
     numberOfInnerIterations <- numberOfInnerIterations + 1
     print(sum(abs(Residual - PageRankSaver[2,])))
   }
   PageRankSaver[2,] <- Residual + BETA * PageRankSaver[3,]
   PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
   currentError <- sum(abs(PageRankSaver[2,] - PageRankSaver[1,]))
   print(currentError)
   
   errorVector <- cbind(errorVector, currentError)
   numberOfIterations <- numberOfIterations + 1
 }
  
  timeTracker.finish()
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH TWO STEP SPLITTING ITERATION COMPLETE")
  return(errorVector)
}

# 
recordErrorForRTSS <- function(LinkMatrix, ALPHA, BETA, GAMMA, ERROR, ErrorInner) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  timeTracker.start()
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankSaver <- matrix(1/numberOfPages, 3, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  
  PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[1,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
  numberOfIterations <- 0
  numberOfInnerIterations <- 0
  
  errorVector <- c()
  currentError <- ERROR + 1
  
  while(ERROR <= currentError) {
    PageRankSaver[1,] <- PageRankSaver[2,]
      Residual <- ErrorInner + 1
    while(ErrorInner <= sum(abs(Residual - PageRankSaver[2,]))) {
      PageRankSaver[2,] <- (ALPHA * PageRankSaver[3,] 
      + ALPHA * PageRankSaver[3,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
      + (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages))
      PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
      PageRankSaver[3,] <- PageRankSaver[3,] / sum(abs(PageRankSaver[3,]))
      Residual <- ((GAMMA - 1) / GAMMA * PageRankSaver[2,] 
      + (ALPHA - BETA) / GAMMA * PageRankSaver[3,] 
      + (1 - ALPHA) / GAMMA * matrix(1/numberOfPages, 1, numberOfPages))
      numberOfInnerIterations <- numberOfInnerIterations + 1
    }
    PageRankSaver[2,] <- Residual + BETA / GAMMA * PageRankSaver[3,]
    PageRankSaver[3,] <- PageRankSaver[2,] %*% LinkMatrix + PageRankSaver[2,] %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages)
    currentError <- sum(abs(PageRankSaver[2,] - PageRankSaver[1,]))
    
    errorVector <- cbind(errorVector, currentError)
    numberOfIterations <- numberOfIterations + 1
  }
  
  timeTracker.finish()
  
  print("INNER ITERATION STEPS:")
  print(numberOfInnerIterations)
  print("OUTER ITERATION STEPS: ")
  print(numberOfIterations)
  print("ITERATION WITH RELAXED TWO STEP SPLITTING ITERATION COMPLETE")
  return(errorVector)
}

recordErrorArnoldiExtrapoaltion <- function(LinkMatrix, ALPHA, ERROR, BETA, restart, tolNorm){
  # TODO
}

