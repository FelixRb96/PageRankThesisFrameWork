#
# R-file for the assisting methods in the project to clean up the code and keep things tidy.
# contains the implementations for the extrapolation steps and the PageRank step.
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#-----------------------------------------------------------------------------------------
# basic PageRank Iteration.

# assuming xn has dimension (m x n) and is normalized. As well as LinkMatrix has dimension (n x n) this returns
# the effective implemented way of a PageRank iteration, with the personalization vector 1/n * e.
# returns the normed value
pageRankStep <- function(xn, LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA) {
  iteratedStepAlpha <- ALPHA * (xn %*% LinkMatrix + xn %*% danglingPagesIndicator %*% matrix(1/numberOfPages, 1, numberOfPages))
  iteratedStepLeftover <- (1 - ALPHA) * matrix(1/numberOfPages, 1, numberOfPages)
  iteratedStep <- iteratedStepAlpha + iteratedStepLeftover
  iteratedStep <- iteratedStep / sum(iteratedStep)
  return(iteratedStep)
}

#----------------------------------------------------------------------------------------
# the implementation of the Aitken delta Squared

# method for the Aitken extrapolation step with three iterates in ascending order (of course x1 is the first x3 the last iterate)
aitkenExtrapolation <- function(x1, x2, x3) {
  extrapolatedVector <- x3 - ((x3 - x2)^2) / (x3 - 2 * x2 + x1)
  print(sum(abs(extrapolatedVector)))
  return(extrapolatedVector)
}

#-----------------------------------------------------------------------------------------
# Implementation of the least square extrapolation step and helping methods

# computes the extrapolated vector from thwe approximated minimal polynomial, with the coefficients obtained by
# least square. 
# WARNING: this function contains an inversion an will get increasingly expensive for a larger number of saved vectors, 
# thus keeping the number of stored iterations reasonably small is advised.
quadraticExtrapolation <- function(SavedIterations) {
  numberOfSaves <- length(SavedIterations[,1])
  dataLength <- length(SavedIterations[1,])
  Differences <- matrix(0, numberOfSaves - 1, dataLength)
  for(i in 2 : numberOfSaves) {
    Differences[i-1,] <- SavedIterations[i,] - SavedIterations[1,]
  }
  lsMatrix <- t(Differences[1:(numberOfSaves - 2), ]) %*% qr.solve(Differences[1:(numberOfSaves -2),] %*% t(Differences[1:(numberOfSaves - 2),]))
  params <- Differences[numberOfSaves - 1,] %*% lsMatrix
  params <- cbind(params, 1)
  
  print("SUM OF COEFFICIENTS IN EXTRAPOLATION STEP:")
  print(sum(params))

  
  extrapolatedVector <- params %*% lowerTriangular(1, numberOfSaves - 1, numberOfSaves - 1) %*% SavedIterations[2:numberOfSaves,]
  return(extrapolatedVector)
}

# A function to construct a lower triangular matrix on size (m x n) filled with the same entry.
lowerTriangular <- function(numberToFill, m, n) {
  result <- matrix(numberToFill, m, n)
  for(i in 1 : (m - 1)) {
    for(j in (i + 1) : n) {
      result[i,j] <- 0
    }
  }
  return(result)
}

#-------------------------------------------------------------------
# implementation of the Arnoldi and supporting methods

# returns an orthogonal base of the induced Krylov subspace by 
# the transformed LinkMatrix 
arnoldi <- function(LinkMatrix, ALPHA, startVec, Dim) {
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  nPages <- length(startVec)
  resultBase <- rbind(startVec, matrix(0, Dim, nPages))
  represLM <- matrix(0,Dim + 1, Dim + 1)
  for(j in 1:Dim) {
    tmp <- pageRankStep(resultBase[j,], LinkMatrix, danglingPagesIndicator, nPages, ALPHA)
    print(tmp[,1:10])
    for(i in 1:j) {
      orthogonalizer <- sum(tmp * resultBase[i,])
      tmp <- tmp - orthogonalizer * resultBase[i,]
      represLM[i,j] <- orthogonalizer
    }
    orthogonalizer <- norm2(tmp)
    represLM[i + 1,j] <- orthogonalizer
    if(orthogonalizer == 0) {
      stop("INVALID OPERATION ERROR")
    }
    resultBase[j + 1, ] <- tmp/orthogonalizer
  }
  setClass(Class = "arnoldiReturnObj", representation(Base = "matrix", SubSpacePres = "matrix"))
  result <- new("arnoldiReturnObj" , Base = resultBase, SubSpacePres = represLM)
  return(result)
}

# simple implementation of the euklidean norm
# for easier code purposes.
norm2 <- function(x) {
  returnVal <- sqrt(sum(x^2))
  return(returnVal) 
}


#--------------------------------------------------------------------
# supporting methods for plotting and performance tracking

# starts the timer for timercording
timeTracker.start <- function() {
 
  if("devtools" %in% rownames(installed.packages()) == FALSE) {
    install.packages("devtools")
  }
  if("tictoc" %in% rownames(installed.packages()) == FALSE) {
    install.packages("tictoc")
  }
  tictoc::tic()
}

# ends the timer and prints the elapsed time.
timeTracker.finish <- function() {
  print("===============================")
  tictoc::toc()
}

# binds together the data of the Errors that are already saved and
# adds a new vector of errors.
fitDataForPlotting <- function(DataSaver, DataToAdd) {
  if(length(DataSaver[1,]) < length(DataToAdd)) {
    lastCOlVector <- matrix(DataSaver[,length(DataSaver[1,])], length(DataSaver[,1]), 1)
    DataSaver <- cbind(DataSaver, matrix(lastCOlVector, length(lastCOlVector), length(DataToAdd) - length(DataSaver[1,])))
    DataSaver <- rbind(DataSaver, DataToAdd)
  } else {
    lastEntry <- DataToAdd[length(DataToAdd)]
    DataToAdd <- cbind(DataToAdd, matrix(lastEntry, 1, length(DataSaver[1,]) - length(DataToAdd)))
    DataSaver <- rbind(DataSaver, DataToAdd)
  }
  return(DataSaver)
}

