#
# R-file fot the implementation of the Gauss-Seidel Algorithm to compute
# a sparse linear system approximatively.
#
# Proofs to the algorithms are not provided, and this file 
# is basically just there for the comparison to the best applicable
# LES solver for the problem proposed in:
# 
# "Fast PageRank Computation via a Sparse Linear System" 
# by Gianna M. Del Corso, Antonio Gullí, and Francesco Romani
#
# It might be added, that this method did (at least due to constraints in R)
# not perform well compared to Iterative computations
# on a primitive test for PageRank, although it might be an 
# interesting subject with permutations applied to the LinkMatrix.
# 
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#--------------------------------------------------------------------

# implementation of the Gauss-Seidel in a compact form. 
# @param x is the initial value for the Iteration, bad choices may 
# lead to a long convergence time.
# 
# If used for the PageRank-Problem use the following inputs:
# @param A = setUpIdMatrix(numberOfPages) - ALPHA * t(LinkMatrix)
# @param x = matrix(1/numberofPages, numberOfPages, 1)
# @param b = (1 - ALPHA) * matrix(1/numberofPages, numberOfPages, 1)
#
gauss_seidel <- function(A, x, b, ERROR, maxItrerations) {
  
  source('~/BachelorThesisData/RScripts/PageRankUtil.R')
  
  timeTracker.start()
  
  currentError <- ERROR + 1
  numberOfPages <- length(b)
  nIterations <- 0
  while(ERROR <= currentError && nIterations < maxItrerations) {
    tmp <- x
    x[1] <- 1 / A[1,1] * (b[1] - sum(A[1, 2:numberOfPages] *  x[2:numberOfPages]))
    for(i in 2: (numberOfPages - 1)) {
      x[i] <- 1/A[i,i] * (b[i] - sum(A[i,(1:i-1)] * x[1:(i-1)]) - sum(A[i, (i+1):numberOfPages] * x[(i+1):numberOfPages]))
    }
    x[numberOfPages] <- 1 / A[numberOfPages, numberOfPages] * (b[numberOfPages] - sum(A[numberOfPages, 1:(numberOfPages - 1)] * x[1:(numberOfPages - 1)]))
    currentError <- sum(abs(x - tmp))
    print(currentError)
    nIterations <- nIterations + 1
  }

  timeTracker.finish()
  print("ITERATION STEPS:")
  print(nIterations)
  
  # if not used for PageRank comment this out !
  x <- x / sum(abs(x))
  
  return(x)
}


# sets up a identity matrix with the 
# dimension n x n. 
# Of course only intger choices for
# n make sense. 
setUpIdMatrix <- function(n) {
  mat <- matrix(0, n, n)
  for(i in 1:n) {
    mat[i,i] <- 1
  }
  return(mat)
}
