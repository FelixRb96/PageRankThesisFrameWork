#
# Same implementations of the PageRank algorithms as in @PageRankComputation and 
# @PageRankErrorComputation.
# 
# But these implementations provide special measurements to further visualize 
# PageRank behaviour for different implementations.
#
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#---------------------------------------------------------------------------------

# returns the error for the standard PageRank algorithm.
recordPagesConverged <- function(LinkMatrix, ALPHA, ERROR) {
  
  source('~/BachelorThesisData/RScripts/PageRankSupportLib.R')
  
  numberOfPages <- length(LinkMatrix[,1])
  PageRankVector <- matrix(1/numberOfPages, 1, numberOfPages)
  danglingPagesIndicator <- abs(rowSums(LinkMatrix) - 1)
  pagesConvergedVector <- c()
  currentError <- ERROR + 1
  pagesConverged <- 0
  
  while(ERROR <= currentError) {
    tmp <- PageRankVector
    PageRankVector <- pageRankStep(PageRankVector, LinkMatrix, danglingPagesIndicator, numberOfPages, ALPHA)
    currentError <- sum(abs(tmp - PageRankVector))
    
    pagesConverged <- sum(abs(PageRankVector - tmp) <  ERROR / numberOfPages)
    pagesConvergedVector <- cbind(pagesConvergedVector, pagesConverged)
  }
  
  return(pagesConvergedVector)
}