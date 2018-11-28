#
# Tests the implemented methods and can/should be executed after changing code, 
# to see if all functions in the framework can deal with the canges.
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#--------------------------------------------------------------------

testMatrix <- function() {
  
  source('~/BachelorThesisData/RScripts/DataConverter.R')
  
  web.NotreDame.txt <- read.delim("~/BachelorThesisData/Data/web-NotreDame.txt.gz", header=FALSE, comment.char="#")
  
  matrix <- transformDataToSubWebLinkMatrix(web.NotreDame.txt, 5000)
  linkMatrixStats(matrix)
  
  return(matrix)

}

test_PageRankComputation <- function(m = testMatrix()) {
  
  source('~/BachelorThesisData/RScripts/PageRankComputation.R')
  
  err <- 10^(-3)
  al <- 0.8
  
  computePageRankByError(m ,al, err)
  computeAdaptivePageRankByError(m, al, err)
  computePageRankAitken(m ,al , err, 5)
  computePageRankLeastSquare(m, al, err, 4, 8, 2)
  computePageRankTSS(m, al, 0.5, err, 0.01)
  computePageRankRTSS(m, al, 0.5, 0.9, err, 0.01)
  computePageRankByIterations(m, al, 10)
  
}

test_ErrorComputation <- function(m = testMatrix()) {
  
  source('~/BachelorThesisData/RScripts/PageRankErrorComputation.R')
  
  err <- 10^(-3)
  al <- 0.8
  
  recordErrorsForpageRankSteps(m ,al, err)
  recordErrorAdaptivePageRank(m ,al, err)
  recordErrorForAitken(m ,al, err, 5)
  recordErrorForIterativeAitken(m ,al, err, 5)
  recordErrorForEpsilonVector2(m ,al, err, 5)
  recordErrorForLeastSquare(m ,al, err, 4, 8, 2)
  recordErrorsTSS(m ,al , 0.5, err, 0.01)
  recordErrorsRTSS(m ,al , 0.5, 0.9, err, 0.01)
  
}

test_all <- function() {
  test_PageRankComputation()
  test_ErrorComputation()
  
  print("ALL TESTS COMPLETE")
}