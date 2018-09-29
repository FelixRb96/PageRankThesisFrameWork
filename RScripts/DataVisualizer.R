#
# method provider for analyzing and visualizing the PageRank algorithms 
# properties, convergence speed, and data presentation for different ALPHAS.
# 
# Automatically loads the needed R - files for executing the framework, and some example data
# is provided in the ./BachelorThesisData/Data folder, for recreating the numerical results
# shown in the actual thesis. 
# 
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#---------------------------------------------PerformanceReporter-----------------------------------------------------------------
#
# Records the Errors for the different algorithms on a LinkMatrix and
# creates a plot to compare the errors on a logarithmic scale.
#
# Functions are given the default parameters from the thesis.
#
performanceReport <- function(LinkMatrix, alpha, beta, gamma, ERROR, errorInner) {
  
  source('~/BachelorThesisData/RScripts/PageRankErrorComputation.R')
  
  print("STARTING CALCULATION...")
  DataSaver <- recordErrorsForpageRankSteps(LinkMatrix, alpha, ERROR)
  ErrorSaverAitken <- recordErrorForAitken(LinkMatrix, alpha, ERROR, 6)
  ErrorSaverLeastSquare <- recordErrorForLeastSquare(LinkMatrix, alpha, ERROR, 4, 10, 2)
  ErrorSaverTSS <- recordErrorForTSS(LinkMatrix, alpha, beta, ERROR, errorInner)
  ErrorSaverRTSS <- recordErrorForRTSS(LinkMatrix, alpha, beta, gamma, ERROR, errorInner)
  
  DataSaver <- fitDataForPlotting(DataSaver, ErrorSaverAitken)
  DataSaver <- fitDataForPlotting(DataSaver, ErrorSaverLeastSquare)
  DataSaver <- fitDataForPlotting(DataSaver, ErrorSaverTSS)
  DataSaver <- fitDataForPlotting(DataSaver, ErrorSaverRTSS)
  x <- length(DataSaver[1,])
  
  print("FINISHED CALCULATION...")
  print("PRINTING PLOT...")
  plot(1:x, DataSaver[1,], col = "red", log = "y", type = "l")
  points(1:x, DataSaver[2,], col = "blue", type = "l")
  points(1:x, DataSaver[3,], col = "green", type = "l")
  points(1:x, DataSaver[4,], col = "black", type = "l")
  points(1:x, DataSaver[5,], col = "turquoise", type = "l")
}

performanceReportForAlphas <- function(LinkMatrix, ERROR) {
  
  source('~/BachelorThesisData/RScripts/PageRankErrorComputation.R')
  
  ALPHAS <- 1:9 / 10
  ErrorsCrude <- c()
  ErrorsAitken <- c()
  ErrorsLS <- c()
  for(i in 1:length(ALPHAS)) {
    if(i > 1) {
      ErrorsCrudeCurrent <- recordErrorsForpageRankSteps(LinkMatrix, ALPHAS[i], ERROR)
      ErrorsCrude <- fitDataForPlotting(ErrorsCrude, ErrorsCrudeCurrent)
      
      ErrorsAitkenCurrent <- recordErrorForAitken(LinkMatrix, ALPHAS[i], ERROR, 5)
      ErrorsAitken <- fitDataForPlotting(ErrorsAitken, ErrorsAitkenCurrent)
      
      ErrorsLSCurrent <- recordErrorForLeastSquare(LinkMatrix, ALPHAS[i], ERROR, 4, 10, 2)
      ErrorsLS <- fitDataForPlotting(ErrorsLS, ErrorsLSCurrent)
      
    } else {
      ErrorsCrude <- recordErrorsForpageRankSteps(LinkMatrix, ALPHAS[i], ERROR)
      ErrorsAitken <- recordErrorForAitken(LinkMatrix, ALPHAS[i], ERROR, 10)
      ErrorsLS <- recordErrorForLeastSquare(LinkMatrix, ALPHAS[i], ERROR, 4, 10)
    }
  }
  
  ErrorLine <- matrix(ERROR, 1, lengthOverall)
  
  for(i in 1:length(ALPHAS)) {
    
    DataSaver <- ErrorsCrude[i,]
    DataSaver <- fitDataForPlotting(DataSaver,ErrorsAitken[i,])
    DataSaver <- fitDataForPlotting(DataSaver, ErrorsLS[i,])
    
    plot(1:lengthOverall,DataSaver[1,], col = "red", log = "y", type = "l")
    points(1:lengthOverall, DataSaver[2,], col = "blue", type = "l")
    points(1:lengthOverall, DataSaver[3,], col = "green ", type = "l")
    points(1:lengthOverall, ErrorLine, col = "black", type = "l")
  }
}

#---------------------------visualization PageRank vector.---------------------------------

# visualizes the changes in the pagerank PageRank vector for different ALPHAS.
# 
# this method is set up with default parameters and is for the sole purpose of visualizing
# the right convergence across all implementations.
#
PageRankChangeReportForAlpha <-function(LinkMatrix, ERROR, numberOfPages){
  
  source('~/BachelorThesisData/RScripts/PageRankComputation.R')
  
  ALPHAS <- 2:9 / 10
  len <- length(LinkMatrix[1,])
  for(i in 1:8){
    PageRankCrude <- computePageRankByError(LinkMatrix, ALPHAS[i], ERROR)
    PageRankAitken <- computePageRankAitken(LinkMatrix, ALPHAS[i], ERROR, 5)
    PageRankLS <- computePageRankLeastSquare(LinkMatrix, ALPHAS[i], ERROR, 4, 10)
    PageRankTSS <- computePageRankTSS(LinkMatrix, ALPHAS[i], ALPHAS[i] + 0.05, ERROR, 0.97)
    # PageRankRTSS <- computePageRankRTSS(LinkMatrix, ALPHAS[i], ALPHAS[i] + 0.05, ALPHAS[i] - 0.1, ERROR, 1.5)
    
    sortedPageRankCrude <- sort(PageRankCrude, decreasing = TRUE)[1:numberOfPages]
    sortedPageRankAitken <- sort(PageRankAitken, decreasing = TRUE)[1:numberOfPages]
    sortedPageRankLS <- sort(PageRankLS, decreasing = TRUE)[1:numberOfPages]
    sortedPageRankTSS <- sort(PageRankTSS, decreasing = TRUE)[1:numberOfPages]
    # sortedPageRankTRSS <- sort(PageRankRTSS)[1:numberOfPages]
    
    print("==================CRUDE====================")
    print(which(PageRankCrude %in% sortedPageRankCrude))
    print(sortedPageRankCrude)
    print("==================AITKEN====================")
    print(which(PageRankAitken %in% sortedPageRankAitken))
    print(sortedPageRankAitken)
    print("==================LEAST SQUARE====================")
    print(which(PageRankLS %in% sortedPageRankLS))
    print(sortedPageRankLS)
    print("==================TSS====================")
    print(which(PageRankTSS %in% sortedPageRankTSS))
    print(sortedPageRankTSS)
    print("==================RTSS====================")
    # print(which(PageRankRTSS %in% sortedPageRankRTSS))
    # print(sortedPageRankRTSS)
  }
}

#-------------------------------DataConverter------------------------------------

# calls the DataConverter to set up a LinkMatrix.
setUpMatrix <- function(Data, MatrixSize) {
  source('~/BachelorThesisData/RScripts/DataConverter.R')
  returnMatrix <- transformDataToSubWebLinkMatrix(Data, MatrixSize)
  linkMatrixStats(returnMatrix)
  return(returnMatrix)
}


