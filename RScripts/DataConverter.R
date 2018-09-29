#
# Library for transforming web hyperlink structure data retrieved from *.txt documents to
# link matrices, and general information about the retrieved linkmatrix.
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

#----------------------------------------------------------------------------------------------------------
# Function to transform the given Data Frames to a link matrix, where m[i,j] represents, if
# page i links to page j. The given LinkMatrix is row stochastic so left side multiplication has
# to be applied to compute the PageRank.
#
transformDataToSubWebLinkMatrix <- function(x, subWebSize) {
  dataLength <- length(x[,"V1"])
  resultMatrix <- matrix(0, subWebSize, subWebSize)
  for (i in  1:dataLength){
    if(x[i, "V1"] <= subWebSize && x[i, "V2"]<= subWebSize){
      resultMatrix[x[i, "V1"], x[i, "V2"]] <- 1
      print("FOUND ENTRY... Continuing")
    }
    if(i == dataLength) {
      print("ADJACENCY MATRIX CONSTRUCTED")
    }
  }
  rowSumsLinkMatrix <- rowSums(resultMatrix)
  for(i in 1:subWebSize){
    if(rowSumsLinkMatrix[i] != 0) {
      resultMatrix[i,] <- resultMatrix[i,] / rowSumsLinkMatrix[i]
    }
  }
  print("LINK MATRIX CONSTRUCTED")
  return(resultMatrix)
}

# Same as @transformDataToSubWebMatrix, but this time with the selection of random pages from the web.
# WARNING: is very very slow and usually only gives a matrix of nearly only zeros.
transformDataToRandomSubWebLinkMatrix <- function(x, subWebSize) {
  dataLength <- length(x[,"V1"])
  RandomPageIndex <- sort(round(runif(subWebSize, 0.5, dataLength + 0.5)))
  print("RANDOM PAGES SELECTED...")
  resultMatrix <- matrix(0, subWebSize, subWebSize)
  for(i in 1:subWebSize) {
    for(j in 1:dataLength) {
      if(x[j,"V1"] == RandomPageIndex[i]) {
        for(k in 1:subWebSize) {
          if(x[j,"V2"] == RandomPageIndex[k]) {
            resultMatrix[i,k] <- 1
            print("ENTRY FOUND CONTINUING...")
          }
        }
      }
    }
  }
  print("ADJACENCY MATRIX CONSTRUCTED")
  rowSumsLinkMatrix <- rowSums(resultMatrix)
  for(i in 1:subWebSize){
    if(rowSumsLinkMatrix[i] != 0) {
      resultMatrix[i,] <- resultMatrix[i,] / rowSumsLinkMatrix[i]
    }
  }
  print("LINK MATRIX CONSTRUCTED")
  return(resultMatrix)
}

#------------------------------------------------------------------

# creates an overview over the properties of the conatructed LinkMatrix.
linkMatrixStats <- function(LinkMatrix) {
  
  numberOfWebPages <- length(LinkMatrix[,1])
  outlinksPerPage <- rowSums(LinkMatrix)
  totalOutLinks <- sum(outlinksPerPage)
  averageOutlinks <- totalOutLinks / numberOfWebPages
  
  danglingPages <- 0
  for(i in 1:numberOfWebPages) {
    if(outlinksPerPage[i] == 0) {
      danglingPages <- danglingPages + 1
    }
  }
  
  print("OUTPUT")
  print("=======================")
  print("NUMBER OF WEB PAGES:")
  print(numberOfWebPages)
  print("=======================")
  print("TOTAL OUTLINKS:")
  print(totalOutLinks)
  print("=======================")
  print("AVERAGE OUTLINKS PER PAGE:")
  print(averageOutlinks)
  print("=======================")
  print("NUMBER OF DANGLING PAGES")
  print(danglingPages)
}
