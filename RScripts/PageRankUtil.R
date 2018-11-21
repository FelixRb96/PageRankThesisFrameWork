#
# R-file for the assisting methods in the project to clean up the code and keep things tidy.
#
# @Author: Felix Reinbott
#          Felix.Ka@web.de
#

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
  print("==============TIME ELAPSED=================")
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

