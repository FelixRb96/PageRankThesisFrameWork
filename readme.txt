===============================================
@Title of Project:
	Google: PageRank and accelerated methods
@author:
	Felix Reinbott, <Felix.Ka@web.de>
@Date Due:
	03.12.2018
===============================================
Short summary:

This project consisting of R-files and some sample Data
is used to display performance and properties for the 
algorithms presented and proven in the thesis.
===============================================
Requirementents:

So far there is no version control needed, so a 
(reasonably actual) distribution of R and for the ease
of use R-Studio should work.
===============================================
Use:
-------------------------
I: IMPORT SAMPLE DATA.

Sample Data can be accessed by calling the following methods in an
R console:

BerkeleyStanfordWebData <- read.delim("~/BachelorThesisData/Data/BerkeleyStanfordWebData.txt", header=FALSE, comment.char="#")
GoogleWebData2002 <- read.delim("~/BachelorThesisData/Data/GoogleWebData2002.txt", header=FALSE, comment.char="#")

calling those methods will provide the Data used to display the implementation results in the thesis.

DISCLAIMER:
I DO NOT OWN THE DATA AND ALL CREDITS GO TO THE RESPECTIVE OWNERS:

THE GOOGLE DATASET IS FROM:
J. Leskovec, K. Lang, A. Dasgupta, M. Mahoney. 
Community Structure in Large Networks: 
Natural Cluster Sizes and the Absence of Large Well-Defined Clusters. 
Internet Mathematics 6(1) 29--123, 2009.
Google programming contest, 2002

THE BERKELEY-STANFORD DATA IS FROM:

J. Leskovec, K. Lang, A. Dasgupta, M. Mahoney. 
Community Structure in Large Networks: 
Natural Cluster Sizes and the Absence of Large Well-Defined Clusters. 
Internet Mathematics 6(1) 29--123, 2009.

---------------------------------
II: Load the framework:

sourcing the DataVisualizer R-file loads all needed methods
automatically.

source('~/BachelorThesisData/RScripts/DataVisualizer.R')

All methods to get "interesting" results for comparing the algorithms
are executed by calling the respective methods and for further
informations check the methods documentation.

After calling the source(...) method it is necessary to load a LinkMatrix
by calling:

setUpMatrix(Data, MatrixSize)

and using that matrix for the analyzing methods.
============================================================
For Questions or Ideas please send an E-mail to the provoded 
E-mail adress.

--Felix Reinbott





