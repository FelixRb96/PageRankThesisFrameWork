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

Download Data from:

https://snap.stanford.edu/data/index.html
(e.g. web-NotreDame)

and Save in the ./Data folder. 

Sample Data can be accessed by calling the read.delim methods in an
R console with the following settings:

web.NotreDame.txt <- read.delim("~/BachelorThesisData/Data/web-NotreDame.txt.gz", header=FALSE, comment.char="#")
(here example for wb-NotreDame)

calling those methods will provide the Data used to display the implementation results in the thesis.



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
For Questions or Ideas please send an E-mail to the provided 
E-mail adress.

--Felix Reinbott





