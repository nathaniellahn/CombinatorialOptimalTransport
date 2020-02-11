# Combinatorial Optimal Transport Implementation

This repository contains the MATLAB code for the optimal transport algorithm described in:
https://arxiv.org/abs/1905.11830
The code is a joint effort between the authors of that paper. 

Note: Paper has been accepted into NeurIPS 2019.

## Requirements
1. To use this code, obtain MATLAB from https://www.mathworks.com/products/matlab.html
2. Download or clone the git project.
3. The following test files can be used to reproduce results from the paper:
	a.)SmallDelta.m : Test our algorithm on small delta values.
	b.)SinkhornComparison.m : Compare performance to the Sinkhorn algorithm. 
	c.)IterationsAllMethods.m : Compare iteration count to Sinkhorn, Greenkhorn, and APDAGD
See the paper for further details.


The tests are written in MATLAB and call a compiled Java implementation of our algorithm. The Java binary files 
for Mapping.java and GTTransport.java will need to be created and placed in the 
GabowTarjanJavaCode\GTTransport\bin\optimaltransport folder. Here, the optimaltransport directory corresponds to a
Java package that includes both files. We used Java 8, but most likely other Java versions will do just fine.

VERY IMPORTANT: When running the MATLAB test files, make sure to run them from within the Testing folder. 
The testing files use relative paths, and executing with a different path location for your script will affect your results 
(most likely, the Java classes won't be found). If you are not comfortable with this, try replacing the command:
	javaaddpath('..\GabowTarjanJavaCode\GTTransport\bin\');
with an absolute path on your system. You can use the "javaclasspath" command to ensure that the path to the Java class files 
is correctly entered.

The APDAGD implementation and some testing code was obtained from https://github.com/chervud/AGD-vs-Sinkhorn.
The Sinkhorn and Greenkhorn implementations, the file computeot_lp.m and a few other helper functions were 
obtained from https://github.com/JasonAltschuler/OptimalTransportNIPS17/
These files are used to compare our implementation against existing algorithms.
All credit for those files goes to the respective authors. 


Note: The first iteration of this implementation was written in MATLAB and has been preserved in this repository. 
However, it is signiciantly slower than the Java implementation and has not been 'cleaned up' as much.
It is highly recommended to use the Java version, especially since the Java code can be called from MATLAB.
While MATLAB is great for Sinkhorn, etc, Java is a better choice for our algorithm.