# Gabow Tarjan Implementation

This repository contains the MATLAB code for the Gabow-Tarjan optimal transport algorithm described in:
https://arxiv.org/abs/1905.11830

## Requirements
1. To use this code, download matlab from https://www.mathworks.com/products/matlab.html
2. Download or clone the git project.
3. The following test files demonstrate usage.
	a. VaryingC.m : Execute our algorithm on randomly generated supplies, demands, and costs. Costs are drawn uniformly at random from the interval [0,C], for a range of different C values.
	b. LinprogCompare.m : A comparison of our algorithm with MATLAB's linear programming solver, https://www.mathworks.com/help/optim/ug/linprog.html.
	c. SinkhornComparison.m : A comparison of our algorithm with the Sinkhorn implementation from https://github.com/chervud/AGD-vs-Sinkhorn. Their implementation is based on Algorithm 3 of the paper https://papers.nips.cc/paper/6792-near-linear-time-approximation-algorithms-for-optimal-transport-via-sinkhorn-iteration.pdf
4. Note: Relevant files will need to be added to the MATLAB path. To do this, right click the folders and select "Add to Path".

The folder Sinkhorn, and some testing code was obtained from https://github.com/chervud/AGD-vs-Sinkhorn.
The file computeot_lp.m was obtained from https://github.com/JasonAltschuler/OptimalTransportNIPS17/
These files are used to compare our implementation against existing algorithms.
All credit for those files goes to the respective authors. 

