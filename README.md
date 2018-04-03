# Melodic_Template
*** COMPUTING MELODIC TEMPLATES IN ORAL MUSIC TRADITIONS ***

+++ ABOUT +++

This repository contains software and data to reproduce the results reported in the publication

S. Bereg, J.-M. Díaz-Báñez, N. Kroher and I. Ventura (2017). Computing melodic templates in flamenco singing. XVII Spanish Meeting on Computational Geometry.

If you use this code in your work, please cite the publication above. The software is provided by the authors for research purposes only, in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

rpiedra at us dot es

+++ DEPENDENCIES +++

* MATLAB (https://es.mathworks.com/products/matlab.html)
* Generate maximally perceptually-distinct colors (https://es.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
* Fast Line Segment Intersection ( https://es.mathworks.com/matlabcentral/fileexchange/27205-fast-line-segment-intersection)

+++ USAGE +++

To use these functions is necessary to download MATLAB 

The problem that you want to solve this software is given n linear functions to pieces find a 2 * epsilon amplitude tube so that at each instant it contains p functions (although it does not have to always be the same functions).

* Decision_Problem3: returns 1 if it is possible to build the tube or 0 if not.
	Input values: 
	* F are n piecewise functions, F is a cell of n matrix where each matrix has by rows each segment of the piecewise function 
  ( [a b c d; c d e f ..], a,b,c,d,e,f in R) 
	* epsilon is half the tube radius
	* p is the number of the functions that have to contain the tube in every time.

* Binary_search_epsilon: return the minimum epsilon for which is possible to construct the tube.
	Input values:
	* F
	* p 

* Tube_Problem: returns by screen, if it is possible to build it, a picture with the functions F, the tube of amplitude 2*epsilon and the function center the tube.
	Input values:
	* F
	* epsilon
	* p
