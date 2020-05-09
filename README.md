# PR-4Pi
PR-4Pi software package for analyzing single molecule blinking dataset from 4Pi-SMSN systems.

# Description
The demo package consists of functions and scripts written in MATLAB (MathWorks, Natick, MA). The code has been tested in MATLAB version R2016b. The current version is only supported on Windows platform with GPU equipped. 

# Required packages
DIP image toolbox (http://www.diplib.org/download) 

NVIDIA GPU Computing Toolkit v7.5 (https://developer.nvidia.com/cuda-75-downloads-archive)

Parallel Computing Toolbox in MATLAB

# Content of the software package
	SR4pi_demo 	-	Matlab class for analyzing 4Pi-SMSN data
	PSF Toolbox	-	A set of Matlab classes for phase retrieval (see manual for PSF Toolbox)
	mex		-	mex-functions used in the software package
	source code	-	c and cuda source code for major mex-functions
	test data	-	test data for the demo script of the software
	SR4pi_demo_example.m	-	demo script of the software
	
# How to run
	1. Change current folder in Matlab to PR-4Pi software.
	2. Run demo script SR4pi_demo_example.m. (Note: there are a few steps requiring user interactions, please see the demo script for details.)
	3. Type ‘help SR4pi_demo’ in Matlab command window for detailed help on SR4pi_demo class.
	
## Citation
Please cite PR-4Pi in your publications if it helps your research:

  

