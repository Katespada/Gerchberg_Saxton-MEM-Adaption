%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
THIS IS MATLAB CODE FOR SIMULATING DIFFRACTION PATTERNS AND RECONSTRUCTING PHASE PROFILES 
WITH GERCHBERG-SAXTON ITERATIVE PHASE RETRIEVAL ALGORITHM AS WELL AS THE MODIFICATION TO INCLUDE THE MAXIMUM ENTROPY MODEL. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The following programs should be placed into the folder where you store your local MATLAB files.

a_simulate_DP.m  	%This simulates a diffraction pattern based on sample amplitude & phase
a_sample_phase.bin 	%This is the simulated sample phase profile as binary image
a_sample_amplitude.bin	%This is the simulated sample amplitude profile as binary image 
a_dp.bin		%This is the diffraction pattern as a result of running a_simulate_DP
a_dp.jpg		%This is the image of the diffraction pattern as a result of a_simulate_DP
FT2Dc.m			%Performs the 2D centered Fourier transform
IFT2Dc.m		%Performs the 2D centered inverse Fourier transform
GS.m			%The Gerchberg-Saxton Iterative Phase Retrieval Algorithm 
MEM.m			%The adapted algorithm for inclusion of the Maximum Entropy Model 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gerchberg-Saxton algorithm recovers phase of the wavefront from two intensity measurements 
in the sample plane and in the detector plane. 

The sample is a complex-valued distribution; the amplitude is a binary cat cartoon 
with maximum of 1, and the phase is the word "cat" with maximum of 0.1 radian. 

The programs should be run as follows:
"a_simulate_DP.m" simulates the diffraction pattern, this step can be ignored with a_dp.bin saved.
"GS.m" iteratively reconstructs the complex-valued wavefront with phase based on Gerchberg-Saxton.
"MEM.m' iteratively reconstructs the complex-valued wavefront with phase based on Maximum Entropy. 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Citation for parts of this algorithm: Tatiana Latychevskaia
"Iterative phase retrieval in coherent diffractive imaging: practical issues",
Applied Optics 57(25), 7187 - 7197 (2018)

Katelyn Spadavecchia 
Katelyn_Spadavecchia@mines.edu
