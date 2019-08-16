!The files in Cantilever_plate folder explain preprocessing and analysis of 
Cantilever beam elstic problem using proposed Abaqus NSFEM UEL
----------------------------------------------------------------------------
# Preprocessing task
Node_AdjNodes_code - MATLAB code for finding attached nodes for the node
under consideration and this also generates correspong UEL files.
UEL definitions are written in seperate dat files. Please refer
U1_nodes, U2_nodes etc. files for reference.

Elements.dat 		- Supporting files for running MATLAB code
Node_COORDS.dat 	- Supporting files for running MATLAB code
get_nod_adjele.m	- Supporting files for running MATLAB code
U1_nodes, U2_nodes etc. - Output files from MATLAB code
----------------------------------------------------------------------------
# Analysis using Abaqus NSFEM UEL
NSFEM_UEL_linear		- UEL for linear elastic analysis
NSFEM_UEL_NLGEOM_PLN_STRN	- UEL for 2D plane strain analysis
NSFEM_UEL_NLGEOM_PLN_STRSS 	- UEL for 2D plane stress analysis
mesh.dat			- Data file for unning Abaqus NSFEM UEL
---------------------------------------------------------------------------