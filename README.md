# COMP 7850 - Project
Agglomerative hierarchical clustering project for COMP 7850: https://github.com/vjayaman/ahc\_project

To run, change NUM\_VECS value in *functions.h* and edit lines 6 and/or 8 in *myjob* file.
(Note: using *myjob* and *makefile* files provided by professor for assignments). 

Runs *hybrid\_clustering.R*, the file that actually runs the agglomerative (parallel and sequential) clustering implementations; calls on functions named and defined in *functions.h, functions.c*.

*for\_demo.R* is a very basic R Shiny app that plots and displays tables for the already-generated output files.

*figs\_and\_tbls/* contain figures and tables used in the project report (as well as a few that were considered and not used). 

*result_processing/* contains the base files used for processing the outputs of the MPI runs - functions of each are saved to rfuncs.R, for easier use during app runtime. Includes a function that runs the *fastcluster* R package on the same matrices as the sequential implementation run in *hybrid\_clustering.R*. (*rfuncs.R* is better commented than the other files).
