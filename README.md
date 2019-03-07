# FACS_NGS

This project focused on using FACS-NGS method to profile the response of thousands of biosensor mutants. Briefly, Based on library consisted of thousands of synthetic sensor mutants, FACS-seq can be adopted to profile the response curve (reporter expressions vs. a series of ligand concentrations) for each individual mutant in living cells.

This script collection is user-friendly for experimental microbiologists with no or limited programming skills. Generally, the user only need to download the script, edit a configure file to set several parameters needed for sgRNA design, and type in one command line in a Linux environment to initiate the process. The output includes the statistics at sgRNA, gene and operon level of each studied phenotype. Meanwhile, diverse visualization files are also presented.

It consisted of two sub packages. 

One is the simulation package, used simulation method to guide the experiment design;

Another is the data processing package, used to convert the raw NGS .fq data from FACS-NGS experiments, based on the customized configure file to sensor response profile.

About the usage of each sub package, please check the readme.txt under relevant directory.
