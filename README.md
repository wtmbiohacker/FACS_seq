# FACS_seq

This project focused on using FACS-seq method to profile the response of thousands of biosensor mutants. Briefly, Based on library consisted of thousands of synthetic sensor mutants, FACS-seq can be adopted to profile the response curve (reporter expressions vs. a series of ligand concentrations) for each individual mutant in living cells. The paper describing this program is currently under review. Please cite this paper if this program is useful to your work.

For biologists interested in quantifying the functional consequences of many mutations (1,000 ~ 100,000) of a genetic regulatory element (transcription factor, riboswitch, DNA binding site, upstream ORF, etc), FACS-seq strategy can be used, together with synthetic library, to quantify the responses of each individual mutant in the library when a series of ligand concentrations are used. This script collection is user-friendly for experimental biologists with no or limited programming skills. Generally, the user only need to download the script, edit a configure file to set several parameters needed for running the code, and type in one command line in a Linux environment to initiate the process. The output includes the final result, namely, the response of each individual mutant in each tested condition (ligand concentration), as well as quality control visualization files such as replicate agreement analysis.

The software consisted of two sub packages. 

One is the simulation package, used simulation method to guide the experiment design;

Another is the data processing package, used to convert the raw NGS .fq data from FACS-NGS experiments, based on the customized configure file to sensor response profile.

About the usage of each sub package, please check the readme under relevant directory.

For any request or question, please contact wtm0217@gmail.com
