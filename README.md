# Cellular-functional decoupling in marmosets

This repository contains code created to conduct the main analyses described in the paper "**Decomposing Cortical Activity through Neuronal Tracing Connectome-eigenmodes in Marmosets**".

Uses an eigenmode-based analysis to study how neuronal tracing connectome (cellular connectome, CC) eigenmodes constrain cortical activity observed in BOLD-fMRI data. 

## System Requirements

Download the repository, and you're good to go. Read the comments and documentation within each code for usage guidance.

*MATLAB* (version R2018b: The MathWorks, Inc.) was used to write the code. 
The toolbox gspbox (please download gspbox from here: [gspbox](https://github.com/epfl-lts2/gspbox)) and the toolbox BrainSpace (please download BrainSpace from here: [BrainSpace](https://github.com/MICA-MNI/BrainSpace/releases)) and set its under the known Matlab paths with Matlab --> Set path --> Add folder with subfolders).

Some of the MATLAB-based scripts depend on packages developed by others. Copies of these packages have been stored in the `functions_matlab/` folder to ensure version compatibility. 

## Instructions

 `CFD_FullPipeline_demo.m`: MATLAB script to demonstrate how to use connectome eigenmodes to analyze fMRI data and compute cellular-functional decoupling (CFD) index.

We only provide  the demo dataset (3 subjects) to perform the calculations described in the paper. Original awake marmosets MRI datasets from the [Marmoset Brain Mapping](www.marmosetbrainmapping.org/data.html). Original human rs-fMRI data from HCP the [Human Connectome Project](https://db.humanconnectome.org/). Please consult the link for detailed information about access, licensing, and terms and conditions of usage.

The code for spatial autocorrelation correction can be implemented through the brainSMASH software (https://github.com/murraylab/brainsmash).

The code for gradient analysis and the Moran spectral randomization can be performed via BrainSpace (http://github.com/MICA-MNI/BrainSpace).

The marmoset brain surfaces are presented using "plot_flatmap", https://github.com/netneurolab/marmoset_connectome.





