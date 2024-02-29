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

## Code
The following analyses were carried out using open-source packages:

1. Surrogate networks that preserve the nodal degree and approximate edge length distribution of the empirical network were generated using the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/) [1].
2. The code for spatial autocorrelation-preserving surrogate brain maps can be implemented through the brainSMASH toolbox (https://github.com/murraylab/brainsmash) [2].
3. The code for gradient analysis and the Moran spectral randomization can be performed via BrainSpace (http://github.com/MICA-MNI/BrainSpace) [3].
4. The marmoset brain surfaces are presented using "plot_flatmap", available at https://github.com/netneurolab/marmoset_connectome [4].
5. The brain surfaces were visualized using Connectome Workbench (v1.5.0, https://www.humanconnectome.org/software/ connectome-workbench) [5]. 

## References
1. Rubinov, M., & Sporns, O. Complex network measures of brain connectivity: uses and interpretations. NeuroImage 52, 1059-1069 (2010).
2. Burt, J. B., Helmer, M., Shinn, M., Anticevic, A. & Murray, J. D. Generative modeling of brain maps with spatial autocorrelation. NeuroImage 220, 117038 (2020).
3. Vos de Wael, R. et al. Brainspace: a toolbox for the analysis of macroscale gradients in neuroimaging and connectomics datasets. Commun. Biol. 3, 103 (2020).
4. Liu, Z.-Q., Zheng, Y.-Q. & Misic, B. Network topology of the marmoset connectome. Netw. Neurosci. 4, 1181–1196 (2020).
5. Marcus, D. S. et al. Informatics and data mining tools and strategies for the human connectome project. Front. Neuroinformatics 5, 4 (2011).


