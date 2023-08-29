### Description:
Demo Code required to train generalised linear model (GLM) for a given session and imaging field of view.

### Hardware requirements:
As size of dataset may range from 500mb - 24gb, package requires standard computer with enough RAM to support the in-memory operations.
Check hardware requirements for MATLAB

#### OS requirements:
This package has been tested on the following systems:
- mac0S: Mojave 10.14.4
Check OS requirements for MATLAB

#### Software requirements:
This package has been tested and on
- MATLAB R2017a
- MATLAB R2018b

### Installation guide:
Download and install MATLAB as instructed on the official website.
Installation time may vary.
Downloading of additional dependencies are not required.
For GLM cross-validation, cvpartitionImpl.m in /Applications/MATLAB_R2018b.app/toolbox/stats/stats/+internal/+stats line 113 needs to be changed to "properties(GetAccess = public, SetAccess = public)"

#### Demo:
Input:

1. Neural activity from Suite2p
2. Processed behaviour recording from Deeplabcut
3. Matched cell indices from RoiMatchPub
4. Experiment data from Wavesurfer
5. Behaviour data from Bpod

4 Steps procedures:

1. get_GLM_parameters.m  
Process raw data to run lassoGLM
Output: Variables required for running GLM. e.g. Spike activity, task predictors

2. run_GLM  
Perform GLM using lassoGLM from MATLAB.  
Output: Results from GLM. E.g. GLM coefficients for task variables

3. get_GLM_significance  
Determine statistical significance of a given task predictors by marginalising predictor  
Output: Pseudo explained variance of task predictors after removal

4. concat_GLM_results.m  
Output: Concatenated variables across sessions and animals
