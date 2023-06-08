# ChiaEtAl2023
Emergence of cortical network motifs for short-term memory during learning

### Hardware requirements:
As size of dataset may range from 500 MB - 24 GB, package requires standard computer with enough RAM to support the in-memory operations.
Check hardware requirements for MATLAB.

### OS requirements:
This package has been tested on the following systems:
	⁃	mac0S: Mojave 10.14.4
Check OS requirements for MATLAB

### Software requirements:
This package has been tested and on
	⁃	MATLAB R2017a
	⁃	MATLAB R2018b

### Installation guide:
Download and install MATLAB as instructed on the official website.
Installation time may vary.
Downloading of additional dependencies is not required as they have been embedded in the functions.

### Demo:
Open Figure_1.m
Set the path by adding with subfolders where the folder ‘code’ was saved.
```
	base = uigetdir ; % Set path folder ‘data’ (downloaded from Zenodo)
	data = get_data(parameters); 
```
In instances where (long runtime) is indicated, a pre-processed dataset will be available. You can choose to load the pre-processed dataset.
```
	load (‘non_matched_dataset.mat’); 
```
Once loaded, run the function to obtain ROC data for ‘stimulus’.
```
	ROC_population_data.stimulus = get_ROC_choice_stim_mode(parameters, data,’stimulus’);
```
Plot the figure.
```
	plot_auc_choice_stim_population(parameters, ROC_population_data.stimulus, parameters.sample_epoch)
```

## Correspondence
Hiroshi Makino, Lee Kong Chian School of Medicine, Nanyang Technological University
Email: hmakino@ntu.edu.sg
