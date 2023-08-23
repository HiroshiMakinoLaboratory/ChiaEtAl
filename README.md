# ChiaEtAl2023NatureCommunications
Emergence of cortical network motifs for short-term memory during learning

### Dataset:
Dataset can be downloaded at: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8031278.svg)](https://doi.org/10.5281/zenodo.8031278)


### Hardware requirements:
As size of dataset may range from 500 MB - 24 GB, package requires standard computer with enough RAM to support the in-memory operations.
Check hardware requirements for MATLAB.

### OS requirements:
This package has been tested on the following systems:
- mac0S: Mojave 10.14.4
Check OS requirements for MATLAB

### Software requirements:
This package has been tested and on
- MATLAB R2017a
- MATLAB R2018b

### Installation guide:
Download and install MATLAB as instructed on the official website.
Installation time may vary.
Downloading of additional dependencies are not required as they have been embedded in the functions.

### Demo:
Open Figure_5.m\
Set the path by adding with subfolders.
```
	parameters.root_folder = uigetdir; % CD to folder where files from Zenodo were downloaded
	cd(parameters.root_folder)
```
In instances where (long runtime) are indicated, a pre-processed dataset will be available. You can choose to load the pre-processed dataset.
```	
	% Load raw dataset containing raw activity
	load 'retained_eliminated_data_raw.mat' 
	
	% Get reconstructed retained and eliminated activity (long runtime)
	retained_elim_data = get_retained_eliminated_activity(parameters,data); 
	
	% Or load pre-processed dataset
	load 'retained_elim_data.mat'
```
Once loaded, run the function to analyse and plot
```
	% Plot average number of retained and lost functional coupling
	[prop_retained,prop_lost] = plot_num_retained_lost_cells(retained_elim_data);
```

## Correspondence
Hiroshi Makino, Lee Kong Chian School of Medicine, Nanyang Technological University
Email: hmakino@ntu.edu.sg
