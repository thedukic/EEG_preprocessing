# EEG Preprocessing Pipeline
This repository contains an automated pipeline for EEG preprocessing, developed at the ALS Centre UMC Utrecht. The pipeline is designed to be user-friendly, reliable, and grounded in best practices for removing common EEG artefacts.

# Get started
1. Set your MATLAB paths to default. This ensures the pipeline starts from a clean slate.
    ```MATLAB
	pathtool; % Then press 'Default' and confirm.
	```

2. Update the folder paths at the top of the preproc_folders.m script
    ```MATLAB
	% Example paths: update these for your system
	myPaths.mycodes = '/home/user/my-eeg-pipeline/';
	myPaths.rootrawdata = '/media/bulk/ALS/raw_data/';
	myPaths.rootpreproc = '/media/bulk/ALS/preprocessed_data/';
	```
	
3. Ensure MATLAB compatibility. The pipeline requires MATLAB R2023b or later.

# Additional information: 
- EEGLAB preferences are automatically set (see preproc_folders.m):
	```MATLAB
	pop_editoptions('option_parallel',1,'option_single',0);
	```
	- option_parallel = 1; % This might be problematic if you have very large files. The script turns off this option for motor data automatically.
	- option_single   = 0; % Do not change this, it is very important for some methods, such as ICA, to uyse double precision.
- By default, CUDAICA is used as it is much faster than other implementations of ICA. This is possible only if an NVIDIA graphics card is available. If this is not possible, the code will automatically switch to RUNICA, which is slower.

# License
This project is licensed under the GNU General Public License v3.0. For details, see the LICENSE.txt file.

# Acknowledgments
This pipeline is based on and relies on existing code from the following repositories:
- [EEGLAB](https://github.com/sccn/eeglab/)
- [RELAX](https://github.com/NeilwBailey/RELAX/)
- [Zapline-plus](https://github.com/MariusKlug/zapline-plus/)
- [MWF](https://github.com/exporl/mwf-artifact-removal/)
- [restingIAF](https://github.com/corcorana/restingIAF/)
- [BrewerMap](https://github.com/DrosteEffect/BrewerMap/)
