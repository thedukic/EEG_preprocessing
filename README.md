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
	myPaths.mycodes     = '/home/user/EEG_preprocessing/';
	myPaths.rootrawdata = '/data/EEG/ALS/raw_data/';
	myPaths.rootpreproc = '/data/EEG/ALS/preprocessed_data/';
	```
	
3. Ensure MATLAB compatibility. The pipeline is tested on MATLAB R2023b and R2024b.

# Additional information: 
- The pipeline is optimised for 128-channel EEG data from [BioSemi systems](https://www.biosemi.com/).
- EEGLAB preferences are automatically set (see preproc_folders.m):
	```MATLAB
	pop_editoptions('option_parallel',1,'option_single',0);
	```
	- option_parallel = 1; % This might be problematic if you have very large files. The script turns off this option for motor data automatically.
	- option_single   = 0; % Do not change this, it is very important for some methods, such as ICA, to use double precision.
- By default, the pipeline uses CUDAICA for faster ICA processing. This is only possible with a dedicated NVIDIA graphics card; if one is not detected, the code will automatically switch to the slower RUNICA.

# License
This project is licensed under the GNU General Public License v3.0. For details, see the LICENSE.txt file.

# Acknowledgments
This pipeline draws inspiration from and builds upon existing knowledge from the following repositories:
- [EEGLAB](https://github.com/sccn/eeglab/)
- [Zapline-plus](https://github.com/MariusKlug/zapline-plus/)
- [Noise Tools](http://audition.ens.fr/adc/NoiseTools/)
- [PREP pipeline](https://vislab.github.io/EEG-Clean-Tools/)
- [MWF](https://github.com/exporl/mwf-artifact-removal/)
- [restingIAF](https://github.com/corcorana/restingIAF/)
- [BrewerMap](https://github.com/DrosteEffect/BrewerMap/)
- [RELAX](https://github.com/NeilwBailey/RELAX/)

