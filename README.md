# EEG_preprocessing
EEG preprocessing pipeline v1, ALS Centre UMC Utrecht.

# Get started
1. Set folder paths in the preproc_folders.m script:
	- Root path of the pipeline folder: myPaths.mycodes
	- Root path of the raw data, which should be organised as on the bulk storage: myPaths.rootrawdata
	- Root path to the folder where preprocessed data will be saved: myPaths.rootpreproc

2. EEGLAB preferences are automatically set (see preproc_main.m)
	```MATLAB
	pop_editoptions('option_parallel',1,'option_single',0);
	```
		- option_parallel = 1; % this might be problematic if you have very large files
		- option_single   = 0; % do not change this
		
3. MATLAB version R2023b or later, on eariler versions the pipeline could fail.

4. By default, CUDAICA is used as it is much faster than other implementations of ICA. This is possible only if an (good, usually nonintegrated) NVIDIA graphics card is available. If this is not possible, the code will automatically switch to RUNICA, which is slower.


# License
Distributed under the GNU General Public License v3.0 License. See LICENSE.txt for more information.

# Acknowledgments
The pipeline is based on an exisiting preprocessing pipeline [RELAX](https://github.com/NeilwBailey/RELAX/).
