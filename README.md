# EEG_preprocessing
EEG preprocessing pipeline, ALS Centre UMC Utrecht.

# Get started
1. MATLAB version R2022b or later, on eariler versions the pipeline might fail
2. By default, CUDAICA is used as it is much faster than other implementations of ICA. This is possible (and useful) only if an (good, usually nonintegrated) NVIDIA graphics card is available. These are the steps:
	a. Install the [CUDA toolkit](https://developer.nvidia.com/cuda-downloads/)
	b. If this is not enough, then, download the [NuGet CLI](https://docs.microsoft.com/en-us/nuget/consume-packages/install-use-packages-nuget-cli/)
	c. Save the file to an empty folder
	d. Open the commmand prompt (cmd)
	e. Using the "cd" command, native to the NuGet folder
	f. Run in the cmd "nuget install intelmkl.devel.win-x64"
	g. In the NuGet folder, navigate to ...\intelmkl.redist.win-x64.2022.0.3.171\runtimes\win-x64\native
	g. Copy these three files to EEGLAB\lugins\CudaICA folder:
		- mkl_core.2.dll
		- mkl_def.2.dll
		- mkl_intel_thread.2.dll
	
	If this does not or cannot work, in the preproc_parameters script you can always select the default version of ICA: RUNICA 
	
3. Set double precision in EEGLAB:
	a. Open EEGLAB and go to File->Preferences
	b. Select "show advanced options", click on OK and reopen EEGLAB
	c. Unselect "use single precision number" and click on OK
	
4. Set folder paths in the preproc_folders script:
	a. Root path of the pipeline folder: myfolders.mycodes
	b. Root path of the raw data, which should be organised as on the bulk storage: myfolders.rootrawdata
	c. Root path to the folder where preprocessed data will be saved: myfolders.rootpreproc

# License
Distributed under the GNU General Public License v3.0 License. See LICENSE.txt for more information.

# Acknowledgments
The pipeline is based on an exisiting preprocessing pipeline [RELAX](https://github.com/NeilwBailey/RELAX/)
