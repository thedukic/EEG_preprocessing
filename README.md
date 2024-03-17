# EEG_preprocessing
EEG preprocessing pipeline, ALS Centre UMC Utrecht.

# Get started
1. MATLAB version R2022b or later, on eariler versions the pipeline might fail

2. By default, CUDAICA is used as it is much faster than other implementations of ICA. This is possible (and useful) only if an (good, usually nonintegrated) NVIDIA graphics card is available. These are the steps:
	- Install the [CUDA toolkit](https://developer.nvidia.com/cuda-downloads/)
	- If this is not enough, then, download the [NuGet CLI](https://docs.microsoft.com/en-us/nuget/consume-packages/install-use-packages-nuget-cli/)
	- Save the file to an empty folder, e.g. C:\NuGet
	- Open the commmand prompt (cmd)
	- Using the "cd" command, native to the NuGet folder
	```cmd
	cd C:\NuGet
	```
	- Run in the cmd
	```cmd
	nuget install intelmkl.devel.win-x64
	```
	- In the NuGet folder, navigate to ...\intelmkl.redist.win-x64.2022.0.3.171\runtimes\win-x64\native
	- Copy these three files to EEGLAB\lugins\CudaICA folder:
		- mkl_core.2.dll
		- mkl_def.2.dll
		- mkl_intel_thread.2.dll
	
	If this does not or cannot work, in the preproc_parameters script you can always set the default version of ICA: cfg.ica.type = 'RUNICA'

3. EEGLAB pop_clean_rawdata edits so that the function is done on the GPU and with less display output/print
	- clean_asr:   change to: if ~exist('usegpu','var') || isempty(usegpu), usegpu = true; end
	- asr_process: change to: if nargin < 9 || isempty(usegpu), usegpu = true; end 
	- asr_process: change to: if usegpu % && length(range) > 1000
	- asr_process: change to: % if splits > 1, fprintf('.'); end

4. Set double precision in EEGLAB:
	- Open EEGLAB and go to File->Preferences
	- Select "show advanced options", click on OK and reopen EEGLAB
	- Unselect "use single precision number" and click on OK
	
5. Set folder paths in the preproc_folders script:
	- Root path of the pipeline folder: myfolders.mycodes
	- Root path of the raw data, which should be organised as on the bulk storage: myfolders.rootrawdata
	- Root path to the folder where preprocessed data will be saved: myfolders.rootpreproc

# License
Distributed under the GNU General Public License v3.0 License. See LICENSE.txt for more information.

# Acknowledgments
The pipeline is based on an exisiting preprocessing pipeline [RELAX](https://github.com/NeilwBailey/RELAX/).
