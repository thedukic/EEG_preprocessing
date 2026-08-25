# Automated EEG Preprocessing & Postprocessing Pipeline

An automated, two-stage EEG analysis pipeline developed at the ALS Centre, University Medical Centre Utrecht. Engineered for high-density 128-channel BioSemi ActiveTwo recordings, this pipeline standardises data curation, multi-modal artefact rejection, individualised spectral metrics, and objective quality assurance.

---

## Pipeline Workflow

```
Raw BioSemi Data (.bdf)
   │
   ▼
[ Part 1: Continuous Cleaning (run_subject_1) ]
   ├── Integrity checks: CMS dropouts, DC offsets, and flat electrodes
   ├── Downsampling (256 Hz) & High-pass filtering
   ├── Dual-stage spectral line-noise suppression
   ├── External lead derivation (Bipolar VEOG, HEOG, ECG)
   ├── Spatial artefact filtering (STAR, GEDAI)
   ├── Bad channel detection & Spherical spline interpolation
   ├── ICA decomposition (CUDAICA / RUNICA) & IC classification and removal
   └── Re-referencing (Average) & Interim export
   │
   ▼
[ Part 2: Postprocessing & Quality Assurance (run_subject_2) ]
   ├── Continuous low-pass filtering
   ├── Individual Alpha Frequency (IAF) estimation
   ├── Data epoching
   ├── Baseline correction & Epoch rejection
   ├── Residual artefact screening (Ocular & Muscular slopes)
   ├── Automagic quality scoring (RBC, OHA, THV, CHV)
   └── BIDS metadata integration & Final export
```

---

## Core Methodology

* **Dual-Stage Spectral Cleaning:** Two-pass line noise attenuation targeting 50 Hz fundamental and harmonic interference without compromising neighbouring neural power.
* **Spatial Filtering:** Early spatial artefact suppression using STAR and GEDAI routines prior to ICA decomposition.
* **Flexible Cardiac Derivation:** True bipolar derivation from dedicated leads (`ECGL`/`ECGR`), with an automatic 2-lead PCA fallback for earlobe references (`LEL`/`REL`).
* **Subspace Rank Verification:** PCA rank verification and quantitative residual checking (1/f spectral slope and excess kurtosis) to prevent loss of neural variance during dimensionality reduction.
* **Accelerated Decomposition:** Native CUDAICA support for GPU-accelerated ICA, with automated CPU fallback to RUNICA.
* **Standardised Quality Control:** Automated computation of objective Automagic metrics (RBC, OHA, THV, CHV) alongside high-frequency EMG leftover quantification.

---

## Data & Montage Specifications

* **Scalp Montage:** Optimised for 128-channel BioSemi ActiveTwo systems using radial ABCD pinout labelling.
* **Auxiliary Channels:** 
  * Ocular: 2x VEOG (superior/inferior) and 2x HEOG (left/right).
  * Cardiac: Dedicated leads (`ECGL`/`ECGR`) or earlobe references (`LEL`/`REL`).
  * Muscular: Optional surface EMG channels for motor tasks.
* **Supported Formats:** BioSemi raw (`.bdf`) and EEGLAB (`.set`).

---

## Getting Started

### 1. Prerequisites
* MATLAB R2023b, R2024b, or R2025b.
* **Required Toolboxes:** Signal Processing, Statistics and Machine Learning, Parallel Computing.
* *(Optional)* NVIDIA GPU with CUDA support for CUDAICA.

### 2. External Dependencies Setup
Create an `external` folder in the root directory of this repository (or at your configured path) and extract the required third-party toolboxes into it:

```
EEG_preprocessing/
├── external/
│   ├── brewermap-3.2.8/
│   ├── eeglab2025.1.0/
│   ├── gedai_05082026/
│   ├── noisetools_29-Apr-2023/
│   ├── restingiaf_20-Jan-2025/
│   └── zaplineplus_14-Apr-2023/
├── files/
├── preproc_folders.m
├── preproc_main.m
├── run_subject_1.m
└── run_subject_2.m
```

The pipeline will automatically index, add, and verify these subfolders during initialisation.

### 3. Reset MATLAB Search Path
To prevent namespace collisions with existing external toolboxes on your system, reset your MATLAB search path to defaults before launching:

```matlab
pathtool; % Click 'Default', then 'Save' and close.
```

### 4. Configure Project Paths
Open `preproc_folders.m` and define your local directory structure:

```matlab
myPaths.mycodes     = 'C:/matlab/codes/EEG_preprocessing/';
myPaths.rootrawdata = '/data/EEG/raw/';
myPaths.rootpreproc = '/data/EEG/preprocessed/';
```

### 5. Running the Pipeline

**Single Subject Execution:**
Run Part 1 and Part 2 sequentially for an individual participant ID:

```matlab
% 1. Load paths and initialise environment
myPaths = preproc_folders;

% 2. Run continuous cleaning and ICA (Part 1)
run_subject_1(myPaths, 'SUBJ001');

% 3. Run epoching, spectral metrics, and QA (Part 2)
run_subject_2(myPaths, 'SUBJ001');
```

**Batch Processing:**
To batch process an entire cohort across groups, visits, and task blocks with automated error tracking:

```matlab
% Configure cohort variables in preproc_folders.m and launch batch runner
preproc_main;
```

---

## Hardware & Execution Notes

* **Parallel Processing:** Multi-core execution is enabled by default via `pop_editoptions('option_parallel', 1)`.
* **CUDA Acceleration:** The pipeline checks for an NVIDIA GPU automatically. If available, matrix decomposition runs via `cudaica`; otherwise, it falls back to standard `runica`.

---

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE.txt](LICENSE.txt) file for details.

---

## Acknowledgments

This pipeline builds upon methods, toolboxes, and algorithms developed across the neuroimaging and electrophysiology community:

* [EEGLAB](https://github.com/sccn/eeglab/) - Core data structures and signal processing routines.
* [Zapline-plus](https://github.com/MariusKlug/zapline-plus/) - Adaptive spectral line-noise removal.
* [Noise Tools](http://audition.ens.fr/adc/NoiseTools/) - Advanced spatial filtering and detrending.
* [GEDAI](https://github.com/neurotuning/GEDAI-master/) - Generalized eigendecomposition spatial filtering and artefact attenuation.
* [PREP Pipeline](https://vislab.github.io/EEG-Clean-Tools/) - Robust referencing and bad channel heuristics.
* [restingIAF](https://github.com/corcorana/restingIAF/) - Individual alpha peak estimation.
* [RELAX](https://github.com/NeilwBailey/RELAX/) - Automated artefact cleaning protocols.
* [Automagic](https://github.com/methlabUZH/automagic) - Standardised objective EEG quality metrics.
* [MWF](https://github.com/exporl/mwf-artifact-removal/) - Spatiotemporal multichannel Wiener filtering.
* [BrewerMap](https://github.com/DrosteEffect/BrewerMap/) - Perceptually uniform colour palettes.