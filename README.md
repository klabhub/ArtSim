# ArtSim
A MATLAB toolbox for simulating and eliminating artifacts caused by electrical stimulation in EEG, LFP, or single-unit recordings.

The companion publication (under review) provides background, interpretation, and illustrative example simulations:

***Artifact Correction for Transcranial Current Stimulation. Bart Krekelberg, 2025.***


_Bart Krekelberg, Rutgers University - Newark, 2025._
[https://vision.rutgers.edu](vision.rutgers.edu)

## Installation
This code was developed and tested in MATLAB R2025a, and is backward compatible only to R2021a  (mainly due to the use of arguments blocks in function definitons). 
To install, clone the repository, _including the submodules_, to your machine :
```
git clone --recurse-submodules https://github.com/klabhub/artSim
```
This will include a fork of the spike detection and sorting pipeline implemnted in [UMS2K](https://github.com/danamics/UMS2K) and the FASTR artifact removal algorithm implemented in [fMRIb](https://github.com/sccn/fMRIb.git). 

## Quick Start
After cloning this repository, open Matlab, `cd` into the `/code` folder and run `startup.m` to add the subfolders to the MATLAB search path.
Then open one of the LiveScripts that produce the main analyses of the companion publication:

- `eegArtifacts.mlx` - Analysis of common stimulation artifacts in EEG or LFP recordings.
- `spikeArtifacts.mlx`  - Analysis of commom stimulation artifacts in single-unit spike recordings.


## Funding
This research was supported by the [National Institute of Neurological Disorders and Stroke](https://www.ninds.nih.gov/) and the [National Institute on Drug Abuse](https://nida.nih.gov/) of the National Institutes of Health under Award Number R01NS120289. The content is solely the responsibility of the author and does not necessarily represent the official views of the National Institutes of Health. The funding institution did not play any role in the study design, data collection and analysis, decision to publish, or preparation of the manuscript. 


## Notes

### fmrib_fastr  has a few peculiarities

1. Segmentation: the segments partially overlap  (this is done to allow realignment).

2. When using slice-triggers, the segments used to determine the mean artifact step in
steps of 2* length of a segment (i.e. skipping one ), when using volume-triggers, the step size is a single segment. (As this changes the effective frequency of the artifact sampling, it afffects which harmonics of the tACS frequency will be removed; see the step parameter in movingAverage.mlx).

3. When using volume-triggers (but not with slice-triggers), the code regresses out the mean (per segment).
