# ArtSim

A MATLAB toolbox for simulating and eliminating artifacts caused by electrical stimulation in EEG, LFP, or single-unit recordings.

The companion publication (under review) provides background, interpretation, and illustrative example simulations:

***ArtSim: Open-source software for evaluating artifact correction in EEG, LFP, and spike recordings acquired with transcranial current stimulation. Bart Krekelberg, 2026.***

*Bart Krekelberg, Rutgers University - Newark, 2025.*
[https://vision.rutgers.edu](vision.rutgers.edu)

## Installation

This code was developed and tested in MATLAB R2025a, and is backward compatible only to R2021a  (mainly due to the use of arguments blocks in function definitions).
To install, clone the repository,*including the submodules*, to your machine:

```
git clone --recurse-submodules https://github.com/klabhub/artSim
```

This will include a fork of the spike detection and sorting pipeline implemented in [UMS2K](https://github.com/danamics/UMS2K) and the FASTR artifact removal algorithm implemented in [fMRIb](https://github.com/sccn/fMRIb.git).

## Quick Start

After cloning this repository, run `startup.m` to add the subfolders to the MATLAB search path.
Then open one of the LiveScripts in the `demo` folder that produce the main analyses of the companion publication:

- `eeg.m` - Analysis of common stimulation artifacts in EEG or LFP recordings.
- `spikes.x`  - Analysis of common stimulation artifacts in single-unit spike recordings.

The other Livescripts provide a deeper dive into some parameter dependencies:

- `eegExploration.m`
- `spikesExploration.m`

or they explore other issues associated with TCS artifacts:

- `movingAverage.m`   - Comb artifacts resulting from moving average templates
- `synchrony.m`       - A closer look at the statistical analysis of spike entrainment

The `restingState.m` file shows an example application of the toolbox to resting state EEG data.

## Funding

This research was supported by the [National Institute of Neurological Disorders and Stroke](https://www.ninds.nih.gov/) and the [National Institute on Drug Abuse](https://nida.nih.gov/) of the National Institutes of Health under Award Number R01NS120289. The content is solely the responsibility of the author and does not necessarily represent the official views of the National Institutes of Health. The funding institution did not play any role in the study design, data collection and analysis, decision to publish, or preparation of the manuscript.

## Notes

### fmrib_fastr  has a few peculiarities

1. Segmentation: the segments partially overlap  (this is done to allow realignment).

2. When using slice-triggers, the segments used to determine the mean artifact step in
steps of 2* length of a segment (i.e. skipping one ), when using volume-triggers, the step size is a single segment. (As this changes the effective frequency of the artifact sampling, it afffects which harmonics of the tACS frequency will be removed; see the step parameter in movingAverage.mlx).

3. When using volume-triggers (but not with slice-triggers), the code regresses out the mean (per segment).

## Developers

The `/test` folder has a unit test suite. When extending the toolbox, use `run('test/runAllTests.m')` to run them all and confirm that the code still generates the same outputs.
