%% Preliminaries
% Use the UMS2K package for spike detection and sorting
here =fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here,'UMS2K')),'-begin');
% Use the fMRIb package for spike detection and sorting
addpath(fullfile(here,"fMRIb"))

%% UMS2k is a bit generous with showing progress bars that take mouse focus.
% Hide those here
global  hideUMS2kProgress
hideUMS2kProgress = true;
