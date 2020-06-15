%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main script that calls the diffusion analysis function for the 5mm
% control analysis
%
% Assumes s1_mrv2nii2fs_5mm.m and tksurfer scripts have already been run
%
% 12/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

control = 1;

%% Step 2 in the pipeline
s1_fs2wmDTI(control)

%% Step 3
% Note: this step takes forever
s1_createFDWM(control)

%% Step 4
s1_createTckmap(control)

%% Step 5
s1_transformTracks2Surface(control)

%% Step 6
s1_extractDensity(control)

%% Plotting and additional analyses done separately