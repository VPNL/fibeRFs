%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main script that calls the diffusion analysis function for the 10mm
% mSTS control analysis
%
% Assumes s1_mrv2nii2fs and tksurfer scripts have already been run
%
% 10/2020 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mSTS_control = 2;

%% Step 2 in the pipeline
s1_fs2wmDTI(mSTS_control)

%% Step 3
% Note: this step takes forever
s1_createFDWM(mSTS_control)

%% Step 4
s1_createTckmap(mSTS_control)

%% Step 5
s1_transformTracks2Surface(mSTS_control)

%% Step 6
s1_extractDensity(mSTS_control)

%% Plotting and additional analyses done separately