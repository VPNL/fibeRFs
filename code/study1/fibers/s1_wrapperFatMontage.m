% This script generates montages of all the pairwise fibers
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all

% get our list of subjects from the Set function:
s1_setAllSessions

hems = {'rh' 'lh'};

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

retDir = [exptDir 'data/study1/toon'];
dtiDir = [exptDir 'data/study1/diffusion'];
fsDir = fullfile(RAID, '3Danat/FreesurferSegmentations'); 


%% Set up ROIs
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

ROIPre = 'fibeRFs_f_'; 

runName = {'96dir_run1'};

% pick radius for ROIs
rad = 1; % in order to extend the fibers within the voxel

for h = 1:length(hems)
    
    ROIs={};
    for r = 1:length(faceROIs) %face ROIs
        ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{r} '_projed_gmwmi_r1.00_fibeRFs_f_', hems{h}, '_EVC_projed_gmwmi']});
    end
    ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi_r1.00_fibeRFs_f_', hems{h}, '_EVC_projed_gmwmi']}); %add CoS places

    for n = 1:length(ROIs)
        ROIName = strsplit(ROIs{n});
        imgName = [ROIName{1}, '_r1.00_WholeBrainFG_19.tiff'];
        fatMontage(dtiDir, dMRI_sessions, runName, imgName, hems{h});
    end
end