% This script is similar to s1_createFDWM but combines the FDWM for the
% EVC ROIs and the functional ROIs to determine the pairwise fibers
%
% It will 1) generate functionally defined fasciculus based on the
% connectome and fROIs and 2) render these fibers as an image using the
% subject anatomy
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

%% Loop through hemis 
for h = 1:length(hems)
    
    %% Get ROIs
    maps={};
    for r = 1:length(faceROIs) %face ROIs
        maps = horzcat(maps,{[ROIPre hems{h} '_' faceROIs{r} '_projed_gmwmi.mat']});
    end
    ROIs = horzcat(maps,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi.mat']}); %add CoS places
    
    inputROI = {[ROIPre hems{h} '_EVC_projed_gmwmi']}; %EVC ROI
    
    
    for s = 1:length(dMRI_sessions)
        
        if strcmp(dMRI_sessions{s}, 'KM25') || strcmp(dMRI_sessions{s}, 'MSH28') || ...
                strcmp(dMRI_sessions{s}, 'EM') || strcmp(dMRI_sessions{s}, 'GB23') ...
                || strcmp(dMRI_sessions{s}, 'DF') || strcmp(dMRI_sessions{s}, 'KGS') ...
                || strcmp(dMRI_sessions{s}, 'MG') || strcmp(dMRI_sessions{s}, 'MJH25') ...
                || strcmp(dMRI_sessions{s}, 'MN') || strcmp(dMRI_sessions{s}, 'SP') ...
                || strcmp(dMRI_sessions{s}, 'MBA24')
            t1_name = ['t1.nii.gz'];
        else
            t1_name = ['T1_QMR_1mm.nii.gz'];
        end

        %% 2) Define FDFs and get fiber count for pairwise
        for r=1:length(runName)
            fgName=[inputROI{1}  '_r' num2str(rad) '.00_WholeBrainFG.mat'];
            fgDir=fullfile(dtiDir, dMRI_sessions{s},runName{r},'dti96trilin','fibers','afq');
            for n=1:length(ROIs)
                ROI = ROIs(n);
                %deal with IOG's overlapping fiber problem
                if ~isempty(strfind(ROI{1},'IOG')) || ~isempty(strfind(ROI{1},'pFus')) || ~isempty(strfind(ROI{1},'CoS'))
                    fatFiberIntersectRoi_EVCpairwise(dtiDir, fgDir, dMRI_sessions{s}, runName{r}, fgName, ROI, 1, rad)
                else
                    fatFiberIntersectRoi(dtiDir, fgDir, dMRI_sessions{s}, runName{r}, fgName, ROI, 1, rad)
                end
            end
        end

        %% 3) Display the pairwise fibers
        for r=1:length(runName)
            for n=1:length(ROIs)
                ROIName=strsplit(ROIs{n},'.');
                ROIfg= [ROIName{1} '_r' num2str(rad) '.00_' inputROI{1} '_r' num2str(rad) '.00_WholeBrainFG.mat'];
                foi=1; %as we want all tracts
                fatRenderFibers(dtiDir, dMRI_sessions{s}, runName{r}, ROIfg, foi,t1_name, hems{h});
            end
        end

        close all

    end
end