function s1_createFDWM(control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3
%
% This script is the first step in processing the fiber tracts after
% generating the connectome. The connectome should have already been
% created by seeeding from the gray-white matter boundary and validated by
% ET/LiFE (currently this processing is done in flywheel). The ROIs of
% interest are also required (can still be mrVista ROIs).
%
% It will 1) generate functionally defined fasciculus based on the
% connectome and fROIs and 2) render these fibers as an image using the
% subject anatomy
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if notDefined('control')
    control = 0;
end

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

ROIPre = 'fibeRFsclean_f_'; 

runName = {'96dir_run1'};

% pick radius for ROIs
rad = 1; % in order to extend the fibers within the voxel

%% Loop through hemis 
for h = 1:length(hems)
    
    %% Get ROIs
    maps={};
    if control == 1
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{[ROIPre hems{h} '_' faceROIs{r} '_5mm_projed_gmwmi.mat']});
        end
        ROIs = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_5mm_projed_gmwmi.mat']}); %add CoS places
    elseif control == 2  %10mm control analysis for mSTS
        ROIs = {[ROIPre hems{h} '_mSTS_faces_10mm_projed_gmwmi.mat']};
    else
        maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_EVC_projed_gmwmi.mat']});
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' faceROIs{r} '_projed_gmwmi.mat']});
        end
        ROIs = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_projed_gmwmi.mat']}); %add CoS places
    end
    acpcROI = {['toon_f_' hems{h} '_EVC.mat']}; %an ROI I know is correctly aligned

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
        
        if control ~= 1 %main analysis only
            if strcmp(hems{h},'rh') %these steps only need to be done once (and I always start with the right hemi)

                % The flywheel output has a bunch of hardcoded/incorrect paths
                % We need to fix the T1 path before we can continue
                load(fullfile(dtiDir, dMRI_sessions{s}, runName{1},'dti96trilin', 'dt6.mat'));
                files.t1 = fullfile('t1', t1_name);
                files.aligned = fullfile('dti96trilin/mrtrix/dwi_aligned_trilin_noMEC_gmwmi.nii.gz');
                save(fullfile(dtiDir, dMRI_sessions{s}, runName{1},'dti96trilin', 'dt6.mat'), 'files', 'params', 'adcUnits');

                %% 1) Convert vista ROI to dti ROI to get transformation matrix
                for r=1:length(runName)
                    if strcmp(dMRI_sessions{s}, 'CR24')
                        fatVistaRoi2DtiRoi(dtiDir, dMRI_sessions{s}, runName{r}, acpcROI, t1_name,2) %need dims from gmwmi file
                    else
                        fatVistaRoi2DtiRoi(dtiDir, dMRI_sessions{s}, runName{r}, acpcROI, t1_name,1) %normal
                    end
                end

            end
        end

        % pick radius for ROIs
        rad = 1; % in order to extend the fibers within the voxel

        %% 2) Define FDFs and get fiber count
        for r = 1:length(runName)
            fgName = ['WholeBrainFG.mat'];
            fgDir = fullfile(dtiDir, dMRI_sessions{s}, runName{r}, 'dti96trilin', 'fibers');
            fatFiberIntersectRoi(dtiDir, fgDir, dMRI_sessions{s}, runName{r}, fgName, ROIs, 1, rad)
        end

        %% 3) Display the fibers
        for r = 1:length(runName)
            if control == 0
                for n = 2:length(ROIs) %don't plot EVC (way too many fibers)
                    ROIName = strsplit(ROIs{n}, '.');
                    ROIfg = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG.mat'];
                    foi = 1; %as we want all tracts
                    fatRenderFibers(dtiDir, dMRI_sessions{s}, runName{r}, ROIfg, foi, t1_name, hems{h});
                end
            else
                for n = 1:length(ROIs) %don't plot EVC (way too many fibers)
                    ROIName = strsplit(ROIs{n}, '.');
                    ROIfg = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG.mat'];
                    foi = 1; %as we want all tracts
                    fatRenderFibers(dtiDir, dMRI_sessions{s}, runName{r}, ROIfg, foi, t1_name, hems{h});
                end
            end
        end

        close all

    end
end

clear all
close all