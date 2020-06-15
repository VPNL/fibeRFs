% Step 1: 
% Preparing the mrVista ROIs for use in FDWM
% Transforms them from the mrVista space to the freesurfer space
%
% Note that after this step you will need to open the ROIs in tksurfer 
% and save with _surf.label before you can run Step 1 (see
% github.com/VPNLtools/tksurfer_scripting for how to do this)
%
% This script is kept separate between main analysis and control because of
% the tksurfer step in between
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

fROIPre = 'fiberRFsclean_f_'; 

runName = {'96dir_run1'};

%% Loop through hemis 
for h = 1:length(hems)
    
    %% Get ROIs
    maps={};
    for r = 1:length(faceROIs) %face ROIs
        maps = horzcat(maps,{[fROIPre hems{h} '_' faceROIs{r} '_5mm']});
    end
    maps = horzcat(maps,{[fROIPre hems{h} '_' placeROIs{1} '_5mm']}); %add CoS places


    for s = 1:length(dMRI_sessions)
        
        subjectDir = fullfile(retDir, dMRI_sessions{s}, runName, 't1');

        %% First, convert the mrVista ROI into a nifti file
        session = fullfile(retDir, dMRI_sessions{s});
        load(fullfile(session, 'mrSESSION.mat'));
        [t1folder, name, ext] = fileparts(vANATOMYPATH);
        t1name = [name, ext]; clear name ext vANATOMYPATH mrSESSION dataTYPES
        dataPath = fullfile(session, '3DAnatomy', 'niftiROIs');

        mrVistaROI2niftiROI(session, maps, t1name, t1folder, dataPath) %convert them all to niftis

        %% Now we will convert the nifti ROI into a freesurfer label
        subject = fs_sessions{s};
        cd(dataPath)

        for r = 1:length(maps) %loop through the ROIs
            if exist([maps{r}, '.nii.gz'], 'file')
                ni = readFileNifti(maps{r});
            else
                continue
            end

            roiVal = 1;
            labelPath = fullfile('/biac2/kgs/3Danat/FreesurferSegmentations/', subject, '/label/');
            
            %transform to fs label
            nii2labelProj(subject, dataPath, maps(r), roiVal, labelPath, hems{h})
            
        end

        % after this you need to open stuff in tksurfer and save it with _surf.label before you can run the rest of the code
    end
end