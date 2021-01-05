function s1_transformTracks2Surface(control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: 
% Transforms the trackmaps to the cortical surface
%
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    elseif control == 2
        ROIs = {[ROIPre hems{h} '_mSTS_faces_10mm_projed_gmwmi.mat']}; 
    else
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' faceROIs{r} '_projed_gmwmi.mat']});
        end
        ROIs = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_projed_gmwmi.mat']}); %add CoS places
    end

    % loop through sessions and transform maps to fsaverage surfaces using CBA
    for s = 1:length(dMRI_sessions)
        
        fs_id = fs_sessions{s};
        % path to subject data in 3Danat
        anat_dir = fullfile(dtiDir, dMRI_sessions{s}, '96dir_run1/t1');
        % path to subject data in FreesurferSegmentations
        fs_dir = fullfile(fsDir, fs_id);
        % paths to subject mri and surf directories
        mri_dir = fullfile(fs_dir, 'mri');
        surf_dir = fullfile(fs_dir, 'surf');

        track_dir = fullfile(anat_dir, 'tracks');

        for n = 1:length(ROIs)
            ROIName = strsplit(ROIs{n}, '.');

            map_name = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG_track'];
            trkName = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG.tck'];
            fg = fullfile(dtiDir, dMRI_sessions{s}, runName{1}, 'dti96trilin/fibers/afq', trkName)

            if exist(fg)
                map_path = fullfile(anat_dir, 'tracks', [map_name, '_resliced.nii.gz']);
                cd(fullfile(anat_dir, 'tracks'));
                unix(['mri_convert -ns 1 -odt float -rt interpolate -rl ', mri_dir, '/orig.mgz ', ...
                    track_dir, '/', map_name, '.nii.gz ', track_dir, '/', map_name, '_resliced.nii.gz --conform']);

                % generate freesurfer-compatible surface files for each hemisphere
                cd(surf_dir);

                proj_values = [-0.5:0.1:0.5];
                for p = 1:length(proj_values)
                    unix(['mri_vol2surf --mov ', map_path, ' ', ...
                            '--reg register.dat --hemi ' hems{h} ' --interp trilin --o ', ...
                            fullfile(surf_dir, strcat(map_name, '_', hems{h}, '_proj_', num2str(proj_values(p)), '.mgh')), ...
                            ' --projfrac ', num2str(proj_values(p))]); % left hemi
                end
                
                if exist(strcat(map_name, '_', hems{h}, '_proj_max.mgh'))
                    prev_max_map = strcat(map_name, '_', hems{h}, '_proj_max.mgh');
                    delete(prev_max_map)
                end
                
                unix(['mri_concat --i ', strcat(map_name, '_', hems{h}, '_proj_*'), ...
                        ' --o ', strcat(map_name, '_', hems{h}, '_proj_max.mgh'), ...
                        ' --max']);
            end

        end
    end
end

clear all
close all
