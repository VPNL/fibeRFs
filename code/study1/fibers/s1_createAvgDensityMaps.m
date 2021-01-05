% This creates an average density map on the fs average surface
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
fsaDir = fullfile(fsDir, 'fsaverage-bkup', 'surf');
%% Set up ROIs
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

ROIPre = 'fibeRFs_f_'; 

%% loop through sessions and transform maps to fsaverage surfaces using CBA
for h = 1:length(hems)

    ROIs={};
    r=5;
    ROIs = horzcat(ROIs,{['fibeRFsclean_f_' hems{h} '_' faceROIs{r} '_10mm_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max']});
%     for r = 1:length(faceROIs) %face ROIs
%         ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{r} '_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max']});
%     end
%     ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max']}); %add CoS places

    for ss = 1:length(fs_sessions)
        
        fs_id = fs_sessions{ss}; 

        % path to subject data in FreesurferSegmentations
        subjDir = fullfile(fsDir, fs_id);
        % paths to subject mri and surf directories
        mriDir = fullfile(subjDir, 'mri'); surfDir = fullfile(subjDir, 'surf');

        for mm = 1:length(ROIs)
            map_name = ROIs{mm}; 

            cd(surfDir);

            % transform surface files to fsaverage
            map_stem = fullfile(fsaDir, map_name);
            unix(['mri_surf2surf --srcsubject ' fs_id ' --srcsurfval ' ...
                map_name '.mgh --trgsubject fsaverage-bkup --trgsurfval ' ...
                map_stem '_regFrom_' fs_id '.mgh --hemi ' hems{h}]);

        end
    end
end

%% average surface maps across all sessions
cd(fullfile(fsDir, 'fsaverage-bkup', 'surf'));
for h = 1:length(hems)
    
    ROIs={};
%     for r = 1:length(faceROIs) %face ROIs
%         ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{r} '_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max']});
%     end
%     ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max']}); %add CoS places
r=5;    
ROIs = horzcat(ROIs,{['fibeRFsclean_f_' hems{h} '_' faceROIs{r} '_10mm_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max']});

    for mm = 1:length(ROIs)
        
        map_name = ROIs{mm};
        
        unix(['mri_concat --i ' map_name '_regFrom_*.mgh --o ' ...
            map_name '_concat.mgh --mean']); 

    end
end