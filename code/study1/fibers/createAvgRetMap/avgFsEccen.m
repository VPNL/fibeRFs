% This creates an average retinotopic map on the fs average surface
%
% Updated 12/2019 by DF
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
fsDir = fullfile(RAID, '3Danat/FreesurferSegmentations'); 
fsaDir = fullfile(fsDir, 'fsaverage-bkup', 'surf');

% parameter map
mapNames = {'retModel-testingDoubleSigsWithLin-cssFit-fFit'};

%% loop through sessions and transform maps to fsaverage surfaces using CBA
for h = 1:length(hems)
    
    map_name = [mapNames{1}, '_' hems{h} '_proj_max'];
    for ss = 1:length(fs_sessions)
        
        fs_id = fs_sessions{ss}; 

        % path to subject data in FreesurferSegmentations
        subjDir = fullfile(fsDir, fs_id);
        % paths to subject mri and surf directories
        mriDir = fullfile(subjDir, 'mri'); surfDir = fullfile(subjDir, 'surf');

        cd(surfDir);

        % transform surface files to fsaverage
        map_stem = fullfile(fsaDir, map_name);
        unix(['mri_surf2surf --srcsubject ' fs_id ' --srcsurfval ' ...
            map_name '.mgh --trgsubject fsaverage-bkup --trgsurfval ' ...
            map_stem '_regFrom_' fs_id '.mgh --hemi ' hems{h}]);

    end
end

%% average surface maps across all sessions
cd(fullfile(fsDir, 'fsaverage-bkup', 'surf'));
for h = 1:length(hems)
    
    map_name = [mapNames{1}, '_' hems{h} '_proj_max'];

    unix(['mri_concat --i ' map_name '_regFrom_*.mgh --o ' ...
        map_name '_mean_concat.mgh --mean']); 
end
