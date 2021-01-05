% This restricts the average density map to EVC and normalizes
% (uses the Benson map EVC ROI for the mask)
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
fsaDir = fullfile(fsDir, 'fsaverage-bkup');
%% Set up ROIs
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

ROIPre = 'fibeRFs_f_'; 

%% Main
for h = 1:length(hems)
    
    ROIs={};
    for r = 1:length(faceROIs) %face ROIs
        if r == 5
            ROIs = horzcat(ROIs,{['fibeRFsclean_f_' hems{h} '_' faceROIs{r} '_10mm_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max_concat']});
        else
            ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{r} '_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max_concat']});
        end
    end
    ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi_r1.00_WholeBrainFG_track_' hems{h} '_proj_max_concat']}); %add CoS places
    EVC = {[hems{h} '.benson_EVC.label']};

    % Read in vertices from EVC label
    EVC_file=fullfile(fsaDir, 'label', EVC);
    EVC_label = read_label_kgs(EVC_file{1}); 
    EVC_label = sort(EVC_label);
    %use the vertices from the label (first column) to select only EVC vals
    idx = EVC_label(2:end,1);

    for ROI = 1:length(ROIs)

        roi_density = fullfile(fsaDir, 'surf',[ROIs{ROI} '.mgh']); 

        if exist(roi_density)
            [density_map, M, mr_parms, volsz] = load_mgh(roi_density); %get the density map
            for i = 1:length(density_map)
                if isempty(find(idx == i)) %if not an EVC idx
                    density_map(i)=0; %set to 0
                end
            end
            normed_density_map = density_map/(max(density_map)); %normalize by max value in map

            new_name = fullfile(fsaDir, 'surf',[ROIs{ROI} '_EVCtrim_normed.mgh']);
            save_mgh2(normed_density_map, new_name, M);
        end
    end
end