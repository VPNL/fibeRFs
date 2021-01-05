function s1_createTckmap(control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: 
% Creates trackmaps (mrTrix3) for each of the ROIs (i.e. determines the
% fibers which intersect with these ROIs)
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


        % write trk files
        for r = 1:length(runName)
            for n = 1:length(ROIs)
                ROIName = strsplit(ROIs{n}, '.');
                ROIfg = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG.mat'];
                fgname = [fullfile(dtiDir, dMRI_sessions{s}, runName{r}, 'dti96trilin', 'fibers', 'afq', ROIfg)];
                if exist(fgname)
                    disp(ROIName)
                    fg = load(fgname)
                    trkName = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG.tck'];
                    myfg = fg.roifg(1);
                    fpn = fullfile(dtiDir, dMRI_sessions{s}, runName{r}, 'dti96trilin/fibers/afq', trkName);
                    dr_fwWriteMrtrixTck(myfg, fpn)
                end
            end
        end

        for r = 1:length(runName)
            for n = 1:length(ROIs)
                ROIName = strsplit(ROIs{n}, '.');
                t1 = fullfile(dtiDir, dMRI_sessions{s}, '96dir_run1/t1/', t1_name);
                trkName = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG.tck'];
                fg = fullfile(dtiDir, dMRI_sessions{s}, runName{r}, 'dti96trilin/fibers/afq', trkName)
                if exist(fg)
                    disp(ROIName)
                    outname = [ROIName{1}, '_r', num2str(rad), '.00_WholeBrainFG_track.nii.gz'];
                    outdir = fullfile(dtiDir, dMRI_sessions{s}, '96dir_run1/t1/tracks')
                    mkdir(outdir)
                    output = fullfile(outdir, outname)
                    bkgrnd = false;
                    verbose = true;
                    mrtrixVersion = 3;

                    cmd_str = ['tckmap -template ', t1, ...
                        ' -ends_only -info', ...
                        ' -contrast tdi -force ', fg, ' ', output]

                    [status, results] = AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose, mrtrixVersion);
                end
            end
        end
    end
end

clear all
close all
