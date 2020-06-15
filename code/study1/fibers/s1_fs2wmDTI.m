function s1_fs2wmDTI(control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: 
% Preparing the mrVista ROIs for use in FDWM
% Transforms them from the freesurfer space to the dti space - extending
% them into the white matter so they can be intersected with fibers later
%
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Loop through hemis 
for h = 1:length(hems)
    
    %% Get ROIs
    if control == 1 %5mm control analysis
        maps={};
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{[ROIPre hems{h} '_' faceROIs{r} '_5mm']});
        end
        maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_5mm']}); %add CoS places
    else
        maps={};
        maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_EVC']});
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' faceROIs{r}]});
        end
        maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1}]}); %add CoS places
    end
     
    %% now extend the gray matter labels into gmwmi
    for s = 1:length(dMRI_sessions)
        
        session = fullfile(retDir, dMRI_sessions{s});
        fs_subject = fs_sessions{s};
        labelPath = fullfile('/biac2/kgs/3Danat/FreesurferSegmentations/', fs_subject, '/label/');

        dataPath = fullfile(session, '3DAnatomy', 'niftiROIs');

        for r = 1:length(maps)

            %dilate the labels and save them as niftis as well as .mat
            label = strcat(labelPath, maps{r}, '_surf.label');
            nifti = strcat(dataPath, '/', maps{r}, '_projed');

            ref = fullfile(dtiDir, dMRI_sessions{s}, runName{1}, 'dti96trilin', 'mrtrix', 'dwi_aligned_trilin_noMEC_gmwmi.nii.gz');
            % NOTE: you need to create 'dwi_aligned_trilin_noMEC_gmwmi.nii.gz'
            % run "mrconvert dwi_aligned_trilin_noMEC_gmwmi.mif dwi_aligned_trilin_noMEC_gmwmi.nii.gz" and save the nifti with your rois

            cd(fullfile(dataPath))

            if exist(label) > 0

                direction = -1; %extend into white, not gray matter
                fs_labelFileToNiftiRoiGWMI(fs_subject, label, nifti, hems{h}, ref, 0, [], direction) %this dilates the label using freesurfer

                nii = strcat(dataPath, '/', maps{r}, '_projed.nii.gz');
                gmwm = ref;

                load(strcat(dtiDir, '/', dMRI_sessions{s}, '/', runName{1}, '/dti96trilin/ROIs/ATR_roi1_R.mat')) %loads an roi so we have a template for the different mat fields

                gmwmi = niftiRead(gmwm); % this crops the dilated rois to the gmwmi created by mrtrix
                imgROI = niftiRead(nii);
                imgBoth = niftiRead(nii);
                imgBoth.data(:, :, :) = 0;
                idx1 = (find(gmwmi.data(:, :, :) > 0));
                idx2 = (find(imgROI.data(:, :, :) > 0));
                idx = intersect(idx1, idx2);

                imgBoth.data(idx) = 1;

                name = strcat(maps{r}, '_projed_gmwmi.nii.gz'); %saves the new roi as nii
                niftiWrite(imgBoth, name);

                imgCoords = find(imgBoth.data);
                [I, J, K] = ind2sub(imgBoth.dim, imgCoords);
                roi.coords = mrAnatXformCoords(imgBoth.qto_xyz, [I, J, K]);

                cd(fullfile(dtiDir, '/', dMRI_sessions{s}, '/', runName{1}, '/dti96trilin/ROIs/'));
                outname = strcat(maps{r}, '_projed_gmwmi.mat')
                save(outname, 'roi', 'versionNum', 'coordinateSpace'); %also saves it as .mat dti roi
                niftiWrite(imgBoth, name); %save the nifti ROIs in the dti folder too
            end
        end
    end
end

clear all
close all
