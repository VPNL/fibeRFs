% This calculates the number of fibers per ROI and then the number of
% fibers between each face ROI and EVC
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% get our list of subjects from the Set function:
s1_setAllSessions

hems = {'rh' 'lh'};

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

dtiDir = [exptDir 'data/study1/diffusion'];
outDir = [exptDir 'results/study1/fibers/10mm'];
fsDir = fullfile(RAID, '3Danat/FreesurferSegmentations'); 
fsaDir = fullfile(fsDir, 'fsaverage-bkup', 'surf');
%% Set up ROIs
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

ROIPre = 'fibeRFs_f_'; 

cd(outDir);

runName={'96dir_run1'}
r=1; 

%% Fiber count for each ROI
for h = 1:length(hems)

    ROIs={};
    for roi = 1:length(faceROIs) %face ROIs
        if roi == 5
            ROIs = horzcat(ROIs,{['fibeRFsclean_f_' hems{h} '_' faceROIs{roi} '_projed_gmwmi']});
        else
            ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{roi} '_projed_gmwmi']});
        end
    end
    ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi']}); %add CoS places
    ROIs = horzcat(ROIs, {[ROIPre hems{h} '_EVC_projed_gmwmi']});

    for n=1:length(ROIs)
        fibercount=[];
        subjcnt=1;
        for s=1:length(dMRI_sessions)  % Ok, here we go
            ROIName=ROIs{n};
            ROIfgname=[ROIName '_r1.00_WholeBrainFG.mat'];
            roifgFile = fullfile(dtiDir,dMRI_sessions{s},runName{r},'dti96trilin','fibers','afq',ROIfgname);
            if exist(roifgFile,'file')
                roifg=load(roifgFile);

                fibercount(s,1) = s;
                fibercount(s,2) = size(roifg.roifg.fibers,1);

                subjcnt=subjcnt+1;
            else
                fibercount(s,1) = s;
                fibercount(s,2) = NaN;
                subjcnt=subjcnt+1;
            end

        end

        outName = [ROIName '_fibers_in_ROI'];
        save(fullfile(outDir,outName),'fibercount')
    end
end

%% Fiber count between each face ROI and EVC

for h = 1:length(hems)
    
    input_ROI = [ROIPre hems{h}, '_EVC_projed_gmwmi'];
    
    ROIs={};
    for roi = 1:length(faceROIs) %face ROIs
        if roi == 5
            ROIs = horzcat(ROIs,{['fibeRFsclean_f_' hems{h} '_' faceROIs{roi} '_projed_gmwmi']});
        else
            ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{roi} '_projed_gmwmi']});
        end
    end
    ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi']}); %add CoS places


    for n=1:length(ROIs)
        fibercount=[];
        subjcnt=1;
        
        for s=1:length(dMRI_sessions)  % Ok, here we go
            ROIfgname=[ROIs{n} '_r1.00_' input_ROI '_r1.00_WholeBrainFG.mat'];
            roifgFile = fullfile(dtiDir,dMRI_sessions{s},runName{r},'dti96trilin','fibers','afq',ROIfgname);
            if exist(roifgFile,'file')
                roifg=load(roifgFile);

                fibercount(s,1) = s;
                fibercount(s,2) = size(roifg.roifg.fibers,1);

                subjcnt=subjcnt+1;
            else
                fibercount(s,1) = s;
                fibercount(s,2) = NaN;
                subjcnt=subjcnt+1;
            end
        end

        fixEVCName = [hems{h} '_EVC'];

        outName = [ROIs{n} '_' fixEVCName '_pairwise'];
        save(fullfile(outDir,outName),'fibercount')
    end

end