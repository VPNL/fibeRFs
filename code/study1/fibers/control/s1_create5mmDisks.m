% Create control 5mm disk ROIs from center of each fROI
%
% Adapted from code by MN
% Updated 12/2019 by DF
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

radius = 5;

%% Loop through hemis 
for h = 1:length(hems) 
    
    %% Get ROIs
    ROIs={};
    for r = 1:length(faceROIs) %face ROIs
        ROIs = horzcat(ROIs,{[fROIPre hems{h} '_' faceROIs{r}]});
    end
    ROIs = horzcat(ROIs,{[fROIPre hems{h} '_' placeROIs{1}]}); %add CoS places
    

    %% Loop through ROIs and sessions
    for s = 1:length(dMRI_sessions)
        session = dMRI_sessions{s};

        sessionPath = [retDir '/' session];
        % move into the session dir
        cd(sessionPath)
        
        for r=1:length(ROIs)

            roiName = ROIs{r};

            % init hidden gray and load ROI
            hG = initHiddenGray('GLMs',1,roiName);

            % CHeck if ROI exists
            if isempty(hG.ROIs)
                sprintf('%s ROI not found for session %s', roiName, session)
            else 
                %% Create both a center ROI 
                diskNameStr = sprintf('%dmm', radius);
                diskROIname = [roiName '_' diskNameStr];
                [hG, newDiskROI, layers] = makeROIdiskGray(hG, radius, diskROIname, [], [], mean(hG.ROIs.coords, 2));
                [vw, status, forceSave] = saveROI(hG, newDiskROI);

            end

        close all;
        clear hG;
        clear global 
        end
    end
end