% Create control 5mm disk ROIs from center of each fROI for mSTS
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

%% Set up ROIs
faceROI = 'mSTS_faces';
fROIPre = 'fiberRFsclean_f_'; 

radius = 10;

%% Loop through hemis 
for h = 1:length(hems) 
    
    roiName = [fROIPre hems{h} '_' faceROI];
    
    %% Loop through ROIs and sessions
    for s = 1:length(sessions)
        session = sessions{s};

        sessionPath = [retDir '/' session];
        % move into the session dir
        cd(sessionPath)
        
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