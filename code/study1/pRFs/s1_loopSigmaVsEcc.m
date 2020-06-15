% s1_loopSigmaVsEcc.m
%
% This script will plot the ecc vs. pRF size for a list of ROIs and loop
% through a list of subjects, saving out each subject's figures into a
% figure directory. It will also save out the regression information for
% each ROI that a subject (slope, intercept, age, variance, etc).
%
% JG 05/2016
% edited by DF 06/2018
% more edits by DF 10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all

hems = {'rh' 'lh'};

% get our list of subjects from the Set function:
s1_setAllSessions
total_N = length(sessions);

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

sessionDir = [exptDir 'data/study1/toon'];
saveDir = [exptDir 'results/study1/pRFs'];
figDir = [exptDir 'results/study1/figs'];

%% Set up ROIs
EVCROIs = standardROIs('EVC');
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

fROIPre = 'fiberRFsclean_f_'; 
retROIPre = 'fibeRFs_f_';
allROIs = standardROIs;

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

%% Set threshold variables here:
vethresh = 0.10;
fieldRange = 40;
eccthresh = [0.5, fieldRange - 0.5]; %avoid edging effects
sigthresh = [0.21, Inf]; 
voxthresh = 25; %min number of voxels

% The structure in which we will store all subject data to save out:
lineData = {};

for h = 1:length(hems)
    
    %% get ROIs
    maps={};
    for r = 1:length(EVCROIs) % V1 through V3
        maps = horzcat(maps,{[retROIPre hems{h} '_' EVCROIs{r}]});
    end
    for r = 1:length(faceROIs) %face ROIs
        maps = horzcat(maps,{[fROIPre hems{h} '_' faceROIs{r}]});
    end
    roiList = horzcat(maps,{[fROIPre hems{h} '_' placeROIs{1}]}); %add CoS places
    
    saveFile = [hems{h} '_ve', num2str(vethresh*100), 'sigthresh', num2str(sigthresh(2)), 'eccen', num2str(round(eccthresh(2))), '_EccVsSigma_lineData'];

    for s = 1:length(sessions)

        cd(fullfile(sessionDir, sessions{s}))
        load mrSESSION.mat;

        %clear dataTYPES mrSESSION vANATOMYPATH;

        % Now load the appropriate datatype
        view = initHiddenGray('Averages', 1);
        rmPath = fullfile('Gray', 'Averages', 'retModel-testingDoubleSigsWithLin-cssFit-fFit.mat');
        view = rmSelect(view, 1, rmPath);

        fprintf('\n\n\nProcessing %s\n\n\n', sessions{s})

        for r = 1:length(roiList)

            if exist(fullfile(sessionDir, sessions{s}, '3Danatomy', 'ROIs', [roiList{r}, '.mat']), 'file') || exist(fullfile(sessionDir, sessions{s}, '3DAnatomy', 'ROIs', [roiList{r}, '.mat']), 'file')
                roiListNew = roiList{r};

                view = loadROI(view, roiListNew);

                list = 1:length(viewGet(view, 'ROIs'));
                data = s_rmPlotMultiEccSigma(view, list, vethresh, eccthresh, sigthresh, voxthresh);

                lineData{s} = data;

            end
        end

        mrvCleanWorkspace;
        close all;

    end

    save(fullfile(saveDir, saveFile), 'lineData')

end