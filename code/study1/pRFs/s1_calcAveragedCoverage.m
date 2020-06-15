% s1_calcAveragedCoverage.m
%
% This script will loop through subject retinotopy data and
% will then create averaged coverage maps
%
% Adapted from JG 05/2016
% DF 07/2018
% DF 10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all

hems = {'rh' 'lh'};

% get our list of subjects from the Set function:
s1_setAllSessions

% where do the subjects live
expt = '/projects/fibeRFs/'; 
%exptDir = '/Volumes/kgs/projects/fibeRFs/';
exptDir = fullfile(RAID,expt);

sessionDir = [exptDir 'data/study1/toon'];
savePath = [exptDir 'results/study1/pRFs'];

% params for rmPlotCoverage_flips function
ve_cutoff = .10;
fieldRange = 40; 
sigThresh = 0.21;
eThresh = [0, fieldRange];
norm = 0;
thresh = 25; %threshold for number of voxels needed to include

%% Set up ROIs
EVCROIs = standardROIs('EVC');
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

fROIPre = 'fiberRFsclean_f_'; 
retROIPre = 'fibeRFs_f_';

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

%% Loop through hemis and calc coverage
for h = 1:length(hems)
    
    %% Calculate map coverage
    maps={};
    for r = 1:length(EVCROIs) % V1 through V3
        maps = horzcat(maps,{[retROIPre hems{h} '_' EVCROIs{r}]});
    end
    for r = 1:length(faceROIs) %face ROIs
        maps = horzcat(maps,{[fROIPre hems{h} '_' faceROIs{r}]});
    end
    maps = horzcat(maps,{[fROIPre hems{h} '_' placeROIs{1}]}); %add CoS places

    for i = 1:length(sessions)

        cd(fullfile(sessionDir, sessions{i}))
        load mrSESSION.mat;

        clear dataTYPES mrSESSION vANATOMYPATH;

        % Now load the appropriate datatype
        view = initHiddenGray('Averages', 1);
        rmPath = fullfile('Gray', 'Averages', 'retModel-testingDoubleSigsWithLin-cssFit-fFit.mat');
        view = rmSelect(view, 1, rmPath);

        for m = 1:length(maps) % if sub has map, put pRF_cov into data, otherwise fill with NaNs
            
            clear pRF_COV coIndices;
            
            if exist(fullfile('/share', 'kalanit', 'biac2', 'kgs', 'projects', 'fibeRFs', 'data', 'study1', 'toon', sessions{i}, '3DAnatomy', 'ROIs', [maps{m}, '.mat'])) || exist(fullfile('/share', 'kalanit', 'biac2', 'kgs', 'projects', 'fibeRFs', 'data', 'study1', 'toon', sessions{i}, '3Danatomy', 'ROIs', [maps{m}, '.mat']))
                
                view = loadROI(view, maps{m});
                [pRF_COV, coIndices] = rmPlotCoverage_flips(view, 'prf_size', 1, 'fieldRange', fieldRange, 'method', 'density', 'nboot', 50, 'normalizeRange', 0, 'smoothSigma', 1, 'cothresh', ve_cutoff, 'weight', 'variance explained', 'sigmathresh', sigThresh, 'eccthresh', eThresh, 'addcenters', 1, 'newfig', -1, 'flips', flips(i), 'css', 1);
                
                fits(m).includedvox(:, i) = sum(coIndices);
                
                %only included coverage data for ROIs with 10 vox or more
                if fits(m).includedvox(:, i) >= thresh
                    fits(m).coverage(:, :, i) = pRF_COV;
                else 
                    fits(m).coverage(:, :, i) = NaN(128, 128);
                end
                
                fits(m).totalvox(:, i) = numel(coIndices);
                fits(m).ROIname = maps{m};
                
            else
                
                fits(m).coverage(:, :, i) = NaN(128, 128);
                fits(m).includedvox(:, i) = NaN;
                fits(m).totalvox(:, i) = NaN;
                
            end
        end

        mrvCleanWorkspace;

    end

    %Save the coverage info 
    saveFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100) '_' num2str(fieldRange)]);
    save(saveFile, 'fits');
end
