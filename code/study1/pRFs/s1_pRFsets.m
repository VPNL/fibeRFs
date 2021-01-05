% Gather all pRF model data to mimic Sonia's pRFsets plotting code
%
% This also stores all the data on individual voxel pRFs centers
%
% Updated 10/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all

% get our list of subjects from the Set function:
s1_setAllSessions

hems = {'rh' 'lh'};

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

sessionDir = [exptDir 'data/study1/toon'];
dataDir = [exptDir 'results/study1/pRFs'];

% params
ve_cutoff = .20; 
fieldRange = 20; %40
sigThresh = 0.21;
eThresh = [0 fieldRange]; 
norm = 0; 
voxThresh = 10; %minimum number of voxels per ROI

%% Set up ROIs
EVCROIs = standardROIs('EVC');
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

fROIPre = 'fiberRFsclean_f_'; 
retROIPre = 'fibeRFs_f_';

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;


%% Loop through hemis and calc centers
for h = 1:length(hems)
    
    clear info subj
    
    %% Set up the pRF sets format
    % 1. info struct
    info.subjs = sessions;
    info.expt = 'toon';
    info.minR2 = ve_cutoff*100;
    info.ROIs = maps;
    info.hems = hems;
    % not sure if these are necc...
    info.task = '';
    info.whichStim = 'sweep';
    info.whichModel ='kayCSS';
    info.setNotes = '';
    info.fitSuffix = ''; 
    
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

       cd(fullfile(sessionDir,sessions{i}))
       load mrSESSION.mat;

       clear dataTYPES mrSESSION vANATOMYPATH;

       % Now load the appropriate datatype
       view = initHiddenGray('Averages',1);
       rmPath = fullfile('Gray','Averages','retModel-testingDoubleSigsWithLin-cssFit-fFit.mat');
       view = rmSelect(view,1,rmPath);

        for m = 1:length(maps) % if sub has map, put pRF_cov into data, otherwise fill with NaNs
            clear pRF_COV coIndices figHandle all_models weight data; 

            subj(i).roi(m).fits.session = sessions{i};
            subj(i).roi(m).fits.ROIname = maps{m};
            
            if exist(fullfile('/share', 'kalanit', 'biac2', 'kgs', 'projects', 'fibeRFs', 'data', 'study1', 'toon', sessions{i}, '3DAnatomy', 'ROIs', [maps{m}, '.mat'])) || exist(fullfile('/share', 'kalanit', 'biac2', 'kgs', 'projects', 'fibeRFs', 'data', 'study1', 'toon', sessions{i}, '3Danatomy', 'ROIs', [maps{m}, '.mat']))

                view = loadROI(view,maps{m});
                % easy way to get all the data
                [pRF_COV, coIndices, figHandle, all_models, weight, data] = rmPlotCoverage_flips(view,'prf_size',1,'fieldRange',fieldRange,'method','maximum profile','nboot',50,'normalizeRange',0,'smoothSigma',1,'cothresh',ve_cutoff,'weight','variance explained','sigmathresh',sigThresh,'eccthresh',eThresh,'addcenters',1,'newfig',0, 'flips', flips(i), 'css', 1);

                if ~isempty(data) %only if there are above threshold voxels
                    if length(data.subSize1) < voxThresh %and more of them than our set threshold
                        subj(i).roi(m).fits(1).vox(1).size = NaN; %prf size for that voxel
                        subj(i).roi(m).fits(1).vox(v).sigma = NaN;
                        subj(i).roi(m).fits(1).vox(1).XYdeg = NaN; %prf center
                        subj(i).roi(m).fits(1).vox(1).eccen = NaN; % eccentricity
                        subj(i).roi(m).fits(1).vox(1).r2 = NaN; % variance explained
                    else
                        for v = 1:length(data.subSize1)
                            if data.subSize1(v) > sigThresh
                                subj(i).roi(m).fits(1).vox(v).size = data.subSigSize(v); %prf size for that voxel - already takes into account diving by sqrt(expt) in rmGet.m
                                subj(i).roi(m).fits(1).vox(v).sigma = data.subSize1(v); %this one is actually just sigma
                                subj(i).roi(m).fits(1).vox(v).XYdeg = [data.subx0(v), data.suby0(v)]; %prf center
                                subj(i).roi(m).fits(1).vox(v).eccen = data.subEcc(v); % eccentricity
                                subj(i).roi(m).fits(1).vox(v).r2 = data.subCo(v); % variance explained
                            else
                                subj(i).roi(m).fits(1).vox(v).size = NaN;
                                subj(i).roi(m).fits(1).vox(v).sigma = NaN;
                                subj(i).roi(m).fits(1).vox(v).XYdeg = NaN;
                                subj(i).roi(m).fits(1).vox(v).eccen = NaN;
                                subj(i).roi(m).fits(1).vox(v).r2 = NaN;
                            end
                        end
                    end
                else
                    subj(i).roi(m).fits(1).vox(1).size = NaN; %prf size for that voxel
                    subj(i).roi(m).fits(1).vox(v).sigma = NaN;
                    subj(i).roi(m).fits(1).vox(1).XYdeg = NaN; %prf center
                    subj(i).roi(m).fits(1).vox(1).eccen = NaN; % eccentricity
                    subj(i).roi(m).fits(1).vox(1).r2 = NaN; % variance explained
                end    

            else
                subj(i).roi(m).fits(1).vox(1).size = NaN; %prf size for that voxel
                subj(i).roi(m).fits(1).vox(v).sigma = NaN;
                subj(i).roi(m).fits(1).vox(1).XYdeg = NaN; %prf center
                subj(i).roi(m).fits(1).vox(1).eccen = NaN; % eccentricity
                subj(i).roi(m).fits(1).vox(1).r2 = NaN; % variance explained
            end
        end

        mrvCleanWorkspace; 

    end

    saveFile = fullfile(dataDir,[hems{h} '_pRFset_' num2str(fieldRange) '_ve' num2str(ve_cutoff*100) '_voxthresh' num2str(voxThresh) '_10mmcontrol']); 
    save(saveFile,'info', 'subj');
end

close all

