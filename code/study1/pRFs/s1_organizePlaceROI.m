% unifyROIs - combines (for fibeRFs study) ROIs across sessions in a single
% naming convention, eliminates overlap, saves information about original
% studies

% adapted from sonia's roiPrep.m to finalize ROIs for an experiment. it will
% A) rename ROIs from different sessions/experiments to a standard preFix,
% which avoids overwriting ROIs that you share with others. it does this by
% looking for a softlink in each experimental session folder called
% 'faceLoc' or 'Retinotopy'
% B) exclude your EVC ROIs or any other set
% C) exclude overlapping voxels within this set in a prescribed order
% D) save a summary struct of info for each subject, number of vox in each ROI, which loc, and number of loc runs

% 10/31/19 adapted from SP's code by DF
% transfers and cleans the previous latprf CoS places ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill in this info:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

expt = '/projects/fibeRFs/data/study1/loc'; 
exptDir = fullfile(RAID,expt);

s1_setAllSessions; %loads sessions for this study

%%%% which ROIs
hems = {'rh' 'lh'};
ROIs = standardROIs('place');
colors = {'g'};

%%%% exclusion paramenters
exclude = standardROIs('EVC');
removeOverlap = 0; % remove overlap between the faceROIs themselves (separate from exclude EVC); 0 = no, -1 = reverse standardROIs order (recommended), 1 = standardROIs order
%removeOverlap currently only works if you have all ROIS

%%%% naming convention of output
newPre = 'fiberRFsclean_f_';     % preFix to add to the ROI name - usually expt_f/a_.

%%%%% comment field information - session & anat will be auto-filled
comment.name = 'Dawn Finzi';
comment.note = 'Oct 2020';
outStr = []; % aggregate output information, display at end (since mrV output will disrupt readability)

%%%%% output info
outputMat = [exptDir '/placeROIinfo.mat'];
outputTxt = [exptDir '/placeROIinfo.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start processing ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% housekeeping
if removeOverlap == -1 ROIs = fliplr(ROIs); colors = fliplr(colors); end % clean in reverse-heirarchy order

for s = 1:length(sessions)
    
% initialize info/comment formatting for this subject    
session = sessions{s}; 
cd(fullfile(exptDir, session));
info{s}.missing = []; % track which ROIs we don't locate
info{s}.numVox = zeros(length(hems),length(ROIs));
outStr = [outStr sprintf('\n***%s***\n',sessions{s})];
comment.sess = pwd;

% initialize hidden gray
vw = initHiddenGray('Averages');

for h = 1:length(hems) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create a single exclusion ROI for this hemisphere
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 1:length(exclude) vw = loadROI(vw,['fibeRFs_f_' hems{h} '_' exclude{n}],[],[],[],0); end
    vw = combineROIs(vw,1:length(vw.ROIs), 'union', 'exclude','k');
    vw = deleteROI(vw,1:length(vw.ROIs)-1); % exclude is vw.ROIs(1) for the rest of this processing
    
    missing = 0; %how many fROIs is this subj missing
    
    for r = 1:length(ROIs) 
        
        % input and output ROI names for this hemisphere
        inROI = ['latprf_f_' hems{h},'_' ROIs{r}];
        backupROI = ['fibeRFs_f_' hems{h},'_' ROIs{r}];
        outROI = strcat(newPre,hems{h},'_',ROIs{r});
        info{s}.outROIs{h,r} = outROI;
        
        % check if they exist in the current loc shared path
        if (~exist(fullfile('3DAnatomy','ROIs',[inROI '.mat'])) && ~exist(fullfile('3Danatomy','ROIs',[inROI '.mat']))) %case issues
            if (~exist(fullfile('3DAnatomy','ROIs',[backupROI '.mat'])) && ~exist(fullfile('3Danatomy','ROIs',[backupROI '.mat']))) %case issues
                fprintf('Missing %s...\n',inROI);
                info{s}.missing = [info{s}.missing inROI];
                missing = missing + 1;
                outStr = [outStr sprintf('Missing %s %s...\n',sessions{s},inROI)];
            else
                inROI = backupROI;
            end
        end
        
        if (exist(fullfile('3DAnatomy','ROIs',[inROI '.mat'])) || exist(fullfile('3Danatomy','ROIs',[inROI '.mat']))) %if determined backup or main naming ROI exists

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scrape info about localizer & add comments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % how many runs?
            if exist([pwd '/Stimuli']) pars = dir([pwd '/Stimuli/parfiles/*par']); else pars = dir([pwd '/stim/parfiles/*par']); end
            comment.numRuns = length(pars);
            
            % which localizer was this?
            lookFor = {'fLoc' 'CategoryChannels' 'kidLoc'};
            for n = 1:length(lookFor)
                if containsTxt(pars(1).name,lookFor{n}) comment.loc = lookFor{n}; end
            end
            
            % aggregate info into ROI comment field in mrVista, and our
            % info struct
            comment.txt = sprintf('Session: %s\nExperiment: %s\nNum runs: %d\nOriginal name: %s\nAuto-excluding: %s\nDrawn by: %s\nNote: %s\n',...
                comment.sess,comment.loc,comment.numRuns,inROI,strTogether(exclude),comment.name,comment.note);
            if ~(isfield(info{s}, 'subj')) %have we written the basics yet
                info{s}.subj = sessions{s};
                info{s}.comment = comment;
                info{s}.loc = comment.loc;
                info{s}.numRuns = comment.numRuns;
            end
            
            % load shared/input-ROI
            vw = loadROI(vw,inROI,[],[],[],0);
            
            % save ROI with new name, comment field (overwrites existing)
            vw.ROIs(end).name = outROI;
            vw.ROIs(end).color = colors{r};
            vw.ROIs(end).comment = comment.txt;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % exclude EVC (or other set of regions)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % exclude our 'exclusion' ROI
            vw = combineROIs(vw, {vw.ROIs(end).name 'exclude'}, 'AnotB', outROI, colors{r}, [comment.txt]);
            
            % now save this cleaned and renamed ROI in shared
            [vw, ~, ~] = saveROI(vw,length(vw.ROIs), 0, 1);    % 1 = forceSave
            info{s}.numVox(h,r)= size(vw.ROIs(end).coords,2);
            outStr = sprintf('%sExcluded %s and saved %s as %s...\n',outStr,strTogether(exclude),inROI,outROI);
            vw = deleteROI(vw,2:length(vw.ROIs)); % deletes everything except for 'exclude' ROI
        end
    end
    
    vw = deleteROI(vw,1:length(vw.ROIs));
    
end

mrvCleanWorkspace;
end

fprintf(['\n\n\n' outStr]);
save(outputMat,'info');

% text pRF info
ROIs = reshape(info{1}.outROIs',1,numel(info{1}.outROIs));

fid = fopen(outputTxt,'w');
fprintf(fid,[sprintf('subj\tloc\tnumRuns') sprintf('\t%s',ROIs{:})]); % header

data  = [];
for s = 1:length(info)
    data = sprintf('%s\n%s\t%s\t%d',data,info{s}.subj,info{s}.loc,info{s}.numRuns);
    voxNum = reshape(info{s}.numVox',1,numel(info{1}.outROIs));
    data=sprintf('%s%s',data,sprintf('\t%d',voxNum(:)));
end

fprintf(fid,data);
cd(dirOf(outputTxt));
