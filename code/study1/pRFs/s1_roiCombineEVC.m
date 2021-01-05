% Combine the V1-V3 ROIs drawn from toonotopy into EVC ROIs

clear all

expt = 'fibeRFs/data/study1/toon';
exptDir = fullfile('/share/kalanit/biac2/kgs/projects/', expt);

s1_setAllSessions

mapName = 'retModel-cssFit-fFit';

hems = {'rh', 'lh'};
ROIs = {'V1', 'V2', 'V3'};
currPre = 'fibeRFs_f_';

for s = 1:length(sessions)

    cd(fullfile(exptDir, sessions{s}));

    for h = 1:length(hems)

        inROIs = strcat(currPre, hems{h}, '_', ROIs);
        nameEVC = strcat(currPre, hems{h}, '_EVC');

        vw = initHiddenGray('Averages'); %open a view

        % load and combine ROIs
        for r = 1:length(inROIs)
            vw = loadROI(vw, inROIs{r});
        end

        vw = combineROIs(vw, {vw.ROIs(1:end-1).name, vw.ROIs(end).name}, ...
            'union', nameEVC, 'b', 'Combined V1-V3 via roiCombineEVC.m');

        % now save the combined ROI in shared dir
        [vw, status, forceSave] = saveROI(vw, 'selected', 0, 1);
        outStr = sprintf('Successfully made joint EVC!\n');
    end

    fprintf(['\n', outStr]);

    mrvCleanWorkspace;
end
