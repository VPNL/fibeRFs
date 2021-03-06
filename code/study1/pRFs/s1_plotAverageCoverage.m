% s1_plotAveragedCoverage.m
%
% This script will loop through subject retinotopy data and
% will then plot averaged coverage maps
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
exptDir = fullfile(RAID,expt);

sessionDir = [exptDir 'data/study1/pRFs'];
savePath = [exptDir 'results/study1/pRFs'];
figDir = [exptDir 'results/study1/figs/manuscript'];

% params
ve_cutoff = .20;
fieldRange = 20; %30
norm = 0;
thresh=10;

%% Set up ROIs
allROIs = standardROIs;

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

%% Loop through hemis and plot coverage
for h = 1:length(hems)

    %% Let's load previously saved coverage
    dataFile = fullfile(savePath, [hems{h} '_coverage_data_ve' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_voxthresh' num2str(thresh) '_10mmcontrol']);
    load(dataFile);
    
    avg_centers = fullfile(savePath, ['average_prf_centers_' num2str(fieldRange) '_ve' num2str(ve_cutoff*100) '_10mmcontrol']);
    load(avg_centers)
    %% Now plot
    
    % make internally consistent
    for i = 1:length(sessions)
        for r = IOG:CoS
            if isnan(fits(r).coverage(1,1,i))
                fits(r).includedvox(:, i) = NaN;
            end
        end
    end
    
    % Now we will plot the coverages!

    % functional ROIs
    f = figure('Position', [100, 100, 1500, 600], 'color', 'w');
    
    i = 1;
    for r = IOG:CoS
        subplot_tight(1, 6, i, [0.03, 0.03]);
        s_createCoveragePlot_averaged(nanmean(fits(r).coverage,3), allROIs{r}, fieldRange, norm, num2str(sum(~isnan(fits(r).includedvox))));
        scatter(nanmean(x_val(h,:,r)),nanmean(y_val(h,:,r)),50,'w','*'); 
        colorbar off;
        i=i+1; 
    end
    
    set(gcf, 'PaperPositionMode', 'auto');
    brighten(gcf,-0.2)
    saveFigFile = fullfile(figDir, [hems{h} '_f_averageCoverage_ve' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_voxthresh' num2str(thresh) '_10mmcontrol.fig']); 
    print('-r300', '-dpng', fullfile(figDir, [hems{h} '_f_averageCoverage_ve' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_voxthresh' num2str(thresh) '_10mmcontrol']));
    saveas(gcf, saveFigFile)


end

close all

