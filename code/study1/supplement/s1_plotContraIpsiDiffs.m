% pRF_plotContraIpsiDiffs.m
%
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
fieldRange = 40;
norm = 0;
voxthresh = 10;

%% Set up ROIs
allROIs = standardROIs;

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

for h = 1:length(hems)
    
    %% Let's load previously saved coverage
    dataFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100), '_' num2str(fieldRange) '_voxthresh' num2str(voxthresh) '_10mmcontrol']);
    load(dataFile);

    % average
    for r = V1:CoS
        average_coverage(r,:,:) = nanmean(fits(r).coverage,3);
    end

    % Now we will plot the coverage difference
    
    % do EVC separately
    f = figure('Position', [100, 100, 1000, 400], 'color', 'w');
    
    for r = V1:V3
        if strcmp(hems{h}, 'rh')
            contra = squeeze(average_coverage(r,:,1:64));
            ipsi = squeeze(average_coverage(r,:,65:128));
            diff = contra-fliplr(ipsi);
            total = [diff zeros(128,64)];
        else
            contra = squeeze(average_coverage(r,:,65:128));
            ipsi = squeeze(average_coverage(r,:,1:64));
            diff = fliplr(contra)-ipsi;
            total = [zeros(128,64) fliplr(diff)];
        end
        
        subplot_tight(1,3,r,[0.03,0.03]);
        s_createCoverageDifferencePlot(total,allROIs{r},fieldRange); colorbar off;
    end
    set(gcf, 'PaperPositionMode', 'auto');
    brighten(gcf,-0.2)
    saveFigFile = fullfile(figDir, [hems{h} '_EVC_contraIpsi_ve' num2str(ve_cutoff*100) '_' num2str(fieldRange) '.fig']); 
    print('-r300', '-dpng', fullfile(figDir, [hems{h} '_EVC_contraIpsi_ve' num2str(ve_cutoff*100) '_' num2str(fieldRange)]));
    saveas(gcf, saveFigFile)
    
    % functional ROIs
    f = figure('Position', [100, 100, 1500, 400], 'color', 'w');
    
    i = 1;
    for r = IOG:CoS
        if strcmp(hems{h}, 'rh')
            contra = squeeze(average_coverage(r,:,1:64));
            ipsi = squeeze(average_coverage(r,:,65:128));
            diff = contra-fliplr(ipsi);
            total = [diff zeros(128,64)];
        else
            contra = squeeze(average_coverage(r,:,65:128));
            ipsi = squeeze(average_coverage(r,:,1:64));
            diff = fliplr(contra)-ipsi;
            total = [zeros(128,64) fliplr(diff)];
        end
        
        subplot_tight(1,6,i,[0.03,0.03]);
        s_createCoverageDifferencePlot(total,allROIs{r},fieldRange); colorbar off;
        i = i+1;
    end
    set(gcf, 'PaperPositionMode', 'auto');
    brighten(gcf,-0.2)
    saveFigFile = fullfile(figDir, [hems{h} '_f_contraIpsi_ve', num2str(ve_cutoff*100), '_' num2str(fieldRange) '_10mmcontrol.fig']); 
    print('-r300', '-dpng', fullfile(figDir, [hems{h} '_f_contraIpsi_ve', num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol']));
    saveas(gcf, saveFigFile)
end
