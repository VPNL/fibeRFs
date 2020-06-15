% s1_compareAverageEccSigma
%
% This code will load data files saved out form pRF_loopSigmaVsEcc and the
% _VTC version of it as well in order to compare average sigma and average
% eccentricity for each ROI across groups.
%
% JG 06/2016
% DF 09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all

hems = {'rh' 'lh'};

% params
eThresh = 40; 
vThresh = 0.10;
trim = 0;
siginput = Inf; %to pull the line data
sigthresh = Inf;

% load subject names:
s1_setAllSessions;

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

sessionDir = [exptDir 'data/study1/toon'];
dataDir = [exptDir 'results/study1/pRFs'];
figDir = [exptDir 'results/study1/figs'];

roiNames = {'V1', 'IOG', 'pFus', 'mFus', 'pSTS', 'mSTS', 'CoS'};

% Set up data structure, rows per subject, columns per ROI
V1 = 1;
IOG = 2;
pFus = 3;
mFus = 4;
pSTS = 5;
mSTS = 6;
CoS = 7;

% initialize
ecc = NaN(length(hems), length(sessions), length(roiNames));
sig = NaN(length(hems), length(sessions), length(roiNames));

band = NaN(length(hems), length(sessions), length(roiNames), 3); %divided into 3 eccen bands
control_band = NaN(length(hems), length(sessions), length(roiNames), 4); %divided into 4 eccen bands

sizeByBand = NaN(length(hems), length(sessions), length(roiNames), 3); %divided into 3 eccen bands
control_sizeByBand = NaN(length(hems), length(sessions), length(roiNames), 4); %divided into 4 eccen bands

varByBand = NaN(length(hems), length(sessions), length(roiNames), 3); %divided into 3 eccen bands
control_varByBand = NaN(length(hems), length(sessions), length(roiNames), 4); %divided into 4 eccen bands

for h = 1:length(hems)

    % Which file would you like to visualize?
    fileName = [hems{h}, '_ve', num2str(vThresh*100), 'sigthresh', num2str(siginput), 'eccen', num2str(eThresh), '_EccVsSigma_lineData'];

    load(fullfile(dataDir, fileName));

    %% Now extract  roi information
    for i = 1:numel(lineData)

        for m = 1:numel(lineData{1, i})

            % Check if subject has any ROIs
            if isempty(lineData{1, i})
                continue
            else
                roi_name = lineData{1, i}(1, m).roi;
            end

            if containsTxt(roi_name, 'V1')
                r = V1;
            elseif containsTxt(roi_name, 'IOG')
                r = IOG;
            elseif containsTxt(roi_name, 'pFus')
                r = pFus;
            elseif containsTxt(roi_name, 'mFus')
                r = mFus;
            elseif containsTxt(roi_name, 'pSTS')
                r = pSTS;
            elseif containsTxt(roi_name, 'mSTS')
                r = mSTS;
            elseif containsTxt(roi_name, 'CoS')
                r = CoS;
            end

            variance = lineData{1, i}(1, m).variance;
            eccent = lineData{1, i}(1, m).ecc;
            sigma = lineData{1, i}(1, m).sigma;

            %remove eccen vals where variance < threshold or sigma less
            %than model limits (0.21)
            eccent(variance <= vThresh) = NaN;
            eccent(eccent > eThresh) = NaN;
            eccent(sigma < 0.21) = NaN;
            eccent(sigma > sigthresh) = NaN;

            %now remove these voxels from corresponding sigma and variance
            %variables
            sigma(isnan(eccent)) = NaN;
            variance(isnan(eccent)) = NaN;

            ecc(h, i, r) = nanmedian(eccent);
            sig(h, i, r) = nanmedian(sigma);

            eccent = eccent(~isnan(eccent));
            sigma = sigma(~isnan(sigma));
            variance = variance(~isnan(variance));
            
            if ~isnan(nanmedian(eccent))    
                for b = 1:4
                    if b == 1
                        e = length(find(eccent <= 5)) / length(eccent);
                        s = nanmedian(sigma(eccent <= 5));
                        v = nanmedian(variance(eccent <= 5));
                    elseif b == 2
                        e = length(find(eccent > 5 & eccent <= 10)) / length(eccent);
                        s = nanmedian(sigma(eccent > 5 & eccent <= 10));
                        v = nanmedian(variance(eccent > 5 & eccent <= 10));
                    elseif b ==3 
                        e = length(find(eccent > 10 & eccent <= 20)) / length(eccent);
                        s = nanmedian(sigma(eccent > 10 & eccent <= 20));
                        v = nanmedian(variance(eccent > 10 & eccent <= 20));
                    elseif b ==4
                        e = length(find(eccent > 20)) / length(eccent);
                        s = nanmedian(sigma(eccent > 20));
                        v = nanmedian(variance(eccent > 20));
                    end
                    
                    if b <= 3
                        band(h, i, r, b) = e;
                        control_band(h, i, r, b) = e;
                        
                        sizeByBand(h, i, r, b) = s;
                        control_sizeByBand(h, i, r, b) = s;
                        
                        varByBand(h, i, r, b) = v;
                        control_varByBand(h, i, r, b) = v;
                    else
                        control_band(h, i, r, b) = e;
                        control_sizeByBand(h, i, r, b) = s;
                        control_varByBand(h, i, r, b) = v;
                    end
                end
            end

        end
    end

    %% Reorder band info and save
    for roi = V1:CoS
        eccenByBand(roi,:,:) = [band(h,:,roi ,1); band(h,:,roi,2); band(h,:,roi,3)]';
        eccenByBandControl(roi,:,:) = [control_band(h,:,roi ,1); control_band(h,:,roi,2);...
            control_band(h,:,roi,3); control_band(h,:,roi,4)]';
        
        sigByBand(roi,:,:) = [sizeByBand(h,:,roi ,1); sizeByBand(h,:,roi,2); sizeByBand(h,:,roi,3)]';
        sigByBandControl(roi,:,:) = [control_sizeByBand(h,:,roi ,1); control_sizeByBand(h,:,roi,2);...
            control_sizeByBand(h,:,roi,3); control_sizeByBand(h,:,roi,4)]';
        
        veByBand(roi,:,:) = [varByBand(h,:,roi ,1); varByBand(h,:,roi,2); varByBand(h,:,roi,3)]';
        veByBandControl(roi,:,:) = [control_varByBand(h,:,roi ,1); control_varByBand(h,:,roi,2);...
            control_varByBand(h,:,roi,3); control_varByBand(h,:,roi,4)]';
    end
    
    saveMatFile = fullfile(dataDir, [hems{h}, '_eccenBands_' num2str(eThresh) '_ve', num2str(vThresh*100) '.mat']);
    save(saveMatFile, 'eccenByBand');
    saveMatFile = fullfile(dataDir, [hems{h}, '_eccenBandsControl_' num2str(eThresh) '_ve', num2str(vThresh*100) '.mat']);
    save(saveMatFile, 'eccenByBandControl');
    
    saveMatFile = fullfile(dataDir, [hems{h} '_sigBands_' num2str(eThresh) '_ve', num2str(vThresh*100) '.mat']);
    save(saveMatFile, 'sigByBand');
    saveMatFile = fullfile(dataDir, [hems{h} '_sigBandsControl_' num2str(eThresh) '_ve', num2str(vThresh*100) '.mat']);
    save(saveMatFile, 'sigByBandControl');
    
    saveMatFile = fullfile(dataDir, [hems{h} '_varBands_' num2str(eThresh) '_ve', num2str(vThresh*100) '.mat']);
    save(saveMatFile, 'veByBand');
    saveMatFile = fullfile(dataDir, [hems{h} '_varBandsControl_' num2str(eThresh) '_ve', num2str(vThresh*100) '.mat']);
    save(saveMatFile, 'veByBandControl');


    %% Now let's plot (not for pub - visualization purposes)

    % Let's make a 2 panel plot with pRF size on the top, eccentricity on the
    % bottom.
    Ecc = nanmedian(squeeze(ecc(h, :, :)), 1);
    EccSEM = nanstd(squeeze(ecc(h, :, :)), 1) ./ sqrt(sum(~isnan(squeeze(ecc(h, :, :))), 1));

    Sig = nanmedian(squeeze(sig(h, :, :)), 1);
    SigSEM = nanstd(squeeze(sig(h, :, :)), 1) ./ sqrt(sum(~isnan(squeeze(sig(h, :, :))), 1));

    f = figure('Position', [100, 100, 1400, 800]); %set(gca,'ylim',[0 30]);

    if strcmp(hems{h}, 'rh')
        suptitle('Right Hemisphere', 'fontsize', 18)
    else
        suptitle('Left Hemisphere', 'fontsize', 18)
    end

    subplot_tight(2, 1, 1, [0.08, 0.08])
    x = [1:1:size(Sig, 2)];

    errorbar(Sig, SigSEM, 'ko', 'Color', [0.5273, 0.8047, 0.9792], 'LineWidth', 4);
    hold on;
    scatter(x, Sig, 75, 'MarkerFaceColor', [0.5273, 0.8047, 0.9792], 'MarkerEdgeColor', 'b');

    set(gca, 'xtick', [1:1:7], 'xticklabel', roiNames, 'FontSize', 16) %,'ylim',[0 30])
    ylabel('Median pRF size (deg.)', 'FontSize', 18)

    % Now let's plot mean eccentricity of the pRFs in a given ROI
    subplot_tight(2, 1, 2, [0.04, 0.08])
    errorbar(Ecc, EccSEM, 'ko', 'Color', [0.8555, 0.4375, 0.5742], 'LineWidth', 4); hold on
    scatter(x, Ecc, 75, 'MarkerFaceColor', [0.8555, 0.4375, 0.5742], 'MarkerEdgeColor', 'r');

    set(gca, 'xtick', [1:1:7], 'xticklabel', roiNames, 'FontSize', 16)
    ylabel('Median Eccentricity (deg.)', 'FontSize', 18)

    set(gcf, 'PaperPositionMode', 'auto');
    saveFile = [hems{h}, '_pRF_median_size_and_eccentricity_upTo', num2str(eThresh), '_sigthresh', num2str(sigthresh), '_ve', num2str(vThresh*100), '.fig'];
    saveas(gcf, fullfile(figDir, saveFile))
    print('-r300', '-dpng', fullfile(figDir, [hems{h}, '_pRF_median_size_and_eccentricity_upTo', num2str(eThresh), '_sigthresh', num2str(sigthresh), '_ve', num2str(vThresh*100)]))
    
    %save ecc and sig medians
    saveMatFile = fullfile(dataDir, [hems{h}, '_pRF_median_size_and_eccentricity_upTo', num2str(eThresh) ,'_ve', num2str(vThresh*100), '.mat']);
    fix_ecc = squeeze(ecc(h,:,:))';
    fix_sig = squeeze(sig(h,:,:))';
    save(saveMatFile, 'fix_ecc', 'fix_sig')
end
