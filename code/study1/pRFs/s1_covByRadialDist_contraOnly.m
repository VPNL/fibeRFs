% s1_covByRadialDist_contraOnly.m
%
% This script calculates the proportion of VFC by radial distance from the
% fovea for each ROI of interest and plots
%
% Only includes the contralateral visual field
%
% edited 10/31/19 DF
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
savePath = [exptDir 'results/study1/pRFs'];
figDir = [exptDir 'results/study1/figs/manuscript'];

% params
ve_cutoff = .10;
fieldRange = 40;
norm = 0;

%% Set up ROIs
allROIs = standardROIs;
%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

for h=1:length(hems)
    clear indivs condensed
    
    %% Let's load previously saved coverage
    dataFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100), '_' num2str(fieldRange)]);
    load(dataFile);
    
    %% Now plot
    
    % make internally consistent
    for i = 1:total_N
        for r = V1:CoS
            if isnan(fits(r).coverage(1,1,i))
                fits(r).includedvox(:, i) = NaN;
            end
        end
    end

    %also need to mask the bits outside the central circle for fairness!
    %(like in the coverage plots)
    mask = makecircle(128);

    rows = 128; columns = 128;
    fovea = 64;
    for i = 1:rows
        for j = 1:columns
            % Now let's find radial distance in DVA from center/fovea. Every 3.2
            % pixels is a degree of visual angle (128pixels/40 degree diameter)
            dva = 128/(fieldRange*2); 
            radialDist(i,j) = (sqrt((i-fovea)^2 + (j-fovea)^2))/dva;
        end
    end

    rounding = 1;  
    
    if strcmp(hems{h},'rh')
        radialDist = radialDist(:,1:64);
    else
        radialDist = radialDist(:,65:128);
    end

    %%
    for ROI = IOG:CoS
        clear inds

        N(ROI) =  sum(~isnan(squeeze(fits(ROI).coverage(1,1,:))));
        inds = find(~isnan(fits(ROI).coverage(1,1,:)));

        cov(ROI,:,:,1:N(ROI)) = fits(ROI).coverage(:,:,inds);

        flatDist = radialDist(:);

        for s = 1:N(ROI)
            
            a = squeeze(cov(ROI,:,:,s)).*mask; %mask bits outside of circle!
            
            %restrict coverage to contralateral visual field only
            if strcmp(hems{h},'rh')
                a = a(:,1:64);
            else
                a = a(:,65:128);
            end
            flatA = a(:);

            full = [flatDist, flatA];
            full = sortrows(full);
            edge = find(full(:,1)<fieldRange-0.5); 
            trimmed = full(edge,:);
            trimmed = round(trimmed,rounding); % round to help find out

            indivs = unique(trimmed(:,1));
            
            if ~exist('condensed', 'var')
                condensed = NaN(CoS, total_N, length(indivs));
            end

            for i = 1:length(indivs)
                inds = find(trimmed(:,1) == indivs(i));

                condensed(ROI,s,i) = mean(trimmed(inds,2));
            end

        end
    end


    colors = [[215/255, 48/255, 39/255];...
            [244/255, 109/255, 67/255];...
            [253/255, 174/255, 97/255];...
            [116/255, 173/255, 209/255];...
            [69/255, 117/255, 180/255];...
            [0/255, 153/255, 0/255]];

    %ROI indices
    ROI_names = {'IOG', 'pFus', 'mFus', 'pSTS', 'mSTS','CoS'}; 
    f = figure('Position',[100 100 3500 700]); hold on;
    
    i = 1;
    for ROI = IOG:CoS
        subplot_tight(1,6,i,[0.1, 0.03])
        shadedErrorBarOrig(indivs,nanmean(condensed(ROI,:,:)),(nanstd(condensed(ROI,:,:))/sqrt(N(ROI))),colors(i,:),'patchSaturation',0.6,'lineprops',{'-','color',colors(i,:),'linewidth',2});
        hold on

        if ROI == IOG
            set(gca,'ylim',[0 1],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0],'fontsize',12,'tickdir','out')
        else
           set(gca,'ylim',[0 1],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0],'fontsize',12,'tickdir','out')
           %ff = gca; ff.YAxis.Visible = 'off'; %(not working on 2015a)
        end   
        i = i+1;
    end
    set(gcf, 'PaperPositionMode', 'auto');
    saveFigFile = fullfile(figDir,[hems{h} '_radialDist_contraOnly_'  num2str(ve_cutoff*100) '_' num2str(fieldRange) '.fig']);
    print('-r300','-dpng',fullfile(figDir,[hems{h} '_radialDist_contraOnly_' num2str(ve_cutoff*100) '_' num2str(fieldRange)]))
    saveas(gcf,saveFigFile)
    
    %% Stats
    %let's start by fitting a line to all of them (by subject) and testing
    %slopes
    
    for s = 1:total_N %for the 15 subjects
        for ROI = IOG:CoS
            if s <= N(ROI)
                mdl = fitlm(indivs, squeeze(condensed(ROI,s,:))); %fit linear regression
                Rsquared(ROI-3,s) = mdl.Rsquared.Ordinary;
                RSME(ROI-3,s) = mdl.RMSE;
                intercept(ROI-3,s) = mdl.Coefficients.Estimate(1);
                slope(ROI-3,s) = mdl.Coefficients.Estimate(2);
            else 
                Rsquared(ROI-3,s) = NaN;
                RSME(ROI-3,s) = NaN;
                intercept(ROI-3,s) = NaN;          
                slope(ROI-3,s) = NaN;
            end
        end
    end

    saveMatFile = fullfile(savePath,[hems{h} '_radialFits_linear_contraOnly_' num2str(ve_cutoff*100) '_' num2str(fieldRange) '.mat']);
    save(saveMatFile, 'intercept','slope', 'Rsquared', 'RSME');
end 

close all