% s1_calcUpperLower_contraOnly.m
%
% This script will loop through the pRF coverage data and
% will calculate the upper field bias for each participant and ROI
% (but only for the contralateral visual field)
% where UVB = mean upper - mean lower
%
% DF 07/2018
% DF 11/2019
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
figDir = [exptDir 'results/study1/figs'];

% params
ve_cutoff = .10;
fieldRange = 40;
norm = 0;

%% Set up ROIs
allROIs = standardROIs;

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;


for h = 1:length(hems)
    
    %% Let's load previously saved coverage
    dataFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100), '_' num2str(fieldRange)]);
    load(dataFile);

    %calculate UVF bias from coverage
    for r = V1:CoS
        for s = 1:length(fits(r).coverage(1, 1,:))

            if strcmp(hems{h}, 'rh')
                upper(r,s,:,:)=squeeze(fits(r).coverage(65:128,1:64,s));
                lower(r,s,:,:)=squeeze(fits(r).coverage(1:64,1:64,s));
            else
                upper(r,s,:,:)=squeeze(fits(r).coverage(65:128,65:128,s));
                lower(r,s,:,:)=squeeze(fits(r).coverage(1:64,65:128,s));
            end

            mean_upper(r,s) = mean(mean(upper(r,s,:,:)));
            mean_lower(r,s) = mean(mean(lower(r,s,:,:)));
            total = (mean_upper(r,s) + mean_lower(r,s));

            uvl_index(r,s) = (mean_upper(r,s) - mean_lower(r,s))/total;
        end
    end

    %% compile
    both_uvl_index(h,:,:) = [uvl_index(1,:); uvl_index(2,:); uvl_index(3,:); uvl_index(4,:);...
        uvl_index(5,:); uvl_index(6,:); uvl_index(7,:); uvl_index(8,:); uvl_index(9,:)];
    li(h,:) = nanmean(both_uvl_index(h,:,:),3); 
    liErr(h,:) = nanstd(squeeze(both_uvl_index(h,:,:))',1) ./ sqrt(sum(~isnan(squeeze(both_uvl_index(h,:,:))'),1));

end

%% Plotting 
figure;
vals = [li(2,:); li(1,:)];
errs = [liErr(2,:); liErr(1,:)];
width = .85;
groupnames = {'V1', 'V2', 'V3', 'IOG', 'pFus', 'mFus', 'pSTS', 'mSTS', 'CoS'};
bw_title = 'Upper vs. Lower Index';
bw_xlabel = 'ROI';
bw_ylabel = 'UVF - LVF';
bw_colormap = jet;
gridstatus = 'none';
bw_legend = {'Left', 'Right'};
error_sides = [];
legend_type = 'plot';

barweb(vals', errs', width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
ylim([-0.4 0.4])

set(gcf, 'PaperPositionMode', 'auto');
saveFigFile = fullfile(figDir,['UpperVSLower_contraOnly_'  num2str(ve_cutoff*100) '_' num2str(fieldRange) '.fig']);
print('-r300','-dpng',fullfile(figDir,['UpperVSLower_contraOnly_'   num2str(ve_cutoff*100) '_' num2str(fieldRange)]))
saveas(gcf,saveFigFile)

right_lat_index = squeeze(both_uvl_index(1,:,:));
left_lat_index = squeeze(both_uvl_index(2,:,:));

saveMatFile = fullfile(savePath,['UpperVSLower_contraOnly_'   num2str(ve_cutoff*100) '_' num2str(fieldRange) '.mat']);
save(saveMatFile, 'right_lat_index', 'left_lat_index'); 
