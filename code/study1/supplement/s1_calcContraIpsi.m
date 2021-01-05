% s1_calcContraIpsi.m
%
% This script will loop through the pRF coverage data and
% will calculate the contra field bias for each participant and ROI
% where CVI = mean contra - mean ipsi
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
ve_cutoff = .20;
fieldRange = 30;
norm = 0;
voxthresh = 10;

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

for h = 1:length(hems)
    
    %% Let's load previously saved coverage
    dataFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100), '_' num2str(fieldRange) '_voxthresh' num2str(voxthresh) '_10mmcontrol']);
    load(dataFile);

    %calculate laterality index from coverage
    for r = V1:CoS
        for s = 1:length(fits(r).coverage(1, 1,:))
            
            if strcmp(hems{h}, 'rh')
                contra(r,s,:,:)=squeeze(fits(r).coverage(:,1:64,s));
                ipsi(r,s,:,:)=squeeze(fits(r).coverage(:,65:128,s));
            else
                contra(r,s,:,:)=squeeze(fits(r).coverage(:,65:128,s));
                ipsi(r,s,:,:)=squeeze(fits(r).coverage(:,1:64,s));
            end

            mean_contra(r,s) = mean(mean(contra(r,s,:,:)));
            mean_ipsi(r,s) = mean(mean(ipsi(r,s,:,:)));
            total = (mean_contra(r,s) + mean_ipsi(r,s));

            lat_index(r,s) = (mean_contra(r,s) - mean_ipsi(r,s))/total;
        end
    end

    %% compile
    both_lat_index(h,:,:) = [lat_index(1,:); lat_index(2,:); lat_index(3,:); lat_index(4,:);...
        lat_index(5,:); lat_index(6,:); lat_index(7,:); lat_index(8,:); lat_index(9,:)];
    li(h,:) = nanmean(both_lat_index(h,:,:),3); 
    liErr(h,:) = nanstd(squeeze(both_lat_index(h,:,:))',1) ./ sqrt(sum(~isnan(squeeze(both_lat_index(h,:,:))'),1));

end

%% Plotting 
vals = [li(2,:); li(1,:)];
errs = [liErr(2,:); liErr(1,:)];
width = .85;
groupnames = {'V1', 'V2', 'V3' 'IOG', 'pFus', 'mFus', 'pSTS', 'mSTS', 'CoS'};
bw_title = 'Contra vs. Ipsi Index';
bw_xlabel = 'ROI';
bw_ylabel = 'Laterality';
bw_colormap = jet;
gridstatus = 'none';
bw_legend = {'Left', 'Right'};
error_sides = [];
legend_type = 'plot';

barweb(vals', errs', width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)

set(gcf, 'PaperPositionMode', 'auto');
saveFigFile = fullfile(figDir,['ContraVSIpsi' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol.fig']);
print('-r300','-dpng',fullfile(figDir,['ContraVSIpsi' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol']))
saveas(gcf,saveFigFile)

right_lat_index = squeeze(both_lat_index(1,:,:));
left_lat_index = squeeze(both_lat_index(2,:,:));

saveMatFile = fullfile(savePath,['ContraVSIpsi_' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol.mat']);
save(saveMatFile, 'right_lat_index', 'left_lat_index'); 
