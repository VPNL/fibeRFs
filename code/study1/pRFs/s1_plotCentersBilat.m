% Plots the pRF centers for voxels across all subjects
% (also stores the average centers for subsequent use)
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
%exptDir = '/Volumes/kgs/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

sessionDir = [exptDir 'data/study1/toon'];
dataDir = [exptDir 'results/study1/pRFs'];
figDir = [exptDir 'results/study1/figs/manuscript'];

% params
ve_cutoff = .20;
fieldRange = 40;
fieldRange2plot = 30;
norm = 0;
thresh = 10; %just to pull in files

%% Set up ROIs
EVCROIs = standardROIs('EVC');
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

fROIPre = 'fiberRFsclean_f_'; 
retROIPre = 'fibeRFs_f_';

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

%% load prf sets

set(1) = load(fullfile(dataDir, ['rh_pRFset_' num2str(fieldRange)  '_ve' num2str(ve_cutoff*100) '_voxthresh' num2str(thresh) '_10mmcontrol.mat'])); 
set(2) = load(fullfile(dataDir, ['lh_pRFset_' num2str(fieldRange)  '_ve' num2str(ve_cutoff*100) '_voxthresh' num2str(thresh) '_10mmcontrol.mat'])); 

[~, xydeg_idx] = ismember('XYdeg', fieldnames(set(1).subj(1).roi(1).fits.vox)); %for reorg

%initialize
for h = 1:2
    for r = V1:CoS
        hemi(h).allRoi(r).vox = [];
    end
end

for h = 1:2
    for i = 1:length(sessions)
        for r = V1:CoS
            vox = [];
            %get average center per sub
            for v = 1:length(set(h).subj(i).roi(r).fits.vox)
                if set(h).subj(i).roi(r).fits.vox(v).eccen <= fieldRange2plot
                    vox = [vox set(h).subj(i).roi(r).fits.vox(v)];
                end
            end
            if ~isempty(vox)
                temp = struct2cell(vox.').'; XYdeg = temp(:, xydeg_idx); clear temp;
                XYdeg(cellfun(@(XYdeg) any(isnan(XYdeg)),XYdeg))=[]; %remove nans
                vals = cell2mat(XYdeg); 
                x_val(h,i,r) = nanmean(vals(:,1));
                y_val(h,i,r) = nanmean(vals(:,2));
            else
                x_val(h,i,r) = NaN;
                y_val(h,i,r) = NaN;
            end
            
            %collapse across subs for plotting
            for v = 1:length(set(h).subj(i).roi(r).fits.vox)
                if set(h).subj(i).roi(r).fits.vox(v).eccen <= fieldRange2plot
                    hemi(h).allRoi(r).vox = [hemi(h).allRoi(r).vox set(h).subj(i).roi(r).fits.vox(v)];
                end
            end
        end
    end
end

saveMatFile = fullfile(dataDir,['average_prf_centers_' num2str(fieldRange2plot) '_ve' num2str(ve_cutoff*100) '_10mmcontrol.mat']);
save(saveMatFile, 'x_val','y_val');
%% Now plot the centers


colors_left = [.7*[255/255, 0/255, 0/255];...
        .7*[255/255, 128/255, 0/255];...
        .7*[255/255, 255/255, 0/255];...
        .6*[0/255, 128/255, 255/255];...
        .4*[0/255, 0/255, 255/255];...
        .7*[0/255, 255/255, 0/255]];
colors_right = [[255/255, 0/255, 0/255];...
        [255/255, 128/255, 0/255];...
        [255/255, 255/255, 0/255];...
        [0/255, 128/255, 255/255];...
        .9*[0/255, 0/255, 255/255];...
        [0/255, 255/255, 0/255]];
    
f = figure('Position',[100 100 3000 600],'color','w');
vfc.fieldRange = fieldRange2plot;

i = 1;
for r = IOG:CoS
    
    subplot_tight(1,6,i,[0.02,0.02]);
    
    %set up polar plot
    p.fillColor='w';
    p.ringTicks = (1:3)/3*vfc.fieldRange;
    p.color = 'k';
    polarPlot([], p);
    caxis([0 1]);

    xlim([-vfc.fieldRange vfc.fieldRange])
    ylim([-vfc.fieldRange vfc.fieldRange])

    for rev_h = 1:2 %let's switch the ordering
        
        %lil switcheroo
        if rev_h == 1
            h = 2;
        else
            h = 1;
        end
        
        temp = struct2cell(hemi(h).allRoi(r).vox.').'; XYdeg = temp(:, xydeg_idx); clear temp;
        XYdeg(cellfun(@(XYdeg) any(isnan(XYdeg)),XYdeg))=[]; %remove nans
        vals = cell2mat(XYdeg); 

        
        if h == 1 %right hemi
            color = colors_right(i,:);
        elseif h == 2 %left hemi
            color = colors_left(i,:);
        end
        %if ~(r==mSTS&&h==2) %don't plot left mSTS
        scatter(vals(:,1),vals(:,2),15,color); hold on;
        %end
        
    end
    
    i=i+1; 
end

% Now save the figure
cd(figDir)

saveFigFile = fullfile(figDir,['all_pRF_centers_' num2str(fieldRange2plot) '_ve' num2str(ve_cutoff*100) '_10mmcontrol.fig']); 
print('-r300','-dpng',fullfile(figDir,['all_pRF_centers_' num2str(fieldRange2plot) '_ve' num2str(ve_cutoff*100) '_10mmcontrol'])) 

saveas(gcf,saveFigFile)

