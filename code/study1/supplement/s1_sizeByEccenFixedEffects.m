% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: changes from lsline to fitline2derror
%
% updated for fibeRFs by DF 2020

clear all;
close all

% get our list of subjects from the Set function:
s1_setAllSessions

hems = {'rh', 'lh'};

% where do the subjects live
expt = '/projects/fibeRFs/'; 
%exptDir = fullfile(RAID,expt);
exptDir = '/Volumes/kgs/projects/fibeRFs/';

sessionDir = [exptDir 'data/study1/toon'];
dataDir = [exptDir 'results/study1/pRFs'];
figDir = [exptDir 'results/study1/figs/manuscript'];

ve_cutoff = .20; %.10; 
fieldRange = 40;
thresh = 10;

xmax = 40;
ymax = 40;

%% Set up ROIs
ROIs_txt = {'V1' 'V2' 'V3' 'IOG faces' 'pFus faces' 'mFus faces' 'pSTS faces ' 'mSTS faces' 'CoS places'};
ROIs= [standardROIs('EVC') standardROIs('face') standardROIs('place')];

%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

%% Which size formula to use?
size_form = 2; 
% 1 = standard sigma, same as coverage
% 2 = mrvista option of sigma/sqrt(exp)
% 3 = what kendrick uses (2*sigma/sqrt(exp))

%% Load prf sets
prfset(1) = load(fullfile(dataDir, ['rh_pRFset_' num2str(fieldRange)  '_ve' num2str(ve_cutoff*100) '_voxthresh' num2str(thresh) '.mat'])); 
prfset(2) = load(fullfile(dataDir, ['lh_pRFset_' num2str(fieldRange)  '_ve' num2str(ve_cutoff*100) '_voxthresh' num2str(thresh) '.mat'])); 

%initialize
for h = 1:length(hems)
    for r = V1:CoS
        hemi(h).allRoi(r).vox = [];
    end
end

%collapse across subs for plotting
for h = 1:length(hems)
    for i = 1:length(sessions)
        for r = V1:CoS
             hemi(h).allRoi(r).vox = [hemi(h).allRoi(r).vox prfset(h).subj(i).roi(r).fits.vox];
        end
    end
end

%% Plot
sampleVox = 500; 
saveFig = 1;

fontSize = 11; titleSize = 14;

for h = 1:length(hems)
    
    f(1) = niceFig([.1 .1 .4 .8],fontSize,2);

    for r = 1:length(ROIs)
        figure(f(1));
        subplot(3,3,r);

        ROInum = r; %cellNum(ROIs{r},info.ROIs);
        fits = hemi(h).allRoi(ROInum); % roi(ROInum).fits;
        
        eccen_reshape = [fits.vox.eccen];
        eccen_trim = eccen_reshape(eccen_reshape <= (fieldRange - 0.5)); %remove edge artifacts
        eccen = eccen_trim(~isnan(eccen_trim));
        
        if size_form == 1
            size = double([fits.vox.size]);
        elseif size_form == 2
            size = double([fits.vox.mrvsize]);
        elseif size_form == 3
            size = double([fits.vox.kksize]);
        end
        size = size(eccen_reshape <= (fieldRange - 0.5));
        size = size(~isnan(size));
        
        
        % use all vox if less than sampleVox num of voxels
        if sampleVox < 1 sv = round(length(eccen)*sampleVox); %less than eccen b/c removed NaNs
            elseif sampleVox > length(eccen) sv = length(eccen);
            else sv = sampleVox; 
        end
        sv = datasample([1:length(eccen)],sv,'Replace',false);

        eccen_subset = eccen(sv);
        size_subset = size(sv);
        
        s = scatter(eccen_subset,size_subset,10,roiColors(ROIs{r}),'o'); hold on;
        s.MarkerEdgeAlpha = 0.25;
        
        %only fit the line where we have voxels
        if ROInum >= IOG && ROInum <= pFus
            fitRange = [0:.5:10]; 
        elseif ROInum == mFus
            fitRange = [0:.5:7.5];
        else
            fitRange = [0:.5:xmax];
        end
        
        [N,edges] = histcounts(eccen,fitRange);
        [errorObj,lineObj] = scatterline(eccen,size,fitRange,NaN,500,roiColors(ROIs{r}),0,1); hold on;
        extX = [0:.5:xmax];
        extLine = plot(extX, extendLine(lineObj,extX),'Color',roiColors(ROIs{r}),'LineWidth',1.5);
        set(lineObj,'LineWidth',1.5); set(extLine,'LineStyle',':');
        alpha(errorObj,.75);

        xl = xlim; xlim([0 xmax]); yl = ylim; ylim([0 ymax]);
        set(gca,'TickDir','out');
        xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (dva)'],'FontSize',fontSize); 
        axis square; 
        title(ROIs_txt{r});

    end

    if saveFig == 1
        if size_form == 1
            txt = ['sizeEccen_r'  num2str(ve_cutoff*100) '_' hems{h}];
        elseif size_form == 2
            txt = ['sizeEccen_r'  num2str(ve_cutoff*100) '_mrvSize_' hems{h}];
        elseif size_form == 3
            txt = ['sizeEccen_r'  num2str(ve_cutoff*100) '_kkSize_' hems{h}];
        end
        figure(f(1));
        niceSave(figDir,['/scatter_square_4040_' txt],[],[],{'png' 'svg'});
    end
end
