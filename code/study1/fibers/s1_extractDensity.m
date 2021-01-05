function s1_extractDensity(control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: 
% Extract the endpoint density values from the surface trackmaps
%
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if notDefined('control')
    control = 0;
end

% get our list of subjects from the Set function:
s1_setAllSessions

hems = {'rh' 'lh'};

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

retDir = [exptDir 'data/study1/toon'];
dtiDir = [exptDir 'data/study1/diffusion'];
fsDir = fullfile(RAID, '3Danat/FreesurferSegmentations'); 
if control == 0
    resultsDir = [exptDir 'results/study1/fibers'];
elseif control == 1
    resultsDir = [exptDir 'results/study1/fibers/control'];
elseif control == 2
    resultsDir = [exptDir 'results/study1/fibers/control/10mm'];
end


bins = 40;
num_bands = 4;
%% Set up ROIs
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

ROIPre = 'fibeRFsclean_f_'; 

%corresponding indices - same throughout
IOG = 1; pFus = 2; mFus = 3; pSTS = 4; mSTS = 5; CoS = 6;
%% Loop through hemis 
for h = 1:length(hems)  
    
    %% Get ROIs
    maps={};
    if control == 1
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{[ROIPre hems{h} '_' faceROIs{r} '_5mm_projed_gmwmi']});
        end
        ROIs = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_5mm_projed_gmwmi']}); %add CoS places
    elseif control == 2
        for r = 1:length(faceROIs) %face ROIs
            if r == 5 %mSTS
                maps = horzcat(maps,{[ROIPre hems{h} '_' faceROIs{r} '_10mm_projed_gmwmi']});
            else
                maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' faceROIs{r} '_projed_gmwmi']});
            end
        end
        ROIs = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_projed_gmwmi']}); %add CoS places
    else
        for r = 1:length(faceROIs) %face ROIs
            maps = horzcat(maps,{['fibeRFs_f_' hems{h} '_' faceROIs{r} '_projed_gmwmi']});
        end
        ROIs = horzcat(maps,{['fibeRFs_f_' hems{h} '_' placeROIs{1} '_projed_gmwmi']}); %add CoS places
    end
    
    EVC = ['fibeRFs_f_' hems{h} '_EVC_surf.label']; %fs label
    % retinotopic map
    ret = {['retModel-testingDoubleSigsWithLin-cssFit-fFit_' hems{h} '_proj_max.mgh']};
    
    %% initialize
    medians = [];
    means = [];

    weighted_medians = [];
    weighted_means = [];

    weighted_hist = [];
    weighted_hist_cdf = [];
    
    %% Loop subjects
    for s = 1:length(fs_sessions);

        dir = fullfile(fsDir, fs_sessions(s)); %subj fs dir

        EVC_file = fullfile(dir, 'label', EVC);
        EVC_label = read_label_kgs(EVC_file{1});
        EVC_label = sort(EVC_label);
        %use the vertices from the label (first column) to select only EVC vals
        %idx = EVC_label(2:end,1);
        EVC_lind = EVC_label(:, 1);
        idx = EVC_lind(EVC_lind ~= 0);

        eccen = fullfile(dir, 'surf', ret);
        eccen_map = load_mgh(eccen{1});

        EVC_eccen = eccen_map(idx);

        for ROI = IOG:CoS
            clear vals total band wvals wtotal wband
            
            roi_density = fullfile(dir, 'surf', [ROIs{ROI}, '_r1.00_WholeBrainFG_track_', hems{h}, '_proj_max.mgh']);

            if exist(roi_density{1})
                density_map = load_mgh(roi_density{1});
                EVC_density = density_map(idx);

                combined = [EVC_eccen, EVC_density];

                any_fibers = combined(combined(:, 2) ~= 0, :);
                any_nonzero_fibs = any_fibers(any_fibers(:, 1) ~= 0, :);
                
                if length(any_nonzero_fibs(:,1)) > 1 %need more than one fiber endpoint
                    medians(ROI, s) = nanmedian(any_nonzero_fibs(:, 1));
                    means(ROI, s) = nanmean(any_nonzero_fibs(:, 1));

                    %% weighted
                    x = sortrows(any_nonzero_fibs, 1);
                    x(:, 2) = x(:, 2) * 1000;
                    x(:, 2) = round(x(:, 2));

                    wvals = [];
                    for i = 1:length(x)
                        wvals = vertcat(wvals, repmat(x(i, 1), x(i, 2), 1));
                    end
                    weighted_medians(ROI, s) = nanmedian(wvals);
                    weighted_means(ROI, s) = nanmean(wvals);

                    %% bands analysis
                    vals = any_nonzero_fibs(:, 1);
                    total = length(vals);
                    wtotal = length(wvals);
                    for b = 1:num_bands
                        if b == 1
                            band(b).vals = vals(vals > 0 & vals <= 5); %0 to 5 deg
                            band(b).wvals = wvals(wvals > 0 & wvals <= 5);
                        elseif b == 2
                            band(b).vals = vals(vals > 5 & vals <= 10); %5 to 10 deg
                            band(b).wvals = wvals(wvals > 5 & wvals <= 10);
                        elseif b == 3
                            band(b).vals = vals(vals > 10 & vals <= 20); %10 to 20 deg
                            band(b).wvals = wvals(wvals > 10 & wvals <= 20);
                        elseif b == 4
                            band(b).vals = vals(vals > 20);
                            band(b).wvals = wvals(wvals > 20);
                        end
                        props(ROI, b, s) = length(band(b).vals)/total;
                        wprops(ROI, b, s) = length(band(b).wvals) / wtotal;
                    end


                    %% weighted hists
                    hi = histogram(wvals, bins, 'Normalization', 'cdf', 'BinLimits', [0, 40]);
                    hi.FaceColor = 'c';
                    ylim([0, 0.2]);
                    weighted_hist_cdf(ROI, s, :) = hi.Values;

                    p = histogram(wvals, bins, 'Normalization', 'pdf', 'BinLimits', [0, 40]);
                    ylim([0, 0.2]);
                    weighted_hist(ROI, s, :) = p.Values;
                else %only one fiber endpoint edgecase
                    medians(ROI, s) = NaN;
                    means(ROI, s) = NaN;
                    props(ROI, 1:num_bands, s) = NaN;
                    weighted_medians(ROI, s) = NaN;
                    weighted_means(ROI, s) = NaN;
                    wprops(ROI, 1:num_bands, s) = NaN;
                    weighted_hist_cdf(ROI, s, :) = NaN;
                    weighted_hist(ROI, s, :) = NaN;
                end
            
            else %no ROI for that subj
                medians(ROI, s) = NaN;
                means(ROI, s) = NaN;
                props(ROI, 1:num_bands, s) = NaN;
                weighted_medians(ROI, s) = NaN;
                weighted_means(ROI, s) = NaN;
                wprops(ROI, 1:num_bands, s) = NaN;
                weighted_hist_cdf(ROI, s, :) = NaN;
                weighted_hist(ROI, s, :) = NaN;
            end
        end
    end

    cd(resultsDir)
    
    save([hems{h} '_' num2str(num_bands) '_bands'], 'props');
    save([hems{h} '_weighted_' num2str(num_bands) '_bands'], 'wprops');
    save([hems{h} '_medians'], 'medians', 'weighted_medians');
    save([hems{h} '_means'], 'means', 'weighted_means');
    save([hems{h} '_hists'], 'weighted_hist_cdf', 'weighted_hist');

end

clear all
close all
