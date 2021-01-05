% s1_covByRadialDist_comparingFits.m
%
% This script fits either a line or a sigmoid curve to the coverage by
% radial distance measurements for each subject and ROI separately
%
% Saves out params and goodness of fit measures
%
% 11/20 DF
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
ve_cutoff = .20; 
fieldRange = 40;% 40;
thresh = 10;
norm = 0;

%% Set up ROIs
allROIs = standardROIs;
%corresponding indices - same throughout
V1 = 1; V2 = 2; V3 = 3; 
IOG = 4; pFus = 5; mFus = 6; pSTS = 7; mSTS = 8; CoS = 9;

for h=1:length(hems)
    
    %% Let's load previously saved radial calcs
    dataFile = fullfile(savePath,[hems{h} '_radial_raw_contraOnly_' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol.mat']);
    load(dataFile);
    
    %% Stats
    
    sig = fittype('a + (b-a) ./ (1 + 10.^((c-x)*d))','coeff',{'a','b','c','d'});
    initial_params = [0,1,5,.1]; 
    
    for s = 1:total_N 
        for ROI = IOG:CoS
            if ~isnan(condensed(ROI,s,1))
                %linear fit
                mdl = fitlm(indivs, squeeze(condensed(ROI,s,:))); %fit linear regression'Lower', [0, 0, 0, 0],
                Rsquared(ROI-3,s) = mdl.Rsquared.Ordinary;
                Adj_Rsquared(ROI-3,s) = mdl.Rsquared.Adjusted;
                RSME(ROI-3,s) = mdl.RMSE;
                intercept(ROI-3,s) = mdl.Coefficients.Estimate(1);
                slope(ROI-3,s) = mdl.Coefficients.Estimate(2);
                
                %sigmoid function
                [curve, goodness] = fit(indivs, squeeze(condensed(ROI,s,:)), sig, 'Lower', [0, -1, 0, 0], 'Upper', [2, 1, 250, 1], 'StartPoint', initial_params);
                if ~isinf(goodness.adjrsquare) %only if it actually converged
                    Sig_Rsquared(ROI-3,s) = goodness.rsquare;
                    Sig_Adj_Rsquared(ROI-3,s) = goodness.adjrsquare;
                    Sig_RSME(ROI-3,s) = goodness.rmse;
                    Sig_a(ROI-3,s) = curve.a;
                    Sig_b(ROI-3,s) = curve.b;
                    Sig_c(ROI-3,s) = curve.c;
                    Sig_d(ROI-3,s) = curve.d;
                else
                    Sig_Rsquared(ROI-3,s) = NaN;
                    Sig_Adj_Rsquared(ROI-3,s) = NaN;
                    Sig_RSME(ROI-3,s) = NaN;
                    Sig_a(ROI-3,s) = NaN;
                    Sig_b(ROI-3,s) = NaN;
                    Sig_c(ROI-3,s) = NaN;
                    Sig_d(ROI-3,s) = NaN;
                end

            else 
                %linear fit
                Rsquared(ROI-3,s) = NaN;
                Adj_Rsquared(ROI-3,s) = NaN;
                RSME(ROI-3,s) = NaN;
                intercept(ROI-3,s) = NaN;          
                slope(ROI-3,s) = NaN;
                
                %sigmoid function
                Sig_Rsquared(ROI-3,s) = NaN;
                Sig_Adj_Rsquared(ROI-3,s) = NaN;
                Sig_RSME(ROI-3,s) = NaN;
                Sig_a(ROI-3,s) = NaN;
                Sig_b(ROI-3,s) = NaN;
                Sig_c(ROI-3,s) = NaN;
                Sig_d(ROI-3,s) = NaN;

            end
        end
    end
    
    saveMatFile = fullfile(savePath,[hems{h} '_radialFits_linear_contraOnly_' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol.mat']);
    save(saveMatFile, 'intercept','slope', 'Rsquared', 'Adj_Rsquared', 'RSME');
    
    saveMatFile = fullfile(savePath,[hems{h} '_radialFits_sigmoid_contraOnly_' num2str(ve_cutoff*100) '_' num2str(fieldRange) '_10mmcontrol.mat']);
    save(saveMatFile, 'Sig_a', 'Sig_b', 'Sig_c', 'Sig_d', 'Sig_Rsquared', 'Sig_Adj_Rsquared', 'RSME');
    
    %% summary adjusted r-squared stats
    linear_mean_by_subj = nanmean(Adj_Rsquared);
    linear_mean(h) = nanmean(linear_mean_by_subj);
    linear_ste(h) = nanstd(linear_mean_by_subj)/sqrt(sum(~isnan(linear_mean_by_subj)));
    
    sig_mean_by_subj = nanmean(Sig_Adj_Rsquared);
    sig_mean(h) = nanmean(sig_mean_by_subj);
    sig_ste(h) = nanstd(sig_mean_by_subj)/sqrt(sum(~isnan(linear_mean_by_subj)));
end

