% Calculates proportion of fibers for each face ROI that connect to EVC
%
% Updated 11/2019 by DF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% get our list of subjects from the Set function:
s1_setAllSessions

hems = {'rh' 'lh'};

% where do the subjects live
expt = '/projects/fibeRFs/'; 
exptDir = fullfile(RAID,expt);

dtiDir = [exptDir 'data/study1/diffusion'];
outDir = [exptDir 'results/study1'];

%% Set up ROIs
faceROIs = standardROIs('face');
placeROIs = standardROIs('place');

ROIPre = 'fibeRFs_f_'; 

cd(outDir);

runName={'96dir_run1'};
r=1; 

%% Fiber count for each ROI
for h = 1:length(hems)
    
    input_ROI = [hems{h}, '_EVC'];
    
    ROIs={};
    for roi = 1:length(faceROIs) %face ROIs
        ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' faceROIs{roi} '_projed_gmwmi']});
    end
    ROIs = horzcat(ROIs,{[ROIPre hems{h} '_' placeROIs{1} '_projed_gmwmi']}); %add CoS places

    for n=1:length(ROIs)

        %loads in pariwise
        name_ending=['_pairwise'];
        filename=strcat(ROIs{n}, '_', input_ROI, name_ending,'.mat');
        load(fullfile(outDir,filename))
        
        %check subjcount
        num_subj=length(fibercount);
        subjidx=fibercount(1:num_subj,1);
        
        pair_fibercount = fibercount;
        
        % pull data from non-pairwise
        name_ending=['_fibers_in_ROI'];
        filename=strcat(ROIs{n}, name_ending,'.mat');
        load(fullfile(outDir,filename))
        
        %pull data from pairwise
        for i=1:num_subj
            
            fibercount_all_for_DC(i,1)=subjidx(i,1);
            fibercount_all_for_DC(i,3)=pair_fibercount(i,2); %num of pairwise fibers
            
            idx=find(fibercount(:,1)==subjidx(i,1));
            fibercount_all_for_DC(i,2)=fibercount(idx,2); %total fibers for ROI

            perROI(i,1)=fibercount(i,1);
            perROI(i,2)=fibercount_all_for_DC(i,3)/fibercount_all_for_DC(i,2);
        end

        perROI(i+1,2)=nanmean(perROI(1:i,2));
        perROI(i+2,2)=nanstd(perROI(1:i,2)/sqrt(i));

        outName=[ROIs{n} '_' input_ROI '_perFaceROI'];
        save(fullfile(outDir, outName),'perROI')
        
        overall(:,n) = perROI(1:num_subj,2);

        clear fibercount_all_for_DC
        clear perROI
        
    end
    
    outName=[hems{h} '_fibers_overall_per_ROI'];
    save(fullfile(outDir, outName),'overall')
    
    clear overall


end
