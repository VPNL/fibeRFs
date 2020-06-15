function fatVistaRoi2DtiRoi(fatDir, sessid, runName, roiName, t1_name, type)
% fatVistaRoi2DtiRoi(fatDir, sessid, runName, roiName)
% runName, a cell array for run name
% roiName, a cell array for roi name 

        fprintf('Make ROI for (%s, %s)\n',sessid,runName);
   
        runDir = fullfile(fatDir,sessid,runName);
        
        % Path to the dt6 file
        dtFile = fullfile(runDir,'dti96trilin','dt6.mat');
        
        % Path to anatomy file
        if type == 1 %normal
            vAnatFile = fullfile(runDir,'t1',t1_name);
        elseif type == 2 %weird and needs gmwmi file
            vAnatFile = fullfile(runDir,'dti96trilin','mrtrix','dwi_aligned_trilin_noMEC_gmwmi.nii.gz');
        end
     
        % Path to fROI file        
        roiList = fullfile(runDir,'t1','ROIs',roiName);

        % Do convert
        dtiXformMrVistaVolROIs(dtFile, roiList, vAnatFile)        
    end
