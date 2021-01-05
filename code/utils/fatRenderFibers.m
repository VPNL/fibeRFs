function  fatRenderFibers(fatDir, sessid, runName, fgName,foi,t1name,hemi)
%  fatRenderFibers(fatDir, sessid, runName, fgName, hemi)
% fgName: full name of fg including path and postfix
% foi, a vector to indicate fiber of interest
% hemi, 'rh' or 'lh'
if nargin < 7, hemi = 'lh'; end

[~,fName] = fileparts(fgName);

if strcmp(hemi,'lh')
    cameraView = [-60,10];
    xplane =  [-1, 0, 0];
else strcmp(hemi,'rh')
    cameraView = [60,10];
    xplane =  [1,0,0];
end
zplane = [0, 0, -15];


% colorMap  = linspecer(length(foi));
% colorMap= jet(length(foi));
if foi==1
    if ~isempty(regexp(fgName, 'IOG'))
        colorMap = [215/255, 48/255, 39/255];
    elseif ~isempty(regexp(fgName, 'pFus'))
        colorMap = [244/255, 109/255, 67/255];
    elseif ~isempty(regexp(fgName, 'mFus'))
        colorMap = [253/255, 174/255, 97/255];
    elseif ~isempty(regexp(fgName, 'pSTS'))
        colorMap = [116/255, 173/255, 209/255];
    elseif ~isempty(regexp(fgName, 'mSTS'))
        colorMap = [69/255, 117/255, 180/255];
    elseif ~isempty(regexp(fgName, 'CoS'))
        colorMap = rgb('PaleGreen');
    else
        colorMap = rgb('Cyan');
    end
else

colorMap = [rgb('DarkOliveGreen');rgb('Orange'); rgb('Gold'); rgb('Pink');  ...
    [1 0.5 0];[1 1 0];[1 0 0];rgb('Purple'); rgb('Gold');...
    [0 0 1];rgb('mediumspringgreen');[0 1 1]; rgb('Magenta');rgb('Chocolate')]; 
%colorMap = [1 0 0; 1 1 0; 1 0 1; 0 0 1; 0 1 1; 0 1 0; 0 .4 1; 0 .6 .6; 0 .6 0; .4 0 0; 1 0.5 0; .5 .5 .5; .9 .9 0; 1 1 1; 0 0 0]; %red; yellow; magenta; blue; cyan; green  
    %colorMap = [1 0.5 0;1 1 0;1 0 1;1 0 0;0 0 1;0 1 1;0 1 0];
end

% set criteria
maxDist = 3;maxLen = 2;numNodes = 30;M = 'mean';maxIter = 1;count = false;
numfibers = 10000; %was 100
        fprintf('Plot fiber %s-%s:%s\n',sessid,runName,fgName);
        runDir = fullfile(fatDir,sessid,runName,'dti96trilin');
        afqDir = fullfile(runDir, 'fibers','afq');
        imgDir = fullfile(afqDir,'image');
        if ~exist(imgDir,'dir')
            mkdir(imgDir);
        end
        
        %% Load fg
        fgFile = fullfile(afqDir,fgName);
        if exist(fgFile,'file')
            load(fgFile);
            if exist('roifg')
            fg = roifg(foi);
            else
            fg = bothfg(foi); 
            colorMap = [0.9 0.75 0]; %[1 0 0;0 1 0; 0 1 1; 1 0 1];  
            end
                       
%             for i = 1:length(foi)
%                 fg(i) = AFQ_removeFiberOutliers(fg(i),maxDist,maxLen,numNodes,M,count,maxIter);
%             end
%             
         
            % b0 = readFileNifti(fullfile(runDir,'bin','b0.nii.gz'));
            b0 = readFileNifti(fullfile(fatDir,sessid,runName,'t1',t1name));
            fibers = extractfield(fg, 'fibers');
            I =  find(~cellfun(@isempty,fibers));
            
            if isempty(I)==0
            AFQ_RenderFibers(fg(I(1)),'numfibers',numfibers,...
                'color',colorMap(I(1),:),'camera',cameraView);
            
            for j = 2:length(I)
                AFQ_RenderFibers(fg(I(j)),'numfibers',numfibers,...
                    'color',colorMap(I(j),:),'newfig',false)
            end
            
            AFQ_AddImageTo3dPlot(b0,zplane);
            AFQ_AddImageTo3dPlot(b0,xplane);
            axis off
            axis square
            
            clear fg roifg
            
           print('-depsc',fullfile(imgDir,sprintf('%s_19.eps',fName)));
           print('-dtiff','-r300',fullfile(imgDir,sprintf('%s_19.tiff',fName)));
%             end
            %close all;
        end
end
