function [] = mesh_montage(in)

% Produce images of retinotopy prf model - phase and parameter maps
% (estimated one of a number of
% different ways) on a single-subject mesh
%
% JG 2016 - Changed for new map fit kid study.
% DF 2018 - changed for V1 myelination
% SP 2019 - changed for fixPRF project/sp-utils set-up
% DF 2019 - changed again for fibeRFs

%%%%%%%%%

setpref('mesh', 'layerMapMode', 'all');
setpref('mesh', 'overlayLayerMapMode', 'mean');

% go to the session directory
for i =1:length(in.sess)
     
    cd(fullfile(in.dataDir, in.sess{i}));    
    
    load mrSESSION.mat;
    %clear dataTYPES mrSESSION vANATOMYPATH;
    
    hG = initHiddenGray('GLMs',1);
    
    % load the ROI, get ready for the main loop
    
    % Delete from the list the ROIs the subject doesnt have
    toDelete = [];
    allMaps = in.roi;
    for m = 1:length(allMaps)
        if ~exist(fullfile('3DAnatomy','ROIs',[allMaps{m} '.mat']))
            if ~exist(fullfile('3Danatomy','ROIs',[allMaps{m} '.mat']))
                toDelete(end+1) = m;
                fprintf('Missing %s!\n',[allMaps{m} '.mat']);
            end
        end
    end
    allMaps(toDelete) = [];
    
    [hG , ok] = loadROI(hG, allMaps,[],[],[],0); % now loads shared maps, not local
    
    % Change colors if specified
    if ~isempty(in.colors)
        for y=1:numel(hG.ROIs)
            hG.ROIs(y).color = in.colors{y};
        end
    end
    
    % load the user prefs -- meshColorOverlay needs a 'ui' field :P
    %enforce consistent preprocessing / event-related parameters
    
    %hG.ui = load('Gray/userPrefs.mat');
    %hG.ui.dataTypeName= 'GLMs'; % in case ui preferences are different

    dtOpts = {'Averages' 'GLMs'}; % search space for relevant datatypes for our parameter maps
    for n = 1:length(dtOpts)
        if ~isempty(cellNum(dtOpts{n},{dataTYPES.name}))
        hG.ui.dataTypeName= dtOpts{n}; % in case ui preferences are different
        hG.curDataType = cellNum(dtOpts{n},{dataTYPES.name});
        end
    end
%     
    % Load and clip the parameter map
    switch in.whichMontage
        case 'ret'
            in.map = 'retModel-cssFit-fFit.mat';
            hG = rmSelect(hG,1,in.map);
            if strcmp(in.whatMap,'none')
            else hG = rmLoadDefault(hG); end;
        case 'face'
            % Load and clip the parameter map
            mapVariants = {'face-vs-all' 'Face1Face2' 'FacesVsAll' 'Faces-vs-allnonFace' 'Faces_vs_allnonFace' 'Faces_vs_all' 'faceadultfacechildVBodyLimbCarGuitarPlaceHouseWordNumber', 'facesVSall'};
            %mapVariants = {'faceadultfacechildVBodyLimbCarGuitarPlaceHouseWordNumber'};
            n=0; ok = 0;
            while ~ok && n < length(mapVariants)
            n=n+1;
            [hG,ok] = loadParameterMap(hG,[mapVariants{n} '.mat']);
            end
    end
    %hG=setDisplayMode(hG,'map'); hG=refreshScreen(hG);

    %%%%%%%%% which map do you want to load %%%%
%     cmapFile = '/sni-storage/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Retinotopy/code/Images/kidMap.mat';
%     load(cmapFile);
%     hG.ui.phMode.cmap = cmap;
    switch in.whatMap
        case 'ph'
            hG=setDisplayMode(hG,'ph');
            hG=refreshScreen(hG,1);
            hG.ui.displayMode='ph';
            hG.ui.phMode=setColormap(hG.ui.phMode, 'blueredyellowCmap');
            hG = cmapPolarAngleRGB(hG, 'both');        
        case 'amp'
            hG=setDisplayMode(hG,'amp');
            hG=refreshScreen(hG,1);
            hG.ui.displayMode='amp';
            cmapFile = '/sni-storage/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Retinotopy/code/Images/kidAmp.mat';
            load(cmapFile);
            front = repmat(cmap(1,:),272,1); cmap = vertcat(front,cmap);
            cmap = resample(cmap,16,33); cmap(cmap>1)=1;
            hG.ui.ampMode.cmap = cmap;
        case 'map'
            hG=setDisplayMode(hG,'map');
            hG=refreshScreen(hG,1);
            hG.ui.displayMode='map';
            %cmap =hsvCmap;
            hG.ui.mapMode=setColormap(hG.ui.mapMode,'autumnCmap');
            hG=setMapWindow(hG, [in.windowmin in.windowmax]);     
    end
    %hG=setDisplayMode(hG,'map'); hG=refreshScreen(hG);
    %% load up a mesh for the view
    % first, we set the hemisphere, by setting the cursor position's 3rd dim.
    if isequal(lower(in.hem), 'rh')
        hG.loc = [100 100 1]; else hG.loc = [100 100 200];end
    
    % go into anatomy folder
    if exist('3DAnatomy')
        cd('3DAnatomy')
    else
        cd('3Danatomy')
    end
    
    % load the mesh
    ok = 0; n = 0;
    while ~ ok 
        n= n+1;
    [hG,ok] = meshLoad(hG, [in.hem '_' in.meshNames{n}], 1);
    end
    pause(2);
    
    %cd(fullfile(in.dataDir,in.sess{i},'Retinotopy'));
    
    % set the view settings
    try
    [mesh1, settings]=meshRetrieveSettings(hG.mesh{1}, [in.hem '_' in.meshAngle]);
    catch
    [mesh1, settings]=meshRetrieveSettings(hG.mesh{1}, [in.hem '_Ventralzoom']);
    end
    % Recompute vertex
    vertexGrayMap = mrmMapVerticesToGray(...
        meshGet(mesh1, 'initialvertices'),...
        viewGet(hG, 'nodes'),...
        viewGet(hG, 'mmPerVox'),...
        viewGet(hG, 'edges'));
    hG.mesh{1} = meshSet(hG.mesh{1}, 'vertexgraymap', vertexGrayMap);
     
    %hG.ui.roiDrawMethod = 'filled perimeter'; hG = refreshScreen(hG);
    hG.ui.roiDrawMethod = 'boxes'; hG = refreshScreen(hG);

    % Set the mesh lighting to be the same across meshes
    meshLighting(hG, hG.mesh{1}, in.L, 1);
    meshColorOverlay(hG);
    meshUpdateAll(hG);
    pause(2);
    
    % grab the snapshot and crop it if it is a ventral view
    img_temp= mrmGet(hG.mesh{1}, 'screenshot' ) ./ 255;
    %image_crop=img_temp(100:512, 100:512, :);
    %img{i}=image_crop;
    img{i}=img_temp;
    hG = meshDelete(hG, Inf); %close all meshes
    plotNum = i;
    %subplot_tight(in.nrows,in.ncols,plotNum);
    subplot(in.nrows,in.ncols,plotNum);
    imagesc(img{i});
    set(gca,'CameraViewAngle',get(gca,'CameraViewAngle')-.1);
    set(gcf,'color','w');
    axis image;
    axis off;
    myTitle = title(in.sess{i},'fontsize',10); set(myTitle,'interpreter','none');
    set(myTitle);
end

% set colorbar
if in.colbar
    cbar = cbarCreate(cmap);
    cbar.cmap = cbar.cmap(128:224,:);
    cbar.nColors = 97;
    cbar.clim = [in.windowmin in.windowmax];
    subplot('Position', [0.05 0.05 .1 .04]);
    cbarDraw(cbar);
end

% save
checkDir(in.outDir);
cd(in.outDir);
set(gcf,'PaperPositionMode','auto')
%saveas(f, [in.imagefile(1:end - 4) '.fig'], 'fig');
print ('-dpng', '-r400', in.imagefile);
