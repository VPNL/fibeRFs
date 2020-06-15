function c = roiColors(roiNames)
% now allow call to more than one color returned

if ~iscell(roiNames) roiNames = {roiNames}; end % if we're giving a single ROI as a string

rColors = {
    [0, 0, 0];...                   % V1 
    [0.4, 0.4, 0.4];...             % V2
    [0.6, 0.6, 0.6];...             % V3
    [215/255, 48/255, 39/255];...   % IOG - red
    [244/255, 109/255, 67/255];...  % pFus - orange
    [253/255, 174/255, 97/255];...  % mFus - yellow
    [116/255, 173/255, 209/255];... % pSTS - light blue
    [69/255, 117/255, 180/255];...  % mSTS - dark blue
    [0/255, 153/255, 0/255]         % CoS - green
    };

c = [];
allROIs = standardROIs;

for n = 1:length(roiNames)
    found = 0;
    m = cellNum(roiNames{n},allROIs);
    if ~isempty(m) % exact name match
        c = [c; rColors{m}]; found = 1;
    else
        for p = 1:length(allROIs)
            if ~found && containsTxt(roiNames{n},allROIs{p}) || containsTxt(allROIs{p},roiNames{n}) % partial name match
                c = [c; rColors{p}]; found = 1; end
        end
    end
    if ~found % no name match - assign grey
        fprintf('No color match for ROI %s...\n');
        c = [c; [.5 .5 .5]];
    end
end
