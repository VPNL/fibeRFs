function [ROIs] = standardROIs(which)
% choose all or some of our most-used ROIs (to save time typing)
% credit to Sonia
ROIset = {'V1' 'V2' 'V3' 'IOG_faces' 'pFus_faces' 'mFus_faces' 'pSTS_faces' 'mSTS_faces' 'CoS_places'};

if ~exist('which','var')
    n = 1:length(ROIset);
elseif containsTxt(which,'face')
    n = 4:8;
elseif containsTxt(which,'place')
    n = 9;
elseif containsTxt(which,'EVC')
    n = 1:3;
else n = which;
end

ROIs = {ROIset{n}};

end