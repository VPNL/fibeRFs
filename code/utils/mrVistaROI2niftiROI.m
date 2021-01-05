function mrVistaROI2niftiROI(session,roiList,t1name,t1folder,outPath)
%
% session path to session
% roiList cell array of roi names
% t1 is the fullfilename to the t1 anatomy file that the session is linked
% t1folder is the name of the folder where T1 anatomy is, (I added this because it changes from 3Danatomy to 3DAnatomy in the newer sessions -VN)
% to; the nifti will be generated based on this file
% outpath is the path to where you save your niftiROIs
%
% MB & KGS 1/15
%

cd(session);

hG = initHiddenGray(1,1,roiList);
t1name=fullfile(t1folder, t1name);


t1=readFileNifti(t1name); % we are loading the t1 nifti to get all the nifti structure information correct so when we later use mri_convert it will be aligned properly to orig.mgz
hG=loadAnat(hG); % this loads the anatomy based on the global variable vANATOMTPATH

anatSize=viewGet(hG,'anatomy size'); % get the anatomy size; 
% add a check that it matches the t1 size

if size(t1.data)~=anatSize
    fprintf('Error: size of t1 and vANAT do not match \n');
    return
end

if isempty(t1.data)
     fprintf('Error: t1 is empty\n');
    return
end

for i=1:length(hG.ROIs)
    %get ROI coords
    hG.selectedROI=i;
    coords = getCurROIcoords(hG);
    len = size(coords, 2);
    
    % make a 3D image with all points set to zero except ROI = roiColor
    roiData = zeros(anatSize);
    for ii = 1:len
        roiData(coords(1,ii), coords(2,ii), coords(3,ii)) = 1;
    end
    
    % Convert mrVista format to our preferred axial format for NIFTI
    mmPerVox = viewGet(hG, 'mmPerVox');
    [data, xform, ni] = mrLoadRet2nifti(roiData, mmPerVox);
    % create the nifti structure based on the original anatomy t1 structure
    % so all the transformation information matches up between the
    % functional ROIs and the t1 anatomy
    fROI=t1;
    fROI.data=ni.data; % data is the ROI data
    outname=fullfile(outPath, [hG.ROIs(i).name '.nii.gz']); % outfile for the nifti ROI
    fROI.fname=outname;
    fprintf('writing %s ...\n', outname);
    writeFileNifti(fROI); % write the nifti file
    
    
end


