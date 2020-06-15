%% transforms ROIs from niftis to fs labels
% Edited to remove adding hemisphere prefixes
% DF 11/2019

function nii2label(subject,dataPath,roiList,roiVal,labelPath,hemi)

	for roi = length(roiList)
	
		
		if isempty(dataPath)
			dataPath = fullfile('/biac2/kgs/3Danat/', subject,'/niftiRois/');
		end
		if isempty(labelPath)
			labelPath=fullfile('/biac2/kgs/3Danat/FreesurferSegmentations/',subject,'/label/');
		end
		if isempty(roiVal)
			roiVal(roi) = 1;
		end

		reslicePath = fullfile('/biac2/kgs/3Danat/FreesurferSegmentations/', subject ,'/mri/orig.mgz');
		%inName=strcat(hemi,'_',roiList(roi),'.nii.gz');
        inName=strcat(roiList(roi),'.nii.gz');
		inFile=fullfile(dataPath,inName);
		%outName=strcat(hemi,'_',roiList(roi),'_conformed.nii.gz');
        outName=strcat(roiList(roi),'_conformed.nii.gz');
		outFile=fullfile(dataPath,outName);
		
		cmd = ['mri_convert -rl ' reslicePath ' -rt nearest -ns 1 --conform ' inFile{1} ' ' outFile{1}];
		unix(cmd);
		

		%labelName = strcat(hemi,'.',roiList(roi),'.label');
        labelName = strcat(roiList(roi),'.label');
		saveLabel=fullfile(labelPath,labelName);
		cmd = ['mri_vol2label --c ' outFile{1} ' --id ' num2str(roiVal(roi)) ' --l ' saveLabel{1}];
		unix(cmd);
		
	end

