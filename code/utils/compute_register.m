%% Compute register.dat file manually if it hasn't been created

fsDirs   = {'sp201803_v6','MJH25','MBA24_scn022715_recon0319_v6'};

% make register.dat
for s=1:length(fsDirs)
    % subject = the subject name of the freesurfer directory
    subject = fsDirs{s};
    origPath = ['/biac2/kgs/3Danat/FreesurferSegmentations/' subject '/mri/orig.mgz'];
    outFile = ['/biac2/kgs/3Danat/FreesurferSegmentations/' subject '/surf/register.dat']; %originally went to subject/label/register.dat
    cmd = ['tkregister2 --mov ' origPath ' --noedit --s ' subject ' --regheader --reg ' outFile];
    unix(cmd)
end
