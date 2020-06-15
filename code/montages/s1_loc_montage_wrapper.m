clear all; close all;

s1_setAllSessions
in.sess = sessions;

in.whichMontage = 'face';
in.montName = 'face';
in.hem ='lh';

% inputs for this experiment
faceROIs = standardROIs('face');
fROIPre = 'fiberRFsclean_f_'; 

maps={};
for r = 1:length(faceROIs) %face ROIs
    maps = horzcat(maps,{[fROIPre in.hem '_' faceROIs{r}]});
end

in.roi = maps;
in.expt = 'study1';
in.imagefile = [in.expt '_' in.hem '_' in.montName '_N' num2str(length(in.sess)) '.png'];

%%%%%%%%%
% some defaults

in.dataDir = ['/share/kalanit/Projects/fibeRFs/data/' in.expt '/loc'];
in.outDir = [in.dataDir '/montages'];
in.L.ambient = [.5 .5 .5];
in.L.diffuse = [.3 .3 .3];

in.meshAngle='lat';%
in.meshNames={'inflated_200_1.mat' 'smooth_200_1.mat' 'smoothed_200_1.mat'};
in.clip=1;

in.nrows=3;
in.ncols=7;
in.colbar=0; %0 = no color bar

if strcmp(in.montName,'face') 
    %baseDir = 'faceLoc'; 
    in.whatMap = 'map'; 
    %in.colors = {'y' 'm' 'r' 'w' 'c'}; 
    in.colors = {'k' 'k' 'k' 'k' 'k'}; 
    in.windowmin= 3;
    in.windowmax = 15;
else baseDir = 'Retinotopy'; 
    in.whatMap='ph'; % ph = phase, ec = eccentricity, amp = pRF size; map = parameter map'
    in.colors = {'y' 'k' 'w' 'k' 'w' 'k' 'w' 'k', 'w'}; 
    in.windowmin= 0.5;
    in.windowmax = 80;
end

f = figure('Position', [50, 50, 2700, 1300]);
%f = niceFig([.5 .5 .6 .6]);
pause(1);

mesh_montage(in)
