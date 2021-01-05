% Find which sessions are missing certain views
% Helpful to run before the montage code

s1_setAllSessions
dataDir = ['/share/kalanit/Projects/fibeRFs/data/study1/loc/'];

% Right Hemisphere

rh_missing = []';

for i = 1:length(sessions)
    tf = [];
    absent = [];
    session = sessions{i};
    f = fullfile([dataDir, session, '/3DAnatomy/']);
    if exist(f)
        cd(f)
    else
        cd(fullfile([dataDir, session, '/3Danatomy/']))
    end
    load MeshSettings.mat
    for k = 1:numel(settings)
        tf(k) = strcmp(settings(k).name, 'rh_lat');
    end
    if any(tf)
        continue;
    else
        rh_missing{i} = session;
    end

end
% rh_missing_MZ(cellfun(@isempty,rh_missing_MZ)) = [];
rh_missing = rh_missing'

% Left Hemisphere


lh_missing = []';

for i = 1:length(sessions)
    tf = [];
    absent = [];
    session = sessions{i};
    f = fullfile([dataDir, session, '/3DAnatomy/']);
    if exist(f)
        cd(f)
    else
        cd(fullfile([dataDir, session, '/3Danatomy/']))
    end
    load MeshSettings.mat
    for k = 1:numel(settings)
        tf(k) = strcmp(settings(k).name, 'lh_lat');
    end
    if any(tf)
        continue;
    else
        lh_missing{i} = session;
    end

end
% lh_missing_MZ(cellfun(@isempty,lh_missing_MZ)) = [];
lh_missing = lh_missing'
