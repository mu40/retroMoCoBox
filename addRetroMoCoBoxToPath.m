% add retroMoCoBox subfolders to the path (without including all the git
% subfolders!)

retroMoCoPath = which('addRetroMoCoBoxToPath.m');
[retroMoCoPath, ~] = fileparts(retroMoCoPath);

addpath(retroMoCoPath);
addpath([retroMoCoPath '/export_fig']);
addpath([retroMoCoPath '/fatnavtools']);
addpath([retroMoCoPath '/vnavtools']);
addpath([retroMoCoPath '/generaltools']);
addpath([retroMoCoPath '/images']);
addpath([retroMoCoPath '/niftitools']);
addpath([retroMoCoPath '/mapVBVD_20160905']);
