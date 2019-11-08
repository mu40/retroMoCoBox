run('~/l/sand/retro-moco/retroMoCoBox/addRetroMoCoBoxToPath.m')
addpath('~/l/git/spm12')
run('~/l/sand/retro-moco/mirt/setup.m')

% rawDataFile = '~/l/sand/retro-moco/data/meas_MID123_mp2rage_FN1000b_FatNav_1mm.dat'; % Example.
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190816/meas_MID01868_FID23788_T1w_MPR_vNav_4e.dat'; % HCP

% PMC off
rawDataFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/meas_MID00026_FID37142_ABCD_T1w_MPR_vNavPMCoff_jawOpen.dat';
% motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/jaw_004_pace.mat';
% motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/jaw_004_flirt_bd4.mat';
motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/jaw_004_flirt_bd2_up2.mat';

% PMC on
% rawDataFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/meas_MID00027_FID37143_ABCD_T1w_MPR_vNavPMCon_jawOpen.dat';
% motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/jaw_007_pace.mat';


doReverseCorr = 0; % Zero means correct retrospectively.
% alignMotToCen = 0; % Default alignment to first frame had lower entropy.



% Fastest option: if you have loads of RAM (tested with 12 CPUs and 96 Gb of RAM on both 1 mm and 600 um data)

reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,...
    'motMatFile',motMatFile,...
    'doReverseCorr',doReverseCorr,...
    'alignMotToCen',0,... % Alignment to first frame had lower entropy.
    'doAverageMot',0,...
    'bGRAPPAinRAM',1,...
    'bKeepReconInRAM',1,...
    'bFullParforRecon',1);

%% if you have plenty of RAM:

% reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'FatNavRes_mm',FatNavRes_mm,'bGRAPPAinRAM',1);


%% otherwise do the slower version:

% reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'FatNavRes_mm',FatNavRes_mm,'bGRAPPAinRAM',0);


%% ROB

alignToStartOrEnd=0; % 1=to start; 2=to end
% Set the resolution of the FatNavs that were used (default is 2 for 7T and 4 for 3T)
vNavRes_mm = 8;
FatNavRes_mm = 2;

rawDataFile = '/autofs/space/sand_001/users/mu40/retro-moco/data/meas_MID00027_FID37143_ABCD_T1w_MPR_vNavPMCon_jawOpen.dat';
vNavDicomDir = '/autofs/space/sand_001/users/mu40/retro-moco/data/jaw_dcm/007_ABCD_T1w_MPR_vNav_setter';
reconstructSiemensMP2RAGEwithvNavs(rawDataFile,'vNavDicomDir',vNavDicomDir,'vNavRes_mm',vNavRes_mm,'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',1,'bKeepFatNavs',1,'alignToStartOrEnd',alignToStartOrEnd);

%% HCP subject with high z motion

run('~/l/sand/retro-moco/retroMoCoBox/addRetroMoCoBoxToPath.m')
addpath('~/l/git/spm12')
run('~/l/sand/retro-moco/mirt/setup.m')
hcp_dir = '/autofs/space/sand_001/users/mu40/retro-moco/data_hcp';

rawDataFile = fullfile(hcp_dir, 'meas_MID01868_FID23788_T1w_MPR_vNav_4e.dat'); % running
vNavDicomDir = fullfile(hcp_dir, 'dcm_MID01868');
% motMatFile = fullfile(hcp_dir, 'mot_MID01868_pace.mat');
% motMatFile = fullfile(hcp_dir, 'mot_MID01868_imf1.mat');
motMatFile = fullfile(hcp_dir, 'mot_MID01868_imf005.mat');
% motMatFile = fullfile(hcp_dir, 'mot_MID01868_imf12.mat');

% rawDataFile = fullfile(hcp_dir, 'meas_MID01876_FID23796_T1w_MPR_vNav_4e.dat');
% vNavDicomDir = fullfile(hcp_dir, 'dcm_MID01876');
% motMatFile = fullfile(hcp_dir, 'mot_MID01876_pace.mat');
% motMatFile = fullfile(hcp_dir, 'mot_MID01876_imf1.mat');
% motMatFile = fullfile(hcp_dir, 'mot_MID01876_imf005.mat');
% motMatFile = fullfile(hcp_dir, 'mot_MID01876_imf12.mat');


doReverseCorr = 1; % Zero means correct retrospectively.
% alignMotToCen = 0; % Default alignment to first frame had lower entropy.
reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,...
    'motMatFile',motMatFile,...
    'doReverseCorr',doReverseCorr,...
    'alignMotToCen',0,...
    'doAverageMot',0,...
    'vNavDicomDir', vNavDicomDir,...
    'replaceReacqs', 0,...
    'bGRAPPAinRAM',1,...
    'bKeepReconInRAM',1,...
    'bFullParforRecon',1);

%% Jaw motion at randomized times with Rob Jones

run('~/l/sand/retro-moco/retroMoCoBox/addRetroMoCoBoxToPath.m')
addpath('~/l/git/spm12')
run('~/l/sand/retro-moco/mirt/setup.m')

% % Correction applied to: vNavs.
% rawDataFile = '/autofs/space/sand_001/users/mu40/retro-moco/data_jaw2/meas_MID02617_FID44558_ABCD_T1w_MPR_vNavPMCoff_jawOpen.dat';
% motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data_jaw2/mot_MID02617_pace.mat';

% Correction applied to: none.
rawDataFile = '/autofs/space/sand_001/users/mu40/retro-moco/data_jaw2/meas_MID02618_FID44559_ABCD_T1w_MPR_vNavPMCoff_neither_jawOpen.dat';
vNavDicomDir = '/autofs/space/sand_001/users/mu40/retro-moco/data_jaw2/dcm_MID02618';
% motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data_jawll2/mot_MID02618_pace.mat';
% motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data_jaw2/mot_MID02618_flirt_bd2_up4/4d_fsl.mat';
motMatFile = '/autofs/space/sand_001/users/mu40/retro-moco/data_jaw2/mot_MIeD02618_flirt_up4/4d_fsl.mat';


doReverseCorr = 0; % Zero means correct retrospectively.
reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,...
    'motMatFile',motMatFile,...
    'doReverseCorr',doReverseCorr,...
    'alignMotToCen',0,... % Alignment to first frame had lower entropy.
    'doAverageMot',0,...
    'vNavDicomDir', vNavDicomDir,...
    'replaceReacqs', 0,...
    'bGRAPPAinRAM',1,...
    'bKeepReconInRAM',1,...
    'bFullParforRecon',1);
