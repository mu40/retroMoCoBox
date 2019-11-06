% Obviously this needs to point to your data (Siemens MP(2)RAGE with 3D-FatNavs)
% rawDataFile = '/home/gallicha/temp/autoTransferTemp/tempRaw/meas_MID36_mp2rage_FN600b_FatNav_06mm.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20180427_DanG_exampleMPRAGE/recon20190826-TCLbox/meas_MID123_mp2rage_FN1000b_FatNav_1mm.dat';
rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20180427_DanG_exampleMPRAGE/recon20190915-TCLbox/meas_MID123_mp2rage_FN1000b_FatNav_1mm.dat';

D_hcp = '/autofs/space/nihilus_002/users/HCP/subjects/';
D_v20 = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/';
D_v21 = '/autofs/space/vault_021/users/rfrost/hcp-memprage/';
D_ger = '/cluster/gerenuk/user/rfrost/moco/vnav/hcp-memprage/';

% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190827/mid1404-t1w-mpr/meas_MID01404_FID26123_T1w_MPR_vNav_4e.dat';
rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190826/mid1268-t1w-mpr/meas_MID01268_FID25987_T1w_MPR_vNav_4e.dat';

rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190816/mid1868/meas_MID01868_FID23788_T1w_MPR_vNav_4e.dat';
vNavDicomDir = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190816/016_T1w_MPR_vNav_4e_T1w_setter';
vNavDicomDir = [D_hcp 'HCA6603969/HCA6603969_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190816/mid/meas_MID01876_FID23796_T1w_MPR_vNav_4e.dat';
% vNavDicomDir = '/autofs/space/vault_020/users/rfrost/moco/hcp-memprage/20190816/028_T1w_MPR_vNav_4e_T1w_setter';
% vNavDicomDir = [D_hcp 'HCA6603969/HCA6603969_V1_A/DICOM_ScanDirs/028_T1w_MPR_vNav_4e_T1w_setter/'];

% rawDataFile = [D_v20 '20191004/mid1627/meas_MID01627_FID36031_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA8931693/HCA8931693_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

% rawDataFile = [D_v21 '20191015/mid592/meas_MID00592_FID39239_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA9317175/HCA9317175_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

% good example:
% rawDataFile = [D_v21 '20191015/mid599/meas_MID00599_FID39246_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA9317175/HCA9317175_V1_A/DICOM_ScanDirs/026_T1w_MPR_vNav_4e_T1w_setter/'];
% above recon in progress...done

% good example:
% rawDataFile = [D_v20 '20190724/mid507/meas_MID00507_FID17485_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA9052973/HCA9052973_V1_B/DICOM_ScanDirs/059_T1w_MPR_vNav_4e_T1w_setter/'];
% above recon in progress... done


% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/undo-vnavPMC/20190823_evaluate_vNav_undoPMC/human_repeatedMove/mid1632-vnav-on-mov/meas_MID01632_FID116635_ABCD_T1w_MPR_vNavON_mov.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/undo-vnavPMC/20190823_evaluate_vNav_undoPMC/human_repeatedMove/mid1633-vnav-off-mov/meas_MID01633_FID116636_ABCD_T1w_MPR_vNavOFF_mov.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/undo-vnavPMC/20190823_evaluate_vNav_undoPMC/human_repeatedMove/mid1634-vnav-on-still/meas_MID01634_FID116637_ABCD_T1w_MPR_vNavON_still.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/undo-vnavPMC/20190823_evaluate_vNav_undoPMC/pha_moveAfterRef/mid1616-vnav-on-movAfterRef/meas_MID01616_FID116619_ABCD_T1w_MPR_vNav_movAfterRef.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/undo-vnavPMC/20190823_evaluate_vNav_undoPMC/pha_moveAfterRef/mid1617-vnav-on-pos2/meas_MID01617_FID116620_ABCD_T1w_MPR_vNav_pos2.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/undo-vnavPMC/20190823_evaluate_vNav_undoPMC/pha_moveAfterRef/mid1615-vnav-on-pos1/meas_MID01615_FID116618_ABCD_T1w_MPR_vNav_pos1.dat';

% 20191014: testing retro moco on scans with vnav PMC and TCL to compare
% the motion
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/abcd_vnav_trac_physio/20180427_abcd_physio_03/mid570-vnav';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/abcd_vnav_trac_physio/20180427_abcd_physio_03/mid579-vnav/meas_MID00579_FID10329_mpragetrac_mov_ABCD_match_vNav.dat';
% rawDataFile = '/autofs/space/vault_020/users/rfrost/moco/abcd_vnav_trac_physio/20180427_abcd_physio_03/mid580-vnav/meas_MID00580';
% TCLdir = '/autofs/space/vault_020/users/rfrost/moco/abcd_vnav_trac_physio/20180427_abcd_physio_03/2018-04-27_abcd_physio_03/2018-04-27_exp';
%%
% And wherever you installed the Retro-MoCo-Box
run('/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/code/mocoWorkshopISMRM17/retroMoCoBoxDemo_CapeTown_Full/retroMoCoBoxTCL-vNav/addRetroMoCoBoxToPath.m')
% And the TCL Suite needs to be on the path too
addpath /autofs/cluster/gerenuk/user/rfrost/moco/tracoline/code/TracSuiteMatlabFiles/
% addpath /autofs/space/vault_020/users/rfrost/moco/parse_vNav_Motion/vnav/
alignToStartOrEnd=0; % 1=to start; 2=to end

% And SPM 12
addpath('../spm12')

% And the Michigan Image Reconstruction Toolbox (MIRT) (http://web.eecs.umich.edu/~fessler/code/) for the NUFFT
run('../mirt/setup.m')

% Set the resolution of the FatNavs that were used (default is 2 for 7T and 4 for 3T)
vNavRes_mm = 8;
FatNavRes_mm = 2;

%% Fastest option: if you have loads of RAM (tested with 12 CPUs and 96 Gb of RAM on both 1 mm and 600 um data)

% rawDataFile = [D_v20 '20190906/mid115/meas_MID00115_FID28436_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA8402969/HCA8402969_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];
% 
% rawDataFile = [D_v20 '20190927/mid886/meas_MID00886_FID34249_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA6629280/HCA6629280_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];
% 
% rawDataFile = [D_v20 '20190930/mid340/meas_MID00340_FID34755_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA9482089/HCA9482089_V1_A/DICOM_ScanDirs/022_T1w_MPR_vNav_4e_T1w_setter/'];
% 
% rawDataFile = [D_v20 '20191008/mid121/meas_MID00121_FID37072_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA9586506/HCA9586506_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];
% 
% rawDataFile = [D_v20 '20190830/mid2128/meas_MID02128_FID26843_T1w_MPR_vNav_4e.dat'];
% vNavDicomDir = [D_hcp 'HCA7501259/HCA7501259_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

for SCAN=1

if SCAN==1
rawDataFile = [D_v20 '20190709/mid202/meas_MID00202_FID13695_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA9761801/HCA9761801_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

end

reconstructSiemensMP2RAGEwithvNavs(rawDataFile,'vNavDicomDir',vNavDicomDir,'vNavRes_mm',vNavRes_mm,'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',1,'bKeepFatNavs',1,'alignToStartOrEnd',alignToStartOrEnd);

end
%%

% for SCAN=[15 2:14]
for SCAN=[16]

if SCAN==1
rawDataFile = [D_v20 '20190710/mid624/meas_MID00624_FID14098_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA8846400/HCA8846400_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==2
rawDataFile = [D_v20 '20190715/mid124/meas_MID00124_FID15003_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA9512173/HCA9512173_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==3
% good example:
rawDataFile = [D_v20 '20190725/mid927/meas_MID00927_FID17903_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA6789909/HCA6789909_V1_B/DICOM_ScanDirs/041_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==4
rawDataFile = [D_v20 '20190730/mid2036/meas_MID02036_FID19006_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA8953805/HCA8953805_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==5
rawDataFile = [D_v20 '20190731/mid2711/meas_MID02711_FID19681_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA7154769/HCA7154769_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==6
rawDataFile = [D_v20 '20190802/mid3176/meas_MID03176_FID20146_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA9578406/HCA9578406_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==7
rawDataFile = [D_v20 '20190806/mid3985/meas_MID03985_FID20952_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA8720983/HCA8720983_V2_A/DICOM_ScanDirs/021_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==8
rawDataFile = [D_v20 '20190809/mid4792/meas_MID04792_FID21754_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA8917700/HCA8917700_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==9
rawDataFile = [D_v20 '20190814/mid615/meas_MID00615_FID22559_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA7310353/HCA7310353_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==10
rawDataFile = [D_v20 '20190819/mid60/meas_MID00060_FID24158_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA7719084/HCA7719084_V2_A/DICOM_ScanDirs/024_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==11
% NOT SURE ABOUT THIS ONE *********** check the nii.gz MPR
rawDataFile = [D_v20 '20190823/mid713/meas_MID00713_FID25434_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA6427066/HCA6427066_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==12
rawDataFile = [D_v20 '20190913/mid790/meas_MID00790_FID30639_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA9743900/HCA9743900_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==13
rawDataFile = [D_v21 '20191018a/mid66/meas_MID00066_FID40320_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA8364482/HCA8364482_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==14
rawDataFile = [D_v21 '20191018b/mid58/meas_MID00058_FID40406_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_hcp 'HCA9636494/HCA9636494_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/'];

elseif SCAN==15
% good example:
rawDataFile = [D_ger '20191025/mid1014/meas_MID01014_FID42954_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_ger '20191025/separated-DICOMS/im017-T1w_setter/017/'];

elseif SCAN==16
% good example:
rawDataFile = [D_ger '20191026/mid1316/meas_MID01316_FID43260_T1w_MPR_vNav_4e.dat'];
vNavDicomDir = [D_ger '20191026/separated-DICOMS/im016-T1w_setter/016/'];

end

reconstructSiemensMP2RAGEwithvNavs(rawDataFile,'vNavDicomDir',vNavDicomDir,'vNavRes_mm',vNavRes_mm,'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',1,'bKeepFatNavs',1,'alignToStartOrEnd',alignToStartOrEnd);

end
% 20190709 HCA9761801/HCA9761801_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190710 HCA8846400/HCA8846400_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190715 HCA9512173/HCA9512173_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190725 HCA6789909/HCA6789909_V1_B/DICOM_ScanDirs/041_T1w_MPR_vNav_4e_T1w_setter/
% 20190730 HCA8953805/HCA8953805_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190731 HCA7154769/HCA7154769_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190802 HCA9578406/HCA9578406_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190806 HCA8720983/HCA8720983_V2_A/DICOM_ScanDirs/021_T1w_MPR_vNav_4e_T1w_setter/
% 20190809 HCA8917700/HCA8917700_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190814 HCA7310353/HCA7310353_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190819 HCA7719084/HCA7719084_V2_A/DICOM_ScanDirs/024_T1w_MPR_vNav_4e_T1w_setter/
% 20190823 ?????retroRecon? HCA6427066/HCA6427066_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20190913 HCA9743900/HCA9743900_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20191018a HCA8364482/HCA8364482_V1_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/
% 20191018b HCA9636494/HCA9636494_V2_A/DICOM_ScanDirs/016_T1w_MPR_vNav_4e_T1w_setter/



%%
TCLtimeOffset_ms = 0; % -82229 is from manual logging that day, 900 is observed additional offset in other data
reconstructSiemensMP2RAGEwithTCL(rawDataFile,TCLdir,'TCLtimeOffset_ms',TCLtimeOffset_ms,'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',1);


%% if you have plenty of RAM:

% reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'FatNavRes_mm',FatNavRes_mm,'bGRAPPAinRAM',1);


%% otherwise do the slower version:

reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'FatNavRes_mm',FatNavRes_mm,'bGRAPPAinRAM',0);
