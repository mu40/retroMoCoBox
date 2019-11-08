function reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,varargin)
% function reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,varargin)
% 
% Dependencies:
%%%%%%%%%%%%%%%
%
%  %% Required %%:
%
%     SPM 12 (www.fil.ion.ucl.ac.uk/spm) 
%           For coregistration of FatNavs ('spm_realign').
%           This is obviously core to the concept of FatNavs for
%           motion-correction. I have tried using alternative registration
%           software, but SPM appears to have default parameters which work
%           well for FatNavs - being particular sensitive to sub-voxel
%           movements.
%
%     Michigan Image Reconstruction Toolbox (MIRT)
%     (http://web.eecs.umich.edu/~fessler/code/index.html)
%     from the group of Prof. J Fessler
%           This is used to perform the NUFFT 3D gridding operation used to deal
%           with the non-Cartesian k-space sampling after rotations have
%           been applied.
%
% Usage:
%%%%%%%%
%
% 'rawDataFile' - this is the filename for the raw data you exported from
%                 TWIX (including the full path). Due to the way I chose to
%                 name output files, if you have renamed the datafile it is
%                 required to still contain the text 'MID' followed by a
%                 number and then an underscore ('_') so that multiple
%                 files can be processed from the same directory and not
%                 get overwritten.
%
%  Optional arguments:-
%      'outRoot' - this gets used as the output folder for all processing
%                 (and the default root for the temporary folder which may
%                 generate several Gb of temporary files along the way - 
%                 - override this with 'tempRoot' option). The default 
%                 'outRoot' is to use the same folder as the raw data file.
%    
%      'tempRoot' - location for creating temporary files - can be many Gb.
%                   Default is to use the same location as 'outRoot'.
%
%      'bLinParSwap' - set this to '1' to indicate that the 'LIN/PAR swap'
%                     option was chosen in the MP2RAGE sequence. This
%                     alters which direction in k-space the FatNavs
%                     correspond to.
%
%      'bGRAPPAinRAM' - set this to '1' to perform the GRAPPA recon (if
%                       necessary!) for the host sequence entirely in RAM. This
%                       can be quite a lot faster - but requires sufficient
%                       RAM...! 
%                       (Note that the temporary variables for applying the
%                       MoCo are currently done outside of RAM as this
%                       is currently not the bottleneck for that part of
%                       the recon)
%
%      'bKeepGRAPPArecon' - set this to '1' to prevent deleting of the
%                           GRAPPA recon of the host data (which would be
%                           deleted by default). 
%
%      'bKeepReconInRAM' - set this to '1' to keep all the reconstructed
%                          files (INV1, INV2 and UNI, all with and without
%                          correction) in RAM, or default is to use a
%                          temporary file.
%
%      'bFullParforRecon' - set this to '1' to enable parfor over the NUFFT
%                           loop. For this to work, bKeepReconInRAM must
%                           also be set - and depending how many CPUs you
%                           have available on your matlabpool, you may need
%                           to have LOADS of RAM...  But it is much faster!
%       
%      'coilCombineMethod' - In MP2RAGE the default method is to match what
%                            I believe is done on the scanner - to weight
%                            each separately calculated UNI image by the
%                            square of the INV2 image for that coil.
%                            Slightly better results might be obtained for
%                            low SNR data by using a lower-res version of
%                            the INV2 image. Set this option to 'lowres' to
%                            try this. 
%                            For INV1 and INV2 the default is to combine
%                            the images using root sum-of-squares
%                            - this is not optimal, so the 'lowres' will
%                            also try to use a low-res version of INV2 to
%                            combine both. This corresponds to the coil
%                            combination method of Bydder et al, MRM 2002,
%                            47:539-548 - but as there is a smoothness
%                            parameter which needs tuning - and that it can
%                            also lead to signal voids, I haven't yet felt
%                            confident enough to make this the default
%                            processing.
%
%       'swapDims_xyz' - 3-component row-vector to determine whether to
%                        reverse each of x, y and z directions. Default is
%                        [0 0 1] which seems to work for a lot of parameter
%                        sets, but not all...  Check the orientation
%                        checker feature in the HTML!
%
%       'bZipNIFTIs' - Use '1' to apply gzip at the end to all the NIFTI
%                     files (default), otherwise just leave them uncompressed.
%
%       'bKeepPatientInfo' - Use '1' to keep sensitive info from raw data header  
%                           in HTML output (default). Use '0' to anomyize completely
%                           and use a string e.g. '0019' to insert the ID from
%                           another database.
%
%       'bKeepComplexImageData' - save out the complex data per coil as MATLAB files 
%                                 before and after application of MoCo.
%
%
%     
%   
% Matlab tools which are included (with 'assumed' permission, as I collected them online):
%
%     mapVBVD (from Philip Ehses) 
% *** Note that this code is presumably 'Siemens-sensitive' as it is    ***
% *** not freely available online, but only from the Siemens user forum.***
% *** Consequently the FatNavs recon code must also be considered       ***
% *** 'Siemens-sensitive' while it contains this code                   ***
% *** I have not included this part in the Github repository - please   ***
% *** email me if you would like it.                                    ***
%           For reading in the Siemens raw data format - and able to handle
%           very large datasets elegantly. 
%           I made small changes to the code to allow handling of the
%           FatNav data.
%
%     NIFTI Tools (from Jimmy Shen)
%           For reading in and saving out in the .nii NIFTI format. 
%           This code is provided here unaltered.            
%
%     GCC - coil compression (from Miki Lustig)
%           We work mostly with a 32-channel head coil, and so massive
%           speedups are possible when using coil compression. This code
%           from Miki Lustig's website implements the method described in
%           Zhang et. al MRM 2013;69(2):571-82.
%           For this code I have directly taken snippets and inserted them
%           into performHostGRAPPArecon.m
% 
%     export_fig (from Oliver Woodford)
%           Useful for making figures that you save from Matlab look nice!
%           :)
%
%     process_options (from Mark A. Paskin)
%           Useful for sending optional inputs to this master file.         
%
%     subplot1 (from Eran O. Ofek)
%           Nicer than Matlab's own way of doing subplots
%
% 
% Optional:
%     ImageMagick (www.imagemagick.org)
%           For making animated GIFs of results ('convert') in the html
%           report
%
%     FSL (www.fmrib.ox.ac.uk/fsl) - tested with v5.0
%           Used for brain extraction ('bet') in order to make MIP views 
%           of INV2 images (which show bright arteries) more interpretable
%
%
%
%
%
% To do:
%%%%%%%%
%
%  Coil compression for the FatNavs
%              - this is now implemented, but preliminary tests have shown
%                it doesn't work well... One problem is the weights
%                themselves for the sparse image - the other problem is
%                that acceleration of 16 is quite a lot...
%
%  Coil compression when no GRAPPA used
%              - this might make sense to still be able to speed up the
%                application of the retrospective motion-correction, but as
%                most people are probably scanning with GRAPPA in the host
%                sequence anyway, I haven't implemented this yet.
%
%
% -------------------------------------------------------------------------
% reconstructSiemensMP2RAGEwithFatNavs.m, 
%   v0.1 - daniel.gallichan@epfl.ch - June 2015
%   v0.2 -    -- February 2016 --
%        - added option to specify resolution of FatNavs 
%        - also added option to specify swapDims_xyz 
%        - added Patient info to HTML output
%        - added option to zip the NIFTIs at the end
%        - added option to keep FatNavs (and changed default behaviour to delete them)
%   v0.3 -   -- April 2016 -- trying to speed things up
%        - changed default oversampling from 2 down to 1.375 for NUFFT
%        - use parfor in NUFFT
%        - use parfor in SPM registration
%        - switch to using SPM12 (previously SPM8)
%        - now have option to do the MP2RAGE recon combination entirely in RAM
%   v0.4 -   -- August 2016
%        - changed default NUFFT oversampling from 1.375 to 1.5 (to reduce aliasing artifact)
%   v0.5 -   -- September 2016
%        - updated to latest version of Philipp Ehses' mapVBVD software (from 18/9/15)
%        - now renamed things to make it part of the 'RetroMoCoBox' and put
%          on Github
%        - Added the bFullParforRecon option to really speed things up if
%          you have enough CPUs and RAM available
%        - changed name of bKeepRecoInRAM to bKeepReconInRAM for consistency
%        - output files no longer start with 'a_host_'
%        - fit to versioning for the whole of 'RetroMoCoBox' as 0.5.0
%        
%   0.5.1 -  -- September 2016
%         - Fixed bug in handling of data acquired without GRAPPA
%         - Fixed bug in handling of data acquired with different orientations 
%
%   0.6.0 -  -- February 2017 - new contact email: gallichand@cardiff.ac.uk
%         - *CHANGED* handling of motion estimates - now average temporal
%           neighbours
%         - Add support for VD/VE data
%         - *RENAMED* output 'uniImage' and 'uniImage_corrected' to 'UNI'
%           and 'UNI_corrected' as seems to match INV1 and INV2 better
%         - Add animated GIF with zoom of front of brain where changes are
%           likely to be most noticeable and put in HTML
%
%   0.6.1 - -- May 2017
%         - Added Torben's option to anonymize data
%
%   0.6.2 -   -- July 2017
%         - Improved 'parfor' handling of data without second inversion
%           time
%         - Added manual discard of channels for HeadNeck_64 coil which
%           have a lot of signal in the neck
%
%   0.6.3 - -- August 2017
%         - Include the FatNav resolution in the HTML (and display to
%           screen)
%         - handle the case where 'PatientName' becomes 'tPatientName' for
%           no apparent reason
%
%   0.6.4 - -- Sep 2017
%         - Automatically set FatNavRes_mm based on field strength (7T - 2mm, 3T - 4mm)         
%
%   0.7.0 - Feb 2018
%         - Create sub-function reconstructSiemensVolume.m so that multiple
%           averages or repetitions can be handled directly in this code.
%           NB. Handling this properly would involve updated the sequence
%           code run on the scanner because it currently doesn't label the
%           FatNavs properly beyond the first volume of the host. This
%           feature probably won't be used enough to make that worthwhile
%           though...
%         - WARNING - I took this opportunity to cleanup the names of some
%           of the input options (always putting a 'b' for boolean in front
%           of logical options) so please check your calling code!
%
%
%   0.7.1 - -- Nov 2018
%         - Added new input flag 'bKeepComplexImageData' to allow saving
%           out of the complex data per coil as MATLAB files - before and
%           after application of MoCo.

reconPars.retroMocoBoxVersion = '0.7.1dev'; % put this into the HTML for reference
reconPars.rawDataFile = rawDataFile;

%% Check SPM and Fessler's toolbox are on path

if ~exist('spm.m','file') % could also check version, but that's more effort...
    disp('Error - SPM (ver 12) must be on the path')
    return
end

if ~exist('nufft_init.m','file')
    disp('Error - Fessler toolbox must be on the path for the NUFFT')
    return
end

%%

[reconPars.outRoot, reconPars.tempRoot, reconPars.bLinParSwap, reconPars.bGRAPPAinRAM, reconPars.bKeepGRAPPArecon, reconPars.bKeepReconInRAM, reconPars.bFullParforRecon,...
    reconPars.coilCombineMethod, reconPars.swapDims_xyz, reconPars.bZipNIFTIs,reconPars.bKeepPatientInfo,...
    reconPars.bKeepComplexImageData, reconPars.motMatFile, reconPars.doReverseCorr, reconPars.alignMotToCen, reconPars.doAverageMot,...
    reconPars.vNavDicomDir, reconPars.replaceReacqs] = process_options(varargin,...
    'outRoot',[],'tempRoot',[],'bLinParSwap',0,'bGRAPPAinRAM',0,'bKeepGRAPPArecon',0,'bKeepReconInRAM',0,...
    'bFullParforRecon',0,'coilCombineMethod','default','swapDims_xyz',[0 0 1],'bZipNIFTIs',1,'bKeepPatientInfo',1,...
    'bKeepComplexImageData',0,'motMatFile',[],'doReverseCorr',1,'alignMotToCen',0,'doAverageMot',0,'vNavDicomDir', 0, 'replaceReacqs', 0);


%%

if reconPars.bFullParforRecon && ~reconPars.bKeepReconInRAM
    disp('Error - you asked for the full parfor option (bFullParforRecon), but not to do the recon in RAM (bKeepReconInRAM)')
    return
end

%% Disply pie chart of current breakdown of reconstruction time

% % Here benchmarked on server allowing parpool size of 12 CPUs and 96 Gb
% % RAM, testing 600 um data
% 
% t_TotalTime = 24;
% t_ParseRaw = 86/60;
% t_reconFatNavs = 332/60;
% t_SPMrealign = 61/60;
% t_GRAPPArecon = 6;
% t_NUFFT = 8;
% t_other = t_TotalTime - t_ParsreRaw - t_reconFatNavs - t_SPMrealign - t_NUFFT;
% 
% figure(1001)
% set(gcf,'Position',[    88   415   588   506])
% clf
% pie([t_other t_ParseRaw t_reconFatNavs t_SPMrealign t_GRAPPArecon t_NUFFT],{'Other','Parse raw data file','Reconstruct FatNavs','SPM realign FatNavs','GRAPPA recon for host','NUFFT'})
% title({'Current breakdown of full reconstruction pipeline on 600 um data', ['total = ' num2str(t_TotalTime) ' mins, running on 12 CPUs with 96 Gb RAM'],char(datetime)})
% fontScale(1.4)
% export_fig('processingTimeBreakdown.png')


          

%%

startTime = clock;


%% Make raw data object (parses file, but does not load into RAM)

tic
twix_obj = mapVBVD_fatnavs(rawDataFile,'removeOS',1);
timingReport_parseRawDataFile = toc;

if length(twix_obj)>1 % on VE (and presumably VD as well) the raw data typically also has the adjcoilsens data as scan 1, so skip this
    twix_obj = twix_obj{2};
end

if ~isfield(twix_obj,'RTfeedback')
    disp('Error, no vNavs found in raw data file!')
    twix_obj % this displays the fields of the twix_obj that are present for comparison/debugging
    return
end

%% show the number of reacqs for each data type
f =fields(twix_obj); 
for ff=2:length(f)
    disp([f{ff} ' has ' num2str(twix_obj.(f{ff}).NReacq) ' reacqs'])
end
% in refscan NReacq is wrong it is counting non-refscan chunks that happen
% before the refscan 
% - look at sort(chunkIDs(firstReacqInd:end)) below
%% do a quick check that the meas.dat and vNav dicom were on the same day...
frameRef = split(twix_obj.hdr.Config.FrameOfReference,'.');
measDate = frameRef{11}(1:8);
dicomDirList = dir(reconPars.vNavDicomDir);
testDicom = dicomDirList(3).name;
testDicom_split = split(testDicom,'.');
disp('')
if contains(testDicom,measDate)    
    disp('')
    disp(['k-space data (' measDate ') and vNav DICOMs (' testDicom ') were acquired on the same day! Proceed!'])
elseif testDicom_split{11}(1:5)=='30000'
    disp(testDicom)
    disp(measDate)
    warning('looks like a retro-recon ... date check not performed')
else
    disp(testDicom)
    disp(measDate)
    error('k-space data and vNav DICOMs should be acquired on the same day...')
end
disp('')
%% call the python script to read the motion from the vNav DICOM headers
% script will save a file "rms_scores_vnav.h5" which will be read below
[datadir, ~, ~] = fileparts(reconPars.rawDataFile);
if reconPars.vNavDicomDir==0
    error('need the path to the directory containing the vNav DICOM files')
else
    str = sprintf('python /autofs/space/vault_020/users/rfrost/moco/parse_vNav_Motion/vnav/parse_vNav_Motion.py --input %s/* --tr 2.5  --radius 64 --rms- --output-mats-path %s',...
        reconPars.vNavDicomDir,datadir);
    disp(str)
    [stat,res] = unix(str);
    disp(res)
    if stat~=0
        error('reading vNav motion params failed!')
    end
end
%% use the vNav motion to decide which data to replace with reacqs
% get the scores output from parse_vNav_Motion.py
rms_scores = h5read([datadir '/rms_scores_vnav.h5'],'/rmsscores');
rms_scores = [-1 rms_scores']; % first chunk does not have a score...
%% find lines to replace with reacquisitions
% have to analyse reacqs in twix_obj.imageWithRefscan because the chunkIDs
% correspond to these in the acquisition with integrated GRAPPA refscan
% need to know which reacqs to replace in:
% - twix_obj.refscan
% - twix_obj.image

if twix_obj.imageWithRefscan.NReacq > 1
    % find the indices of when the chunk counter changes
    new_counter_inds = [1 find(diff(twix_obj.imageWithRefscan.iceParam(8,:)))+1];
    
    % extract the list of chunk IDs and chunk Counters:
    chunkIDs=twix_obj.imageWithRefscan.iceParam(7,new_counter_inds);
    chunkCounters=twix_obj.imageWithRefscan.iceParam(8,new_counter_inds);
    % look at which Lines were acquired - this should match the data found in
    % the reacqs
    chunkLines = twix_obj.imageWithRefscan.Lin(new_counter_inds); % e.g. 2,4,6 ..
    chunkLines = chunkLines - 1; % seems that a -1 is necessary..? data must be placed
    % by mapVBVD in the index preceding these Lin numbers ?
    
    lastScanBeforeReacq = find(chunkIDs==max(chunkIDs));
    firstReacqInd = lastScanBeforeReacq + 2; % the chunkID=0 is always first reacq, so need to take the one after
    
    % intialize a list of scores to be updated with reacq scores as we loop
    % through
    rms_scores_current = rms_scores(1:lastScanBeforeReacq);
    reacqs_to_use_for_replacement = -1*ones(1,twix_obj.imageWithRefscan.NLin);
    
    reacq = 0; clear reacq_lines
    for ind = firstReacqInd:length(chunkIDs)
        reacq = reacq + 1;
        
        ID = chunkIDs(ind);
        LINE = chunkLines(ind);
        SCOREreacq = rms_scores(ind);
        SCOREorig = rms_scores_current(ID + 1); % list position is ID+1, IDs start from 0
        
        fprintf('reacq %d: chunk ID = %d, reacq score = %.2f, original score %.2f\n'...
            ,reacq,ID,SCOREreacq,SCOREorig)
        
        thisLine = twix_obj.imageWithRefscan(1,1,:,1,1,1,1,1,1,1,1,1,1,1,1,1,reacq+2);
        line = find(abs(squeeze(thisLine)));
        
        if line~=LINE, error('different ways of finding lines should be the same'), end
        reacq_lines(reacq) = line;
        
        if SCOREreacq<SCOREorig
            fprintf('replace line %d with twix_obj reacq index = %d\n\n',line,reacq+2)
            rms_scores_current(ID + 1) = SCOREreacq;
            % set the line position in replacement array with reacq number from
            % twix_obj.imageWithRefscan         %
            % *** important that for replacement imageWithRefscan is used ***
            % the same line was reacquired at least once, the best one will be
            % stored
            reacqs_to_use_for_replacement(line) = reacq+2;
        end
        
    end
    lines_to_replace = find(reacqs_to_use_for_replacement~=-1);
    reacq_to_use = reacqs_to_use_for_replacement(lines_to_replace);
    
    % find the lines in image and in refscan
    new_counter_inds = [1 find(diff(twix_obj.image.Lin))+1];
    LINESimage = twix_obj.image.Lin(new_counter_inds) - 1; % need a -1 to match the lines found in reacq data
    
    new_counter_inds = [1 find(diff(twix_obj.refscan.Lin))+1];
    LINESrefscan = twix_obj.refscan.Lin(new_counter_inds) - 1;
    
    new_counter_inds = [1 find(diff(twix_obj.imageWithRefscan.Lin))+1];
    LINESimageWithRefscan = twix_obj.imageWithRefscan.Lin(new_counter_inds) - 1;
    
    % split the lines into image and refscan (for GRAPPA recon and calib parts)
    IND = ismember(lines_to_replace,LINESimage);
    lines_IMAGE = lines_to_replace(IND);
    reacq_to_use_IMAGE = reacq_to_use(IND);
    
    disp('image lines to replace and reacqs to use:')
    disp([lines_IMAGE;reacq_to_use_IMAGE])
    
    % refscan lines to replace and reacqs to use
    IND = ismember(lines_to_replace,LINESrefscan);
    lines_REFSCAN = lines_to_replace(IND);
    reacq_to_use_REFSCAN = reacq_to_use(IND);
    
    disp('refscan lines to replace and reacqs to use:')
    disp([lines_REFSCAN;reacq_to_use_REFSCAN])
    
    reconPars.RData_lines_to_replace    = lines_to_replace;
    reconPars.RData_reacq_to_use        = reacq_to_use;
    reconPars.RData_lines_IMAGE         = lines_IMAGE;
    reconPars.RData_reacq_to_use_IMAGE  = reacq_to_use_IMAGE;
    reconPars.RData_lines_REFSCAN       = lines_REFSCAN;
    reconPars.RData_reacq_to_use_REFSCAN= reacq_to_use_REFSCAN;
    
end

%% Run the reconstruction on each volume in the raw data

nAve = twix_obj.image.NAve;
nRep = twix_obj.image.NRep;

if nAve > 1

    disp('Multiple Averages detected in raw data file')
    
    for iAve = 1:nAve % can't use parfor here as then it wouldn't work inside each thread
        
        disp(['Processing average ' num2str(iAve)])
        
        thisReconPars = reconPars;
        thisReconPars.iAve = iAve;
        timingReport{iAve} = reconstructSiemensVolume_mh(twix_obj,thisReconPars);
    end
    
else
    if nRep > 1
        disp('Error - handling of data with multiple repetitions not yet implemented...! Please contact gallichand@cardiff.ac.uk for more info')
        return
    else
        timingReport = reconstructSiemensVolume_mh(twix_obj,reconPars);       
    end
end


%%
stopTime = clock;
totalTime = etime(stopTime,startTime)/60/60;
totalTime_hrs = floor(totalTime);
if totalTime_hrs > 0
    totalTime_mins = round(rem(totalTime,totalTime_hrs)*60);
else
    totalTime_mins = round(totalTime*60);
end

fprintf('*************************************************************\n')
fprintf('***** reconstructSiemensMP2RAGEwithFatNavs.m completed! *****\n')
fprintf('*************************************************************\n')
fprintf(['Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins\n']);

%%


    


end

