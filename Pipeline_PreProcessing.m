clear all
clc
close all

% This is the pipeline that pre-processes the data collected by the AG
% Fanis from 2 Utah Arrays.

% Abhilash Dwarakanath, MPI for biological cybernetics, April 2017.

%% Specify directories

subjectName = 'Anton';

if strcmp(subjectName,'Anton')
    
    directories.recording = 'B:\A11\20200228';
    recDate = '28022020';
    subjID = 'Anton';
    
    cd(directories.recording) % Enter the data directory here % On Node 2
    if exist('PFC') ~=7 && exist('QNX') ~=7
        mkdir PFC
        mkdir QNX
    end
    
    directories.PFC = [directories.recording '\PFC'];
    directories.QNX = [directories.recording '\QNX'];
else
    directories.recording = 'B:\H07\20161105'; % Always follow this format post september '16
    recDate = '05112016';
    subjID = 'Hayo';
    cd(directories.recording) % Enter the data directory here % On Node 2
    if exist('PFC') ~=7 && exist('PPC') ~=7 && exist('QNX') ~=7
        mkdir PFC
        mkdir PPC
        mkdir QNX
    end
    
    directories.PFC = [directories.recording '\PFC'];
    directories.PPC = [directories.recording '\PPC'];
    directories.QNX = [directories.recording '\QNX'];
    
end

%% Get filenames for nev and ns6

% NEV files
nevfiles = dir( [directories.recording '/*.nev']); %get files matching pattern
nevfiles = {nevfiles.name};

% ns6 files

recfiles = dir( [directories.recording '/*.ns6']); %get files matching pattern
recfiles = {recfiles.name};

% ns2 files

emfiles = dir( [directories.recording '/*.ns2']); %get files matching pattern
emfiles = {emfiles.name};


%% Define parameters

global params % use global variable params so that PARFOR doesn't get cranky

% General params
params.fs = 30000; % in Hz
params.chans = 128;
params.elecs = 96;
params.nAnalogIns = 16;
params.clean = 0; % In samples. Usually for NSP1 this is the "Timestamp" value. If parsing a dataset which was recorded using MultiCentral (it would have an id i1 or i2) then set this to 106. Blackrock recording system had a bug which would pad the dataset with 105 zeros.

% For LFP stuff

params.targetFreq = 500;
%params.trialLength = 10; % in seconds!! This is either the length of your trial from the DGZ file or the pieces you use to compute the coherence
%params.tapers = 5/2;
params.decimCoeff = params.fs/params.targetFreq; % For reducing sampling rate
params.decimFacs = computeDecimationFactors(params.decimCoeff);

% For spike stuff

params.passband = [300 6000]/(params.fs/2);
params.low = 300; params.high = 6000;
%params.stopband = [100 6200]/(params.fs/2); % These are the parameters Alex Ecker uses
%params.passripple = 0.002; % Maximum allowed dB of the pass-band ripple
%params.stopatten = 60; % Minimum 50dB of attenuation in the stop-band.
params.minSpkNum = 4; % Currently this is set to number of spikes. We can later threshold it in Hz
params.numFeatures = 3; % THIS HAS TO BE LESS THAN OR EQUAL TO THE NUMBER OF SPIKES!!!!
params.stdmin = -5; % First number is for negative thresholded channels
params.stdmax = -50;                     % maximum threshold for detection
params.refractory = 0.5; %ms
params.beforeSpike = 0.5;
params.afterSpike = 1;
%params.min_spk = 1;
%params.max_spk = 10000;
params.minclus = 5;
params.maxclus = 5;
params.filterorder = 2; % We will use a 2nd order Butterworth filter.
params.offs = 2; %samples for artifact removal
params.thrsh = 60; % number of channels for artifact removal

% For PSTH stuff
%params.kerneltype = 'alpha'; %
%params.psthbinsize = 0.025; % in SECONDS please
%params.conditions = 8;

% for correlations stuff
%params.binsize = 0.015; % in SECONDS please
%params.maxlag = 0.5; % in SECONDS

%max_memo_GB = 8;
notchfilter = 0;

%% Parse the NEV files

% NEV files
nevfiles = dir( [directories.recording '/*.nev']); %get files matching pattern
nevfiles = {nevfiles.name};

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.recording)
    PFCevents = openNEV(['./' nevfiles{1}]);
    PPCevents = openNEV(['./' nevfiles{2}]);
    
else
    
    cd(directories.recording)
    PFCevents = openNEV(['./' nevfiles{1}]);
    
end

%% Get the events for a particular experiment - THIS IS NEEDED FOR THE EM FILE

% For both Bfsgrad1 and OriDisk1, no need to change epoch logical, because
% they are less than eight charachters. For others change epoch logical
% accordingly from the PFCevents.Data.comments.text

tagsNeeded = [{'Bfsgrad1'}, {'Bfsgrad1 end'}]; % INPUT THE REQUIRED TAG FROM THE ABOVE LIST HERE!! CHECK THE EXCEL LOGFILE BEFORE YOU INPUT THE TAGS
%tagsNeeded = [{'Bfsnatur1'}, {'Bfsnatur1 end'}];
%tagsNeeded = [{'Spontaneous1'}, {'Spontaneous1 end'}]; % INPUT THE REQUIRED TAG FROM THE ABOVE LIST HERE!! CHECK THE EXCEL LOGFILE BEFORE YOU INPUT THE TAGS
%tagsNeeded = [{'OriDiskfixoff1'}, {'OriDiskfixoff1 end'}];
%tagsNeeded = [{'AuSeq1_4'}, {'AuSeq1_4 end'}];
%tagsNeeded = [{'ImRivRE75_1'}, {'ImRivRE75_1 end'}];
%tagsNeeded = [{'BfsgradFsFixOff1'}, {'BfsgradFsFixOff1 end'}];

if strcmp(subjectName,'Hayo')==1
    
    % Create a specific task directory
    directories.taskdirPFC = [directories.PFC '\' tagsNeeded{1}];
    directories.taskdirPPC = [directories.PPC '\' tagsNeeded{1}];
    cd(directories.PFC)
    if ~exist(directories.taskdirPFC)
        mkdir(directories.taskdirPFC)
    end
    cd(directories.PPC)
    if ~exist(directories.taskdirPPC)
        mkdir(directories.taskdirPPC)
    end
else
    % Create a specific task directory
    directories.taskdirPFC = [directories.PFC '\' tagsNeeded{1}];
    
    cd(directories.PFC)
    if ~exist(directories.taskdirPFC)
        mkdir(directories.taskdirPFC)
    end
    
end

%Use ismember to get indices of the epochs that did happen in the experiment
expEpochTags = cellstr(PFCevents.Data.Comments.Text);

[ispresent,epochLogical] = ismember(tagsNeeded,expEpochTags);

if epochLogical(2)==0
    epochLogical(2) = epochLogical(1)+1;
end

% DUE TO A RECORDING BUG; WE WILL HAVE TO CHECK THIS AND CHANGE THE INDICES
% MANUALLY HERE. CONSULT THE LAB RECORD IF UNSURE!!!!!

epochinds = PFCevents.Data.Comments.TimeStamp(epochLogical);

%% Parse the data

if strcmp(subjectName,'Hayo')==1
    whichFiles = [1 2]; % Sometimes we did 2 recordings. If the main experiments are in the second recordings, change this to [3 4] for Hayo
    
    cd(directories.PFC)
    nc5exist=~isempty(dir('*.NC5')); 
    
    if nc5exist==0
        cd(directories.recording)
        [lengthPFC,lengthPPC]=write_data_NSx(params,recfiles,whichFiles,directories); % This parses the ns6 files
        write_data_NSx_single(params,emfiles,whichFiles,directories,nevfiles); % This parses the ns2 file. The EM is only recorded on the PFC NSP so we only parse 1
    else
        fprintf('This dataset has already been parsed!\n')
    end
    
else
    
    whichFiles = [1 2];
    
    cd(directories.PFC)
    nc5exist=~isempty(dir('*.NC5'));
    
    if nc5exist==0
        cd(directories.recording)
        [lengthPFC,nParts]=write_data_NSx_single(params,recfiles,whichFiles,directories,nevfiles);
        write_data_NSx_single(params,emfiles,whichFiles,directories,nevfiles);
    else
        fprintf('This dataset has already been parsed!\n')
    end
end

%% Plot and check the power-spectrum

% Just a quick printing of the spectra using a bart-hanning window
% periodogram. There is a bug in the check_lfp_power_NSX_parall.m code that
% doesn't really do the parallelisation inside. So I have parallelised it
% from the outside. This should work.

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.PFC)
    if exist('spectra','dir')==0
        parfor numChans = 1:params.chans
            check_lfp_power_NSX_parall(numChans,[],notchfilter,1)
            clf
        end
        cd(directories.recording)
        
        % Do PPC
        cd(directories.PPC)
        parfor numChans = 1:params.chans
            check_lfp_power_NSX_parall(numChans,[],notchfilter,1)
            clf
        end
    else
        fprintf('The spectra has been computed, moving on...\n')
    end
    cd(directories.recording)
    
else
    cd(directories.PFC)
    if exist('spectra','dir')==0
        parfor numChans = 1:params.chans
            check_lfp_power_NSX_parall(numChans,[],notchfilter,1)
            clf
        end
    else
        fprintf('The spectra has been computed, moving on...\n')
    end
    cd(directories.recording)
end

%% Recording lengths

if strcmp(subjectName,'Hayo')==1
    
    if params.clean==0
        startSamp = 1;
    else
        startSamp = params.clean+1;
    end
    
    if ~exist('lengthPFC')
        cd(directories.PFC)
        load('NSX_TimeStamps.mat');
        if exist('lts')
            lengthPFC = lts;
        else
            lengthPFC = lengthPFC;
        end
        cd(directories.recording);
    end
    
    if isstruct(lengthPFC)
        endSampPFC = lengthPFC.lengthPFC;
        decLengthPFC = ceil(lengthPFC.lengthPFC/params.decimCoeff);
    else
        endSampPFC = lengthPFC;
        decLengthPFC = ceil(lengthPFC/params.decimCoeff);
    end
    
    if ~exist('lengthPPC')
        cd(directories.PPC)
        load('NSX_TimeStamps.mat');
        if exist('lts')
            lengthPPC = lts;
        else
            lengthPPC = lengthPPC;
        end
        cd(directories.recording);
    end
    
    if isstruct(lengthPPC)
        endSampPPC = lengthPPC.lengthPPC;
        decLengthPPC = ceil(lengthPPC.lengthPPC/params.decimCoeff);
    else
        endSampPPC = lengthPPC;
        decLengthPPC = ceil(lengthPPC/params.decimCoeff);
    end
    
else
    
    if params.clean==0
        startSamp = 1;
    else
        startSamp = params.clean+1;
    end
    
    if ~exist('lengthPFC')
        cd(directories.PFC)
        load('NSX_TimeStamps.mat');
        if exist('lts')
            lengthPFC = lts;
        else
            lengthPFC = lengthPFC;
        end
        cd(directories.recording);
    end
    
    if isstruct(lengthPFC)
        endSampPFC = lengthPFC.lengthPFC;
        decLengthPFC = ceil(lengthPFC.lengthPFC/params.decimCoeff);
    else
        endSampPFC = lengthPFC;
        decLengthPFC = ceil(lengthPFC/params.decimCoeff);
    end
end

%% Spike detection

% Check lab notes for channels that have positive going spikes. I removed
% the plot continuous channels function because it is redundant. Its
% filters are different, its bands are different. All it does is that it
% prints out a 10s piece of the filtered signal along with the threshold
% line. If necessary, I can ask my getSpikes.m function to do it. -AD

% If a channel has positive spikes, it getSpikes.m will extract both
% positive and negative spikes because I assume that due to the nature of
% the Utah array, irrespective of a channel picking up positive spikes, it
% will as usual also pick up negative spikes. Its okay because if there are
% no negative spikes, then getSpikes will simply return the positive spikes
% without crashing out due to an error - AD

if strcmp(subjectName,'Hayo')==1
    
    if exist('PFCSpikes.mat')==0
        
        % Do PFC
        % Detection
        chanListPFC = -1.*ones(params.elecs,1);
        %for e.g. - chanListPFC(83) = 1; This will make it detect positive going
        %spikes
        poschans = find(chanListPFC==1);
        negchans = find(chanListPFC==-1);
        cd(directories.PFC)
        parfor chans = 1:params.elecs
            fprintf('Processing negative going spikes channel %d now...\n',chans)
            filenameraw = ['NSX' num2str(chans) '.NC5'];
            spikesPFC{chans} = getSpikes(filenameraw,params,startSamp,endSampPFC,chanListPFC(chans));
        end
        
        [spikesPFC] = removeCommonArtifacts(spikesPFC,params,4);
        save('PFCSpikes.mat','spikesPFC','-v7.3');
        
    else
        
        fprintf('Spike detection for PFC is done. Moving on...\n');
        
    end
    
    if exist('PPCSpikes.mat')==0
        % Do PPC
        %Spike detection
        chanListPPC = -1.*ones(params.elecs,1);
        poschans = find(chanListPPC==1);
        negchans = find(chanListPPC==-1);
        cd(directories.PPC)
        for chans = 1:params.elecs
            fprintf('Processing negative going spikes channel %d now...\n',chans)
            filenameraw = ['NSX' num2str(chans) '.NC5'];
            spikesPPC{chans} = getSpikes(filenameraw,params,startSamp,endSampPPC,chanListPPC(chans));
        end
        
        [spikesPPC] = removeCommonArtifacts(spikesPPC,params,4);
        save('PPCSpikes.mat','spikesPPC','-v7.3');
        
    else
        
        fprintf('Spike detection for PPC is done. Moving on...\n');
        
    end
    
else
    
    if exist('PFCSpikes.mat')==0
        
        % Do PFC
        % Detection
        chanListPFC = -1.*ones(params.elecs,1);
        %for e.g. - chanListPFC(83) = 1; This will make it detect positive going
        %spikes
        poschans = find(chanListPFC==1);
        negchans = find(chanListPFC==-1);
        cd(directories.PFC)
        for chans = 1:params.elecs
            fprintf('Processing negative going spikes channel %d now...\n',chans)
            filenameraw = ['NSX' num2str(chans) '.NC5'];
            spikesPFC{chans} = getSpikes(filenameraw,params,startSamp,endSampPFC,chanListPFC(chans));
        end
        
        [spikesPFC] = removeCommonArtifacts(spikesPFC,params,4);
        save('PFCSpikes.mat','spikesPFC','-v7.3');
        
    else
        
        fprintf('Spike detection for PFC is done. Moving on...\n');
        
    end
    
end

%% Do Klustakwik clustering

if strcmp(subjectName,'Hayo')==1
    
    if exist('autoResult.clu.1')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        cd(directories.PFC)
        
        parfor chans = 1:params.elecs
            clusteringSingleChannel(params,spikesPFC{chans},chans);
        end
        copyfile('L:\projects\AbhilashD\Preprocessing_Abhi_Toolbox\autoResult.xml','autoResult.xml');
        cd(directories.recording)
        
        if exist('spikesPPC') == 0
            
            cd(directories.PPC)
            load('PPCspikes.mat')
            
        end
        
        cd(directories.PPC)
        for chans = 1:129+params.elecs-1
            clusteringSingleChannel(params,spikesPPC{chans},chans);
        end
        
        copyfile('L:\projects\AbhilashD\Preprocessing_Abhi_Toolbox\autoResult.xml','autoResult.xml');
        
    else
        
        fprintf('Automatic clustering is already done, moving on...\n')
        
    end
    cd(directories.recording)
    
else
    if exist('autoResult.clu.1')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        cd(directories.PFC)
        
        parfor chans = 1:params.elecs
            clusteringSingleChannel(params,spikesPFCNegatives{chans},chans);
        end
        copyfile('L:\projects\AbhilashD\Preprocessing_Abhi_Toolbox\autoResult.xml','autoResult.xml');
    else
        
        fprintf('Automatic clustering is already done, moving on...\n')
        
    end
    cd(directories.recording)
end
